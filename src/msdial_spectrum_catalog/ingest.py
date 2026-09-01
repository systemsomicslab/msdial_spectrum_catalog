from __future__ import annotations

import hashlib
import json
import zlib
from dataclasses import dataclass, field
from pathlib import Path

from .identifiers import make_id
from .parsers import parse_number, read_analysis_csv, read_msp, read_mztab, read_tsv, sha256_file
from .storage import initialize, transaction


@dataclass
class IngestReport:
    run_id: str
    samples: int = 0
    features: int = 0
    spectra: int = 0
    alignments: int = 0
    alignment_members: int = 0
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def _artifact_type(path: Path) -> str:
    lower = path.name.lower()
    if lower.endswith(".mdprovenance.tsv"):
        return "alignment_provenance"
    if lower.endswith(".mdalign"):
        return "alignment_matrix"
    if lower.endswith(".mdpeak"):
        return "sample_feature_table"
    if lower.endswith(".mdmsp"):
        return "msp_spectra"
    if lower.endswith(".mztab"):
        return "mztab_m"
    if lower.endswith(".txt"):
        return "parameter_or_text"
    if lower.endswith(".csv"):
        return "analysis_metadata"
    return "other"


def _put_blob(connection, payload: dict) -> str:
    raw = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode("utf-8")
    payload_sha = hashlib.sha256(raw).hexdigest()
    connection.execute(
        "INSERT OR IGNORE INTO spectrum_blob(payload_sha256, compression, uncompressed_bytes, payload) VALUES (?, 'zlib-json', ?, ?)",
        (payload_sha, len(raw), zlib.compress(raw, level=9)),
    )
    return payload_sha


def ingest_run(
    database: str | Path,
    run_directory: str | Path,
    repository: str,
    accession: str,
    analysis_unit: str,
    *,
    study_title: str | None = None,
    separation_type: str | None = None,
    ion_mode: str | None = None,
    acquisition_type: str | None = None,
    msdial_version: str | None = None,
    interactive_version: str | None = None,
    analysis_files_csv: str | Path | None = None,
    parameter_file: str | Path | None = None,
) -> IngestReport:
    initialize(database)
    run_dir = Path(run_directory).resolve()
    files = sorted(path for path in run_dir.iterdir() if path.is_file())
    relevant = [path for path in files if _artifact_type(path) != "other"]
    for optional_path in (analysis_files_csv, parameter_file):
        if optional_path:
            resolved = Path(optional_path).resolve()
            if resolved not in relevant:
                relevant.append(resolved)
    relevant.sort()
    artifact_labels = {
        path: path.name if path.parent == run_dir else f"inputs/{path.name}"
        for path in relevant
    }
    hashes = {path: sha256_file(path) for path in relevant}
    fingerprint_source = "\n".join(f"{artifact_labels[path]}\t{hashes[path]}" for path in relevant)
    fingerprint = hashlib.sha256(fingerprint_source.encode("utf-8")).hexdigest()
    study_id = make_id("study", repository, accession)
    unit_id = make_id("unit", repository, accession, analysis_unit)
    run_id = make_id("run", repository, accession, analysis_unit, fingerprint[:20])
    report = IngestReport(run_id=run_id)

    with transaction(database) as connection:
        connection.execute(
            "INSERT OR IGNORE INTO study(study_id, repository, accession, title) VALUES (?, ?, ?, ?)",
            (study_id, repository, accession, study_title),
        )
        connection.execute(
            "INSERT OR IGNORE INTO analysis_unit(analysis_unit_id, study_id, external_unit_id, separation_type, ion_mode, acquisition_type) VALUES (?, ?, ?, ?, ?, ?)",
            (unit_id, study_id, analysis_unit, separation_type, ion_mode, acquisition_type),
        )
        connection.execute(
            "INSERT OR IGNORE INTO analysis_run(run_id, analysis_unit_id, run_fingerprint, output_directory, msdial_version, interactive_version) VALUES (?, ?, ?, ?, ?, ?)",
            (run_id, unit_id, fingerprint, str(run_dir), msdial_version, interactive_version),
        )

        artifact_ids: dict[Path, int] = {}
        for path in relevant:
            connection.execute(
                "INSERT OR IGNORE INTO artifact(run_id, artifact_type, relative_path, source_path, sha256, byte_size) VALUES (?, ?, ?, ?, ?, ?)",
                (run_id, _artifact_type(path), artifact_labels[path], str(path), hashes[path], path.stat().st_size),
            )
            artifact_ids[path] = connection.execute(
                "SELECT artifact_id FROM artifact WHERE run_id = ? AND relative_path = ?", (run_id, artifact_labels[path])
            ).fetchone()[0]

        alignment_stems = {path.stem for path in relevant if path.suffix.lower() == ".mdalign"}
        samples: dict[str, str] = {}
        features: dict[tuple[str, int], str] = {}
        analysis_metadata = {}
        metadata_path = Path(analysis_files_csv).resolve() if analysis_files_csv else None
        if metadata_path is None:
            metadata_path = next((path for path in relevant if path.suffix.lower() == ".csv" and "analysis" in path.name.lower()), None)
        if metadata_path is not None:
            analysis_metadata = read_analysis_csv(metadata_path)
        for peak_path in (path for path in relevant if path.suffix.lower() == ".mdpeak"):
            sample_name = peak_path.stem
            sample_id = make_id("sample", repository, accession, analysis_unit, sample_name)
            samples[sample_name] = sample_id
            sample_metadata = analysis_metadata.get(sample_name, {})
            raw_path = sample_metadata.get("file_path") or None
            raw_name = raw_path.replace("\\", "/").rstrip("/").rsplit("/", 1)[-1] if raw_path else sample_name
            connection.execute(
                "INSERT OR IGNORE INTO sample(sample_id, analysis_unit_id, sample_name, raw_file_name, raw_file_path) VALUES (?, ?, ?, ?, ?)",
                (sample_id, unit_id, sample_name, raw_name, raw_path),
            )
            for row_number, row in read_tsv(peak_path):
                peak_id = parse_number(row.get("Peak ID"), int)
                if peak_id is None:
                    continue
                feature_id = make_id("feature", run_id, sample_name, peak_id)
                features[(sample_name, peak_id)] = feature_id
                connection.execute(
                    """INSERT OR IGNORE INTO feature(
                        feature_id, run_id, sample_id, master_peak_id, local_peak_id, ms1_scan_index,
                        rt_min, precursor_mz, height, area, name, adduct, source_artifact_id, source_row
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        feature_id, run_id, sample_id, peak_id, peak_id, parse_number(row.get("Scan"), int),
                        parse_number(row.get("RT (min)")), parse_number(row.get("Precursor m/z")),
                        parse_number(row.get("Height")), parse_number(row.get("Area")), row.get("Name"),
                        row.get("Adduct"), artifact_ids[peak_path], row_number,
                    ),
                )

        sample_msp_paths = [
            path for path in relevant
            if path.suffix.lower() == ".mdmsp" and path.stem not in alignment_stems
        ]
        for msp_path in sample_msp_paths:
            sample_name = msp_path.stem
            sample_id = samples.get(sample_name)
            if sample_id is None:
                report.warnings.append(f"No mdpeak table was found for {msp_path.name}")
                continue
            for record in read_msp(msp_path):
                peak_id = parse_number(record.comment_tokens.get("PEAKID"), int)
                if peak_id is None:
                    report.errors.append(f"Missing PEAKID in {msp_path.name} record {record.index}")
                    continue
                feature_id = features.get((sample_name, peak_id))
                if feature_id is None:
                    report.errors.append(f"{msp_path.name} PEAKID={peak_id} has no mdpeak row")
                payload_sha = _put_blob(connection, {"fields": record.fields, "peaks": record.peaks})
                spectrum_id = make_id("spectrum", run_id, sample_name, "deconvoluted", peak_id)
                connection.execute(
                    """INSERT OR IGNORE INTO spectrum(
                        spectrum_id, run_id, sample_id, feature_id, spectrum_kind, source_peak_id,
                        ms1_scan_index, ms2_scan_index, precursor_mz, rt_min, ion_mode, peak_count,
                        payload_sha256, source_artifact_id, source_record
                    ) VALUES (?, ?, ?, ?, 'deconvoluted', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        spectrum_id, run_id, sample_id, feature_id, peak_id,
                        parse_number(record.comment_tokens.get("MS1SCAN"), int),
                        parse_number(record.comment_tokens.get("MS2SCAN"), int),
                        parse_number(record.fields.get("PRECURSORMZ")),
                        parse_number(record.fields.get("RETENTIONTIME")), record.fields.get("IONMODE"),
                        len(record.peaks), payload_sha, artifact_ids[msp_path], record.index,
                    ),
                )

        alignments: dict[int, str] = {}
        for alignment_path in (path for path in relevant if path.suffix.lower() == ".mdalign"):
            for row_number, row in read_tsv(alignment_path, "Alignment ID\t"):
                alignment_id = parse_number(row.get("Alignment ID"), int)
                if alignment_id is None:
                    continue
                alignment_feature_id = make_id("alignment", run_id, alignment_id)
                alignments[alignment_id] = alignment_feature_id
                connection.execute(
                    """INSERT OR IGNORE INTO alignment_feature(
                        alignment_feature_id, run_id, alignment_master_id, alignment_local_id,
                        average_rt_min, average_mz, name, source_artifact_id, source_row
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        alignment_feature_id, run_id, alignment_id, alignment_id,
                        parse_number(row.get("Average Rt(min)")), parse_number(row.get("Average Mz")),
                        row.get("Metabolite name"), artifact_ids[alignment_path], row_number,
                    ),
                )

        alignment_msp_paths = [
            path for path in relevant
            if path.suffix.lower() == ".mdmsp" and path.stem in alignment_stems
        ]
        for msp_path in alignment_msp_paths:
            for record in read_msp(msp_path):
                alignment_id = parse_number(record.comment_tokens.get("PEAKID"), int)
                alignment_feature_id = alignments.get(alignment_id) if alignment_id is not None else None
                if alignment_feature_id is None:
                    report.errors.append(f"{msp_path.name} record {record.index} has no matching alignment row")
                    continue
                payload_sha = _put_blob(connection, {"fields": record.fields, "peaks": record.peaks})
                spectrum_id = make_id("spectrum", run_id, "consensus", alignment_id)
                connection.execute(
                    """INSERT OR IGNORE INTO spectrum(
                        spectrum_id, run_id, alignment_feature_id, spectrum_kind, source_peak_id,
                        precursor_mz, rt_min, ion_mode, peak_count, payload_sha256,
                        source_artifact_id, source_record
                    ) VALUES (?, ?, ?, 'alignment_consensus', ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        spectrum_id, run_id, alignment_feature_id, alignment_id,
                        parse_number(record.fields.get("PRECURSORMZ")),
                        parse_number(record.fields.get("RETENTIONTIME")), record.fields.get("IONMODE"),
                        len(record.peaks), payload_sha, artifact_ids[msp_path], record.index,
                    ),
                )

        provenance_paths = [path for path in relevant if path.name.lower().endswith(".mdprovenance.tsv")]
        if alignments and not provenance_paths:
            report.errors.append("Alignment output exists but .mdprovenance.tsv is missing")
        for provenance_path in provenance_paths:
            for row_number, row in read_tsv(provenance_path):
                alignment_id = parse_number(row.get("alignment_master_id"), int)
                file_id = parse_number(row.get("file_id"), int)
                sample_name = row.get("file_name", "")
                alignment_feature_id = alignments.get(alignment_id) if alignment_id is not None else None
                sample_id = samples.get(sample_name)
                if alignment_feature_id is None:
                    report.errors.append(f"Unknown alignment ID {alignment_id} in {provenance_path.name}")
                    continue
                if sample_id is None:
                    report.errors.append(f"Unknown sample {sample_name!r} in {provenance_path.name}")
                    continue
                source_peak_id = parse_number(row.get("source_master_peak_id"), int)
                feature_id = features.get((sample_name, source_peak_id)) if source_peak_id is not None and source_peak_id >= 0 else None
                has_source_peak = row.get("has_source_peak", "").lower() == "true"
                if has_source_peak and feature_id is None:
                    report.errors.append(f"Alignment {alignment_id}, {sample_name} references missing peak {source_peak_id}")
                member_id = make_id("alignment-member", run_id, alignment_id, sample_name)
                connection.execute(
                    """INSERT OR IGNORE INTO alignment_member(
                        alignment_member_id, alignment_feature_id, sample_id, feature_id, file_id,
                        is_representative, has_source_peak, source_master_peak_id, source_local_peak_id,
                        ms1_scan_index, ms2_scan_index, rt_min, mz, height, source_artifact_id, source_row
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        member_id, alignment_feature_id, sample_id, feature_id, file_id,
                        int(row.get("is_representative", "").lower() == "true"), int(has_source_peak),
                        source_peak_id, parse_number(row.get("source_peak_id"), int),
                        parse_number(row.get("ms1_raw_spectrum_id_top"), int),
                        parse_number(row.get("ms2_raw_spectrum_id"), int),
                        parse_number(row.get("rt_min")), parse_number(row.get("mz")),
                        parse_number(row.get("height")), artifact_ids[provenance_path], row_number,
                    ),
                )
                if row.get("is_representative", "").lower() == "true":
                    connection.execute(
                        "UPDATE alignment_feature SET representative_sample_id = ?, representative_feature_id = ? WHERE alignment_feature_id = ?",
                        (sample_id, feature_id, alignment_feature_id),
                    )

        for mztab_path in (path for path in relevant if path.name.lower().endswith(".mztab")):
            smf_to_alignment: dict[str, str] = {}
            pending: list[tuple[int, str, dict[str, str]]] = []
            for row_number, section, row in read_mztab(mztab_path):
                pending.append((row_number, section, row))
                if section == "SMF":
                    record_id = row.get("SMF_ID", "")
                    alignment_id = parse_number(record_id, int)
                    if alignment_id in alignments:
                        smf_to_alignment[record_id] = alignments[alignment_id]
            for row_number, section, row in pending:
                id_field = f"{section}_ID"
                record_id = row.get(id_field, "")
                refs = row.get("SMF_ID_REFS", "") if section == "SML" else row.get("SME_ID_REFS", "")
                alignment_feature_id = smf_to_alignment.get(record_id) if section == "SMF" else None
                if section == "SML" and refs:
                    first_ref = refs.replace("|", ",").split(",", 1)[0].strip()
                    alignment_feature_id = smf_to_alignment.get(first_ref)
                record_uid = make_id("mztab", run_id, section, record_id)
                connection.execute(
                    """INSERT OR IGNORE INTO mztab_record(
                        mztab_record_id, run_id, section, record_id, parent_refs_json,
                        alignment_feature_id, record_json, source_artifact_id, source_row
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        record_uid, run_id, section, record_id,
                        json.dumps(refs.split("|") if refs else []), alignment_feature_id,
                        json.dumps(row, ensure_ascii=False, separators=(",", ":")),
                        artifact_ids[mztab_path], row_number,
                    ),
                )

        counts = {
            "samples": ("SELECT COUNT(*) FROM sample WHERE analysis_unit_id = ?", unit_id),
            "features": ("SELECT COUNT(*) FROM feature WHERE run_id = ?", run_id),
            "spectra": ("SELECT COUNT(*) FROM spectrum WHERE run_id = ?", run_id),
            "alignments": ("SELECT COUNT(*) FROM alignment_feature WHERE run_id = ?", run_id),
            "alignment_members": (
                "SELECT COUNT(*) FROM alignment_member WHERE alignment_feature_id IN (SELECT alignment_feature_id FROM alignment_feature WHERE run_id = ?)",
                run_id,
            ),
        }
        for field_name, (query, parameter) in counts.items():
            setattr(report, field_name, connection.execute(query, (parameter,)).fetchone()[0])
    return report
