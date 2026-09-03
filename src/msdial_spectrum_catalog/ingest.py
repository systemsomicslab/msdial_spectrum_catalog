from __future__ import annotations

import hashlib
import json
import zlib
from dataclasses import dataclass, field
from pathlib import Path

from .identifiers import make_id
from .parsers import (
    fold_headers,
    folded_value,
    parse_number,
    read_analysis_csv,
    read_msp,
    read_mztab,
    read_tsv,
    sha256_file,
    split_refs,
)
from .storage import initialize, transaction


@dataclass
class IngestReport:
    run_id: str
    samples: int = 0
    features: int = 0
    spectra: int = 0
    alignments: int = 0
    alignment_members: int = 0
    msdial_annotation_results: int = 0
    mztab_sme_linked: int = 0
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


# MS-DIAL writes this into the metabolite name field of rows that carry no annotation claim at all.
NOT_ANNOTATED_NAME_PREFIXES = ("unknown",)

# MsdialCore/Utility/DataAccess.cs SetMoleculeMsPropertyAsSuggested prefixes a suggested annotation with
# "no MS2: " when MS2RawSpectrumID < 0 and with "low score: " otherwise. "w/o MS2: " is the MS-DIAL 4
# spelling, commented out in MS-DIAL 5, but MS-DIAL 4 exports remain ingestable.
# Colon-terminated rather than space-terminated: MS-DIAL 5 writes "no MS2: " with a space, while MS-DIAL 4
# code searches for "w/o MS2:" without one, so both spellings occur in exported names.
MSDIAL_SUGGESTION_PREFIXES = (
    ("no ms2:", "precursor_only"),
    ("w/o ms2:", "precursor_only"),
    ("low score:", "low_score"),
)

# Reference records whose NAME is an in-house identifier rather than a compound name.
UNNAMED_REFERENCE_PREFIXES = ("riken",)


def is_annotated_metabolite_name(name: str | None) -> bool:
    """Report whether an MS-DIAL metabolite name field asserts an annotation."""
    if name is None:
        return False
    text = name.strip()
    if not text or text.lower() == "null":
        return False
    lowered = text.lower()
    return not any(lowered.startswith(prefix) for prefix in NOT_ANNOTATED_NAME_PREFIXES)


def classify_metabolite_name(name: str | None) -> tuple[str, str | None, bool]:
    """Split an MS-DIAL metabolite name into annotation_kind, candidate_name and candidate_is_named.

    A "no MS2: " row is a precursor-only suggestion: no product-ion spectrum was acquired, so it can
    never support spectral-library evidence. A "low score: " row does have a product-ion spectrum, which
    failed the search criteria. Collapsing either into the same bucket as a real MS/MS match would
    overstate the evidence, so the distinction is stored rather than inferred downstream.
    """
    text = (name or "").strip()
    lowered = text.lower()
    kind = "msms_matched"
    for prefix, suggestion_kind in MSDIAL_SUGGESTION_PREFIXES:
        if lowered.startswith(prefix):
            kind = suggestion_kind
            text = text[len(prefix):].strip()
            break
    candidate = text or None
    is_named = bool(candidate) and not candidate.lower().startswith(UNNAMED_REFERENCE_PREFIXES)
    return kind, candidate, is_named


_ANNOTATION_TEXT_COLUMNS = (
    ("metabolite_name", ("Metabolite name", "Name")),
    ("formula", ("Formula",)),
    ("ontology", ("Ontology",)),
    ("inchikey", ("INCHIKEY", "InChIKey")),
    ("smiles", ("SMILES",)),
    ("adduct", ("Adduct type", "Adduct")),
    ("annotation_tag", ("Annotation tag (VS1.0)", "Annotation tag")),
    ("comment", ("Comment",)),
)
_ANNOTATION_FLAG_COLUMNS = (
    ("is_rt_matched", ("RT matched",)),
    ("is_mz_matched", ("m/z matched",)),
    ("is_msms_matched", ("MS/MS matched",)),
)
_ANNOTATION_SCORE_COLUMNS = (
    ("rt_similarity", ("RT similarity",)),
    ("mz_similarity", ("m/z similarity",)),
    ("ccs_similarity", ("CCS similarity",)),
    ("simple_dot_product", ("Simple dot product",)),
    ("weighted_dot_product", ("Weighted dot product",)),
    ("reverse_dot_product", ("Reverse dot product",)),
    ("matched_peaks_count", ("Matched peaks count",)),
    ("matched_peaks_percentage", ("Matched peaks percentage",)),
    ("total_score", ("Total score",)),
)
_DOT_PRODUCT_FIELDS = ("simple_dot_product", "weighted_dot_product", "reverse_dot_product")

# MS-DIAL writes -1 into these two columns to mean "not applicable", never as a measurement. A count and a
# percentage cannot be negative, so -1 is normalized away rather than stored as if it were a value.
_NOT_APPLICABLE_NEGATIVE_FIELDS = ("matched_peaks_count", "matched_peaks_percentage")

_ANNOTATION_RESULT_FIELDS = (
    tuple(name for name, _ in _ANNOTATION_TEXT_COLUMNS)
    + tuple(name for name, _ in _ANNOTATION_FLAG_COLUMNS)
    + tuple(name for name, _ in _ANNOTATION_SCORE_COLUMNS)
    + ("score_convention", "annotation_kind", "candidate_name", "candidate_is_named")
)


def _annotation_block(row: dict[str, str]) -> dict[str, object] | None:
    """Collect MS-DIAL's annotation columns from one .mdpeak or .mdalign row."""
    folded = fold_headers(row)
    block: dict[str, object] = {
        name: folded_value(folded, *headers) for name, headers in _ANNOTATION_TEXT_COLUMNS
    }
    if not is_annotated_metabolite_name(block["metabolite_name"]):
        return None
    for name, headers in _ANNOTATION_FLAG_COLUMNS:
        raw = folded_value(folded, *headers)
        block[name] = None if raw is None else int(raw.lower() in {"true", "1", "yes"})
    for name, headers in _ANNOTATION_SCORE_COLUMNS:
        block[name] = parse_number(folded_value(folded, *headers))
    for name in _NOT_APPLICABLE_NEGATIVE_FIELDS:
        value = block[name]
        if value is not None and value < 0:
            block[name] = None
    kind, candidate, is_named = classify_metabolite_name(block["metabolite_name"])
    block["annotation_kind"] = kind
    block["candidate_name"] = candidate
    block["candidate_is_named"] = int(is_named)
    if kind == "precursor_only":
        # A precursor-only suggestion had no product-ion spectrum to compare (MS2RawSpectrumID < 0), so
        # no dot product was ever computed. .mdalign writes null here, but .mdpeak writes 0.000 -- the
        # unset default of MsScanMatchResult's float fields. Left as 0.0 it reads as "compared, scored
        # zero", which is a different and stronger statement than "never compared". A non-zero value is
        # kept, so a future MS-DIAL change shows up as an anomaly instead of being silently discarded.
        for name in _DOT_PRODUCT_FIELDS:
            if block[name] == 0.0:
                block[name] = None
    # score_convention describes the three dot-product columns only. The exported columns are plain
    # cosines: MsScanMatching's GetSimpleDotProduct/GetWeightedDotProduct return cos-squared and are stored
    # in MsScanMatchResult.Squared*, but those fields are used only for threshold comparison, and both text
    # exporters write the non-squared computed properties (IMetadataAccessor.cs:114-116,
    # IAnalysisMetadataAccessor.cs:123-125). Total score is an unnormalized weighted composite on its own
    # scale -- 2.611 for a real match in the reference demo -- so it is never a cosine.
    has_dot_product = any(block[name] is not None for name in _DOT_PRODUCT_FIELDS)
    block["score_convention"] = "cosine" if has_dot_product else None
    return block


def _insert_annotation_result(
    connection,
    run_id: str,
    subject_kind: str,
    subject_id: str,
    block: dict[str, object],
    artifact_id: int,
    row_number: int,
) -> None:
    result_id = make_id("msdial-annotation", run_id, subject_kind, subject_id, 1)
    columns = ", ".join(_ANNOTATION_RESULT_FIELDS)
    placeholders = ", ".join("?" for _ in _ANNOTATION_RESULT_FIELDS)
    connection.execute(
        f"""INSERT OR IGNORE INTO msdial_annotation_result(
            msdial_annotation_result_id, run_id, subject_kind, subject_id, feature_id,
            alignment_feature_id, rank, annotator_id, {columns}, source_artifact_id, source_row
        ) VALUES (?, ?, ?, ?, ?, ?, 1, '', {placeholders}, ?, ?)""",
        (
            result_id, run_id, subject_kind, subject_id,
            subject_id if subject_kind == "feature" else None,
            subject_id if subject_kind == "alignment_feature" else None,
            *(block[name] for name in _ANNOTATION_RESULT_FIELDS),
            artifact_id, row_number,
        ),
    )


def _artifact_type(path: Path) -> str:
    lower = path.name.lower()
    if lower.endswith(".mdpeakid.tsv"):
        return "alignment_peak_id_matrix"
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
                block = _annotation_block(row)
                if block is not None:
                    _insert_annotation_result(
                        connection, run_id, "feature", feature_id, block,
                        artifact_ids[peak_path], row_number,
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
        representative_sample_names: dict[int, str] = {}
        for alignment_path in (path for path in relevant if path.suffix.lower() == ".mdalign"):
            for row_number, row in read_tsv(alignment_path, "Alignment ID\t"):
                alignment_id = parse_number(row.get("Alignment ID"), int)
                if alignment_id is None:
                    continue
                alignment_feature_id = make_id("alignment", run_id, alignment_id)
                alignments[alignment_id] = alignment_feature_id
                representative_name = row.get("Spectrum reference file name", "").strip()
                if representative_name and representative_name.lower() != "null":
                    representative_sample_names[alignment_id] = representative_name
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
                block = _annotation_block(row)
                if block is not None:
                    _insert_annotation_result(
                        connection, run_id, "alignment_feature", alignment_feature_id, block,
                        artifact_ids[alignment_path], row_number,
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

        peak_id_paths = [path for path in relevant if path.name.lower().endswith(".mdpeakid.tsv")]
        provenance_paths = [path for path in relevant if path.name.lower().endswith(".mdprovenance.tsv")]
        if alignments and not peak_id_paths and not provenance_paths:
            report.errors.append(
                "Alignment output exists but neither .mdpeakid.tsv nor .mdprovenance.tsv is present"
            )

        unknown_matrix_samples: set[str] = set()
        for peak_id_path in peak_id_paths:
            for row_number, row in read_tsv(peak_id_path):
                alignment_id = parse_number(row.get("alignment_master_id"), int)
                alignment_feature_id = alignments.get(alignment_id) if alignment_id is not None else None
                if alignment_feature_id is None:
                    report.errors.append(f"Unknown alignment ID {alignment_id} in {peak_id_path.name}")
                    continue
                representative_name = representative_sample_names.get(alignment_id)
                for file_id, (sample_name, raw_peak_id) in enumerate(
                    (item for item in row.items() if item[0] != "alignment_master_id")
                ):
                    sample_id = samples.get(sample_name)
                    if sample_id is None:
                        if sample_name not in unknown_matrix_samples:
                            report.errors.append(f"Unknown sample {sample_name!r} in {peak_id_path.name}")
                            unknown_matrix_samples.add(sample_name)
                        continue
                    source_peak_id = parse_number(raw_peak_id, int)
                    has_source_peak = source_peak_id is not None and source_peak_id >= 0
                    feature_id = features.get((sample_name, source_peak_id)) if has_source_peak else None
                    if has_source_peak and feature_id is None:
                        report.errors.append(
                            f"Alignment {alignment_id}, {sample_name} references missing peak {source_peak_id}"
                        )
                    is_representative = sample_name == representative_name
                    member_id = make_id("alignment-member", run_id, alignment_id, sample_name)
                    connection.execute(
                        """INSERT INTO alignment_member(
                            alignment_member_id, alignment_feature_id, sample_id, feature_id, file_id,
                            is_representative, has_source_peak, source_master_peak_id,
                            source_artifact_id, source_row
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(alignment_feature_id, sample_id) DO UPDATE SET
                            feature_id = excluded.feature_id,
                            file_id = excluded.file_id,
                            is_representative = excluded.is_representative,
                            has_source_peak = excluded.has_source_peak,
                            source_master_peak_id = excluded.source_master_peak_id,
                            source_artifact_id = excluded.source_artifact_id,
                            source_row = excluded.source_row""",
                        (
                            member_id, alignment_feature_id, sample_id, feature_id, file_id,
                            int(is_representative), int(has_source_peak),
                            source_peak_id if has_source_peak else None,
                            artifact_ids[peak_id_path], row_number,
                        ),
                    )
                    if is_representative:
                        connection.execute(
                            "UPDATE alignment_feature SET representative_sample_id = ?, representative_feature_id = ? WHERE alignment_feature_id = ?",
                            (sample_id, feature_id, alignment_feature_id),
                        )

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
                existing_member = connection.execute(
                    "SELECT source_master_peak_id FROM alignment_member WHERE alignment_feature_id = ? AND sample_id = ?",
                    (alignment_feature_id, sample_id),
                ).fetchone()
                if (
                    existing_member is not None
                    and existing_member["source_master_peak_id"] is not None
                    and source_peak_id is not None
                    and existing_member["source_master_peak_id"] != source_peak_id
                ):
                    report.errors.append(
                        f"Peak ID matrix/provenance mismatch for alignment {alignment_id}, {sample_name}: "
                        f"{existing_member['source_master_peak_id']} vs {source_peak_id}"
                    )
                connection.execute(
                    """INSERT INTO alignment_member(
                        alignment_member_id, alignment_feature_id, sample_id, feature_id, file_id,
                        is_representative, has_source_peak, source_master_peak_id, source_local_peak_id,
                        ms1_scan_index, ms2_scan_index, rt_min, mz, height, source_artifact_id, source_row
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(alignment_feature_id, sample_id) DO UPDATE SET
                        feature_id = excluded.feature_id,
                        file_id = excluded.file_id,
                        is_representative = excluded.is_representative,
                        has_source_peak = excluded.has_source_peak,
                        source_master_peak_id = excluded.source_master_peak_id,
                        source_local_peak_id = excluded.source_local_peak_id,
                        ms1_scan_index = excluded.ms1_scan_index,
                        ms2_scan_index = excluded.ms2_scan_index,
                        rt_min = excluded.rt_min,
                        mz = excluded.mz,
                        height = excluded.height""",
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
            sme_to_alignment: dict[str, str] = {}
            sme_owner: dict[str, str] = {}
            pending: list[tuple[int, str, dict[str, str]]] = []
            for row_number, section, row in read_mztab(mztab_path):
                pending.append((row_number, section, row))
                if section != "SMF":
                    continue
                record_id = row.get("SMF_ID", "")
                alignment_id = parse_number(record_id, int)
                if alignment_id not in alignments:
                    continue
                smf_to_alignment[record_id] = alignments[alignment_id]
                for sme_ref in split_refs(row.get("SME_ID_REFS", "")):
                    owner = sme_owner.get(sme_ref)
                    if owner is None:
                        sme_owner[sme_ref] = record_id
                        sme_to_alignment[sme_ref] = alignments[alignment_id]
                    elif owner != record_id:
                        report.warnings.append(
                            f"{mztab_path.name} SME_ID {sme_ref} is referenced by SMF_ID {owner} and "
                            f"{record_id}; keeping the {owner} link"
                        )
            unlinked_sme = sum(
                1 for _, section, row in pending
                if section == "SME" and row.get("SME_ID", "") not in sme_to_alignment
            )
            if unlinked_sme:
                report.warnings.append(
                    f"{unlinked_sme} SME rows in {mztab_path.name} could not be linked to an alignment feature"
                )
            for row_number, section, row in pending:
                id_field = f"{section}_ID"
                record_id = row.get(id_field, "")
                refs = row.get("SMF_ID_REFS", "") if section == "SML" else row.get("SME_ID_REFS", "")
                if section == "SMF":
                    alignment_feature_id = smf_to_alignment.get(record_id)
                elif section == "SME":
                    alignment_feature_id = sme_to_alignment.get(record_id)
                else:
                    alignment_feature_id = None
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
            "msdial_annotation_results": (
                "SELECT COUNT(*) FROM msdial_annotation_result WHERE run_id = ?",
                run_id,
            ),
            "mztab_sme_linked": (
                "SELECT COUNT(*) FROM mztab_record WHERE run_id = ? AND section = 'SME' AND alignment_feature_id IS NOT NULL",
                run_id,
            ),
        }
        for field_name, (query, parameter) in counts.items():
            setattr(report, field_name, connection.execute(query, (parameter,)).fetchone()[0])
    return report
