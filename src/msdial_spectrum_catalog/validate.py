from __future__ import annotations

import json
import zlib
from dataclasses import dataclass, field

from .storage import connect


@dataclass
class ValidationReport:
    run_id: str
    counts: dict[str, int] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def validate_run(database, run_id: str) -> ValidationReport:
    report = ValidationReport(run_id)
    connection = connect(database)
    try:
        if connection.execute("SELECT 1 FROM analysis_run WHERE run_id = ?", (run_id,)).fetchone() is None:
            report.errors.append("Run does not exist")
            return report
        queries = {
            "samples": "SELECT COUNT(DISTINCT sample_id) FROM feature WHERE run_id = ?",
            "features": "SELECT COUNT(*) FROM feature WHERE run_id = ?",
            "sample_spectra": "SELECT COUNT(*) FROM spectrum WHERE run_id = ? AND spectrum_kind = 'deconvoluted'",
            "sample_ms2_spectra": "SELECT COUNT(*) FROM spectrum WHERE run_id = ? AND spectrum_kind = 'deconvoluted' AND peak_count > 0",
            "alignments": "SELECT COUNT(*) FROM alignment_feature WHERE run_id = ?",
            "consensus_spectra": "SELECT COUNT(*) FROM spectrum WHERE run_id = ? AND spectrum_kind = 'alignment_consensus'",
            "consensus_ms2_spectra": "SELECT COUNT(*) FROM spectrum WHERE run_id = ? AND spectrum_kind = 'alignment_consensus' AND peak_count > 0",
            "alignment_members": "SELECT COUNT(*) FROM alignment_member WHERE alignment_feature_id IN (SELECT alignment_feature_id FROM alignment_feature WHERE run_id = ?)",
            "mztab_small_molecule_features": "SELECT COUNT(*) FROM mztab_record WHERE run_id = ? AND section = 'SMF'",
            "mztab_small_molecule_summaries": "SELECT COUNT(*) FROM mztab_record WHERE run_id = ? AND section = 'SML'",
        }
        report.counts = {name: connection.execute(query, (run_id,)).fetchone()[0] for name, query in queries.items()}

        missing_sample_links = connection.execute(
            "SELECT COUNT(*) FROM spectrum WHERE run_id = ? AND spectrum_kind = 'deconvoluted' AND feature_id IS NULL",
            (run_id,),
        ).fetchone()[0]
        if missing_sample_links:
            report.errors.append(f"{missing_sample_links} sample spectra are not linked to mdpeak features")

        missing_member_links = connection.execute(
            """SELECT COUNT(*) FROM alignment_member
               WHERE has_source_peak = 1 AND feature_id IS NULL
               AND alignment_feature_id IN (SELECT alignment_feature_id FROM alignment_feature WHERE run_id = ?)""",
            (run_id,),
        ).fetchone()[0]
        if missing_member_links:
            report.errors.append(f"{missing_member_links} alignment members reference missing sample features")

        # A member with a source peak must carry a usable m/z. This check exists because an earlier
        # MS-DIAL provenance export read the column from the chromatogram axis instead of the peak mass,
        # so every member that DID have a source peak stored a sentinel and only the gap-filled ones
        # stored a real value. That run validated clean while being systematically wrong, which is the
        # failure this check makes impossible to repeat.
        unusable_member_mz = connection.execute(
            """SELECT COUNT(*) FROM alignment_member
               WHERE has_source_peak = 1 AND (mz IS NULL OR mz <= 0)
               AND alignment_feature_id IN (SELECT alignment_feature_id FROM alignment_feature WHERE run_id = ?)""",
            (run_id,),
        ).fetchone()[0]
        if unusable_member_mz:
            report.errors.append(
                f"{unusable_member_mz} alignment members have a source peak but no usable m/z"
            )

        # The mirror image: a member with no source peak has no source spectrum either, so a scan index
        # on such a row points at a spectrum belonging to something else.
        false_scan_pointers = connection.execute(
            """SELECT COUNT(*) FROM alignment_member
               WHERE has_source_peak = 0
               AND (ms1_scan_index IS NOT NULL OR ms2_scan_index IS NOT NULL)
               AND alignment_feature_id IN (SELECT alignment_feature_id FROM alignment_feature WHERE run_id = ?)""",
            (run_id,),
        ).fetchone()[0]
        if false_scan_pointers:
            report.errors.append(
                f"{false_scan_pointers} alignment members without a source peak carry a raw-spectrum index"
            )

        bad_representatives = connection.execute(
            """SELECT COUNT(*) FROM (
                SELECT a.alignment_feature_id, SUM(CASE WHEN m.is_representative = 1 THEN 1 ELSE 0 END) AS n
                FROM alignment_feature a LEFT JOIN alignment_member m USING(alignment_feature_id)
                WHERE a.run_id = ? GROUP BY a.alignment_feature_id HAVING n != 1
            )""",
            (run_id,),
        ).fetchone()[0]
        if bad_representatives:
            report.errors.append(f"{bad_representatives} alignments do not have exactly one representative member")

        if report.counts["alignments"] != report.counts["consensus_spectra"]:
            report.errors.append(
                f"Alignment/consensus count mismatch: {report.counts['alignments']} vs {report.counts['consensus_spectra']}"
            )
        unmapped_mztab = connection.execute(
            "SELECT COUNT(*) FROM mztab_record WHERE run_id = ? AND section IN ('SMF', 'SML') AND alignment_feature_id IS NULL",
            (run_id,),
        ).fetchone()[0]
        if unmapped_mztab:
            report.errors.append(f"{unmapped_mztab} mzTab-M SMF/SML records are not linked to alignment features")
        if report.counts["mztab_small_molecule_features"] and report.counts["mztab_small_molecule_features"] != report.counts["alignments"]:
            report.errors.append(
                "mzTab-M SMF/alignment count mismatch: "
                f"{report.counts['mztab_small_molecule_features']} vs {report.counts['alignments']}"
            )
        if report.counts["sample_ms2_spectra"] == 0:
            report.warnings.append("No sample deconvoluted spectrum contains fragment peaks")
        if report.counts["consensus_ms2_spectra"] == 0:
            report.warnings.append("No alignment consensus spectrum contains fragment peaks")
        return report
    finally:
        connection.close()


def load_spectrum(database, spectrum_id: str) -> dict | None:
    connection = connect(database)
    try:
        row = connection.execute(
            """SELECT s.*, b.compression, b.payload
               FROM spectrum s JOIN spectrum_blob b USING(payload_sha256)
               WHERE spectrum_id = ?""",
            (spectrum_id,),
        ).fetchone()
        if row is None:
            return None
        payload = bytes(row["payload"])
        if row["compression"] == "zlib-json":
            payload = zlib.decompress(payload)
        result = {key: row[key] for key in row.keys() if key not in {"payload", "compression"}}
        result["payload"] = json.loads(payload.decode("utf-8"))
        return result
    finally:
        connection.close()
