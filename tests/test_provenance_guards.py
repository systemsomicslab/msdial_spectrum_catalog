import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog import storage
from msdial_spectrum_catalog.ingest import ingest_run
from msdial_spectrum_catalog.validate import validate_run


PROVENANCE_HEADER = "\t".join([
    "alignment_master_id", "alignment_local_id", "parent_alignment_id", "file_id", "file_name",
    "is_representative", "has_source_peak", "peak_origin", "source_master_peak_id", "source_peak_id",
    "source_parent_peak_id", "ms1_raw_spectrum_id", "ms1_raw_spectrum_id_top", "ms2_raw_spectrum_id",
    "ms2_raw_spectrum_ids", "ms2_collision_energies", "rt_min", "mz", "height", "area_above_zero",
    "area_above_baseline",
])

# An older MS-DIAL provenance export, before AlignmentProvenanceExporter was fixed: no peak_origin
# column, a gap-filled row carrying scan 0 in every raw-spectrum field, and the mz column read from the
# chromatogram axis so the row that HAS a source peak reports the axis sentinel.
LEGACY_HEADER = "\t".join([
    "alignment_master_id", "file_id", "file_name", "is_representative", "has_source_peak",
    "source_master_peak_id", "source_peak_id", "ms1_raw_spectrum_id_top", "ms2_raw_spectrum_id",
    "rt_min", "mz", "height",
])


def _row(*fields) -> str:
    """One tab-separated provenance row, so a field count can be read off the call."""
    return "\t".join("" if field == "" else str(field) for field in fields) + "\n"


def _build_run(root: Path, provenance_rows: str, *, header: str = PROVENANCE_HEADER) -> None:
    (root / "sample_a.mdpeak").write_text(
        "Peak ID\tName\tScan\tRT (min)\tPrecursor m/z\tHeight\tArea\tAdduct\n"
        "7\tFeature A\t100\t2.5\t300.1\t1000\t4000\t[M-H]-\n",
        encoding="utf-8",
    )
    (root / "AlignResult.mdalign").write_text(
        "\tClass\n\tFile type\n\tInjection order\n\tBatch ID\n"
        "Alignment ID\tAverage Rt(min)\tAverage Mz\tMetabolite name\tsample_a\n"
        "3\t2.5\t300.1\tFeature A\t1000\n",
        encoding="utf-8",
    )
    (root / "AlignResult.mdprovenance.tsv").write_text(
        header + "\n" + provenance_rows,
        encoding="utf-8",
    )


class IngestGuardTests(unittest.TestCase):
    def test_a_member_without_a_source_peak_stores_no_scan_index(self):
        # The gap filler leaves the spectrum ids at 0, and 0 is a real scan number. Stored unguarded it
        # points an auditor at a spectrum belonging to a different peak.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, _row(
                3, 3, -1, 0, "sample_a", "false", "false", "gap_filled", -1, -1, -1,
                # An older export put scan 0 in all five raw-spectrum columns here.
                0, 0, 0, 0, "", 2.5, 300.1, 12, 34, "",
            ))
            report = ingest_run(root / "catalog.sqlite", root, "test", "S1", "unit")
            self.assertTrue(report.valid, report.errors)
            connection = sqlite3.connect(root / "catalog.sqlite")
            connection.row_factory = sqlite3.Row
            try:
                row = connection.execute("SELECT * FROM alignment_member").fetchone()
                self.assertEqual(row["has_source_peak"], 0)
                self.assertEqual(row["peak_origin"], "gap_filled")
                self.assertIsNone(row["ms1_scan_index"])
                self.assertIsNone(row["ms2_scan_index"])
                # The recovered m/z is a real measurement and is kept.
                self.assertAlmostEqual(row["mz"], 300.1)
            finally:
                connection.close()

    def test_a_negative_mz_is_a_sentinel_not_a_measurement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # The legacy layout, which is where a negative m/z on a row that HAS a source peak came from.
            _build_run(
                root,
                _row(3, 0, "sample_a", "true", "true", 7, 7, 100, 101, 2.5, -1, 1000),
                header=LEGACY_HEADER,
            )
            ingest_run(root / "catalog.sqlite", root, "test", "S1", "unit")
            connection = sqlite3.connect(root / "catalog.sqlite")
            try:
                mz = connection.execute("SELECT mz FROM alignment_member").fetchone()[0]
                self.assertIsNone(mz)
            finally:
                connection.close()

    def test_a_legacy_export_without_peak_origin_still_ingests(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(
                root,
                "3\t0\tsample_a\ttrue\ttrue\t7\t7\t100\t101\t2.5\t300.1\t1000\n",
                header=LEGACY_HEADER,
            )
            report = ingest_run(root / "catalog.sqlite", root, "test", "S1", "unit")
            self.assertTrue(report.valid, report.errors)
            connection = sqlite3.connect(root / "catalog.sqlite")
            try:
                origin = connection.execute("SELECT peak_origin FROM alignment_member").fetchone()[0]
                # Absent from the file, so absent from the row: the sentinel is all the old format had.
                self.assertIsNone(origin)
            finally:
                connection.close()


class ValidatorGuardTests(unittest.TestCase):
    """The checks that would have caught a systematically wrong provenance export."""

    RUN = "urn:msdial:run:test"

    def _minimal_run(self, database: Path, *, has_source_peak: int, mz, ms2_scan_index) -> None:
        storage.initialize(database)
        with storage.transaction(database) as connection:
            connection.execute("INSERT INTO study VALUES ('s', 'test', 'A', NULL)")
            connection.execute(
                "INSERT INTO analysis_unit VALUES ('u', 's', 'unit', NULL, NULL, NULL)"
            )
            connection.execute(
                "INSERT INTO analysis_run(run_id, analysis_unit_id, run_fingerprint, output_directory)"
                " VALUES (?, 'u', 'fp', '.')", (self.RUN,)
            )
            connection.execute(
                "INSERT INTO sample(sample_id, analysis_unit_id, sample_name) VALUES ('sm', 'u', 'a')"
            )
            connection.execute(
                "INSERT INTO alignment_feature(alignment_feature_id, run_id, alignment_master_id)"
                " VALUES ('af', ?, 0)", (self.RUN,)
            )
            connection.execute(
                """INSERT INTO alignment_member(
                    alignment_member_id, alignment_feature_id, sample_id, file_id,
                    is_representative, has_source_peak, mz, ms2_scan_index
                ) VALUES ('am', 'af', 'sm', 0, 1, ?, ?, ?)""",
                (has_source_peak, mz, ms2_scan_index),
            )

    def test_a_source_peak_without_a_usable_mz_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            self._minimal_run(database, has_source_peak=1, mz=-1.0, ms2_scan_index=None)
            report = validate_run(database, self.RUN)
            self.assertFalse(report.valid)
            self.assertTrue(
                any("no usable m/z" in error for error in report.errors),
                report.errors,
            )

    def test_a_gap_filled_member_carrying_a_scan_index_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            self._minimal_run(database, has_source_peak=0, mz=300.1, ms2_scan_index=0)
            report = validate_run(database, self.RUN)
            self.assertFalse(report.valid)
            self.assertTrue(
                any("carry a raw-spectrum index" in error for error in report.errors),
                report.errors,
            )

    def test_a_correct_member_passes_both_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            self._minimal_run(database, has_source_peak=0, mz=300.1, ms2_scan_index=None)
            report = validate_run(database, self.RUN)
            self.assertFalse(any("no usable m/z" in error for error in report.errors), report.errors)
            self.assertFalse(
                any("carry a raw-spectrum index" in error for error in report.errors), report.errors
            )


if __name__ == "__main__":
    unittest.main()
