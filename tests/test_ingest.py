import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.ingest import ingest_run
from msdial_spectrum_catalog.validate import load_spectrum, validate_run


class IngestRunTests(unittest.TestCase):
    def test_ingests_feature_spectra_and_alignment_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample_a.mdpeak").write_text(
                "Peak ID\tName\tScan\tRT (min)\tPrecursor m/z\tHeight\tArea\tAdduct\n"
                "7\tFeature A\t100\t2.5\t300.1\t1000\t4000\t[M-H]-\n",
                encoding="utf-8",
            )
            (root / "analysis_files.csv").write_text(
                "file_path,file_name,file_type,class_id,acquisition_type,batch_order,analytical_order\n"
                "D:\\repository\\sample_a.raw,sample_a,Sample,Case,DDA,1,1\n",
                encoding="utf-8",
            )
            (root / "sample_a.mdmsp").write_text(
                "NAME: Feature A\nPRECURSORMZ: 300.1\nIONMODE: Negative\nRETENTIONTIME: 2.5\n"
                "COMMENT: |PEAKID=7|MS1SCAN=100|MS2SCAN=101\nNum Peaks: 2\n100 10\n150 20\n\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdalign").write_text(
                "\tClass\tSample\n\tFile type\tSample\n\tInjection order\t1\n\tBatch ID\t1\n"
                "Alignment ID\tAverage Rt(min)\tAverage Mz\tMetabolite name\tSpectrum reference file name\tsample_a\n"
                "3\t2.5\t300.1\tFeature A\tsample_a\t1000\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdmsp").write_text(
                "NAME: Feature A\nPRECURSORMZ: 300.1\nIONMODE: Negative\nRETENTIONTIME: 2.5\n"
                "COMMENT: |PEAKID=3\nNum Peaks: 2\n100 10\n150 20\n\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdprovenance.tsv").write_text(
                "alignment_master_id\talignment_local_id\tparent_alignment_id\tfile_id\tfile_name\t"
                "is_representative\thas_source_peak\tsource_master_peak_id\tsource_peak_id\t"
                "source_parent_peak_id\tms1_raw_spectrum_id\tms1_raw_spectrum_id_top\t"
                "ms2_raw_spectrum_id\tms2_raw_spectrum_ids\tms2_collision_energies\trt_min\tmz\t"
                "height\tarea_above_zero\tarea_above_baseline\n"
                "3\t3\t-1\t0\tsample_a\ttrue\ttrue\t7\t7\t-1\t100\t100\t101\t101\t101:20\t"
                "2.5\t300.1\t1000\t4000\t3500\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdpeakid.tsv").write_text(
                "alignment_master_id\tsample_a\n3\t7\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mzTab").write_text(
                "SMH\tSML_ID\tSMF_ID_REFS\tchemical_name\n"
                "SML\t3\t3\tFeature A\n"
                "SFH\tSMF_ID\texp_mass_to_charge\tretention_time_in_seconds\n"
                "SMF\t3\t300.1\t150\n",
                encoding="utf-8",
            )
            database = root / "catalog.sqlite"
            report = ingest_run(database, root, "mb_post", "MPST000007", "lcms_neg_dda")
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.samples, 1)
            self.assertEqual(report.features, 1)
            self.assertEqual(report.spectra, 2)
            self.assertEqual(report.alignments, 1)
            self.assertEqual(report.alignment_members, 1)
            validation = validate_run(database, report.run_id)
            self.assertTrue(validation.valid, validation.errors)
            self.assertEqual(validation.counts["mztab_small_molecule_features"], 1)
            spectrum_id = f"urn:msdial:spectrum:{report.run_id.replace(':', '%3A')}:sample_a:deconvoluted:7"
            spectrum = load_spectrum(database, spectrum_id)
            self.assertIsNotNone(spectrum)
            self.assertEqual(spectrum["payload"]["peaks"], [[100.0, 10.0], [150.0, 20.0]])
            connection = sqlite3.connect(database)
            try:
                linked = connection.execute(
                    "SELECT COUNT(*) FROM alignment_member WHERE feature_id IS NOT NULL AND is_representative = 1"
                ).fetchone()[0]
                self.assertEqual(linked, 1)
                scan_ids = connection.execute(
                    "SELECT ms1_scan_index, ms2_scan_index FROM alignment_member"
                ).fetchone()
                self.assertEqual(scan_ids, (100, 101))
                raw_path = connection.execute("SELECT raw_file_path FROM sample").fetchone()[0]
                self.assertEqual(raw_path, "D:\\repository\\sample_a.raw")
            finally:
                connection.close()

    def test_compact_peak_id_matrix_is_sufficient_without_detailed_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample_a.mdpeak").write_text(
                "Peak ID\tName\tScan\tRT (min)\tPrecursor m/z\tHeight\tArea\tAdduct\n"
                "7\tFeature A\t100\t2.5\t300.1\t1000\t4000\t[M-H]-\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdalign").write_text(
                "\tClass\tSample\n\tFile type\tSample\n\tInjection order\t1\n\tBatch ID\t1\n"
                "Alignment ID\tAverage Rt(min)\tAverage Mz\tMetabolite name\tSpectrum reference file name\tsample_a\n"
                "3\t2.5\t300.1\tFeature A\tsample_a\t1000\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdmsp").write_text(
                "NAME: Feature A\nPRECURSORMZ: 300.1\nIONMODE: Negative\nRETENTIONTIME: 2.5\n"
                "COMMENT: |PEAKID=3\nNum Peaks: 1\n100 10\n\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdpeakid.tsv").write_text(
                "alignment_master_id\tsample_a\n3\t7\n",
                encoding="utf-8",
            )

            report = ingest_run(root / "catalog.sqlite", root, "test", "S1", "unit")

            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.alignment_members, 1)
            connection = sqlite3.connect(root / "catalog.sqlite")
            try:
                row = connection.execute(
                    "SELECT source_master_peak_id, is_representative FROM alignment_member"
                ).fetchone()
                self.assertEqual(row, (7, 1))
            finally:
                connection.close()

    def test_detailed_provenance_remains_a_backward_compatible_link(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample_a.mdpeak").write_text(
                "Peak ID\tName\tScan\tRT (min)\tPrecursor m/z\tHeight\tArea\tAdduct\n"
                "7\tFeature A\t100\t2.5\t300.1\t1000\t4000\t[M-H]-\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdalign").write_text(
                "\tClass\tSample\n\tFile type\tSample\n\tInjection order\t1\n\tBatch ID\t1\n"
                "Alignment ID\tAverage Rt(min)\tAverage Mz\tMetabolite name\tsample_a\n"
                "3\t2.5\t300.1\tFeature A\t1000\n",
                encoding="utf-8",
            )
            (root / "AlignResult.mdprovenance.tsv").write_text(
                "alignment_master_id\tfile_id\tfile_name\tis_representative\thas_source_peak\t"
                "source_master_peak_id\tsource_peak_id\tms1_raw_spectrum_id_top\t"
                "ms2_raw_spectrum_id\trt_min\tmz\theight\n"
                "3\t0\tsample_a\ttrue\ttrue\t7\t7\t100\t101\t2.5\t300.1\t1000\n",
                encoding="utf-8",
            )

            report = ingest_run(root / "catalog.sqlite", root, "test", "S2", "unit")

            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.alignment_members, 1)

    def test_alignment_without_peak_id_link_is_invalid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "AlignResult.mdalign").write_text(
                "\tClass\n\tFile type\n\tInjection order\n\tBatch ID\n"
                "Alignment ID\tAverage Rt(min)\tAverage Mz\tMetabolite name\n"
                "0\t1.0\t100.0\tUnknown\n",
                encoding="utf-8",
            )
            report = ingest_run(root / "catalog.sqlite", root, "test", "S1", "unit")
            self.assertFalse(report.valid)
            self.assertIn(
                "Alignment output exists but neither .mdpeakid.tsv nor .mdprovenance.tsv is present",
                report.errors,
            )


if __name__ == "__main__":
    unittest.main()
