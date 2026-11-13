from contextlib import closing
import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.ingest import ingest_run
from msdial_spectrum_catalog.validate import validate_run
from ws1_fixture import write_unit


class RepresentativeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.unit = write_unit(self.root)
        self.db = self.root / "catalog.sqlite"

    def ingest(self, **kwargs):
        return ingest_run(self.db, self.unit / "output", "synthetic", "SYN001", "unit-a", **kwargs)

    def test_representative_features_spectra_members_and_validation(self):
        report = self.ingest(granularity="representative")
        self.assertTrue(report.valid, report.errors)
        self.assertEqual((report.samples, report.features, report.spectra, report.alignment_members), (2, 2, 2, 4))
        self.assertEqual(report.representative_spectra_nonempty, 1)
        self.assertGreater(report.db_bytes_after, 0)
        validation = validate_run(self.db, report.run_id)
        self.assertTrue(validation.valid, validation.errors)
        self.assertEqual(validation.counts["sample_spectra"], 0)
        self.assertNotIn("No sample deconvoluted spectrum contains fragment peaks", validation.warnings)
        with closing(sqlite3.connect(self.db)) as c, c:
            self.assertEqual(c.execute("SELECT COUNT(*) FROM msdial_annotation_result WHERE subject_kind='feature'").fetchone()[0], 0)
            self.assertEqual(c.execute("SELECT COUNT(*) FROM alignment_feature WHERE representative_feature_id IS NOT NULL").fetchone()[0], 2)

    def test_upgrade_same_run_and_reverse_noop(self):
        rep = self.ingest(granularity="representative")
        full = self.ingest(granularity="full")
        self.assertEqual(rep.run_id, full.run_id)
        self.assertEqual((full.features, full.spectra), (4, 4))
        with closing(sqlite3.connect(self.db)) as c, c:
            before = list(c.iterdump())
        reverse = self.ingest(granularity="representative")
        self.assertEqual(reverse.granularity, "full")
        self.assertIn("already full", reverse.warnings)
        with closing(sqlite3.connect(self.db)) as c, c:
            self.assertEqual(before, list(c.iterdump()))
        self.assertTrue(validate_run(self.db, rep.run_id).valid)

    def test_matrix_only_members_have_unexported_mz_coverage(self):
        (self.unit / "output/AlignResult.mdprovenance.tsv").unlink()
        report = self.ingest(granularity="representative")
        validation = validate_run(self.db, report.run_id)
        self.assertTrue(validation.valid, validation.errors)
        self.assertEqual(validation.counts["member_mz_not_exported"], 4)
        self.assertIn("4 member m/z not exported", validation.warnings)

    def test_dual_lipid_winner_and_supported_name_priority(self):
        sidecar = self.unit / "output/AlignResult.mdcandidate.tsv"
        original = sidecar.read_text()
        for index, (name, supported, expected) in enumerate((
            ("Invented alpha|Invented beta", "", True),
            ("Invented beta|Invented alpha", "", True),
            ("Invented alpha|Invented beta", "Invented beta", False),
            ("Invented beta|Invented gamma", "Invented alpha", True),
            ("Invented alphabet|Invented beta", "", False),
        )):
            with self.subTest(name=name, supported=supported):
                lines = original.replace("Invented alpha", name).splitlines()
                sidecar.write_text(lines[0] + "\tname_at_supported_level\tis_lipid_chains_match\n"
                                   + lines[1] + "\t" + supported + "\ttrue\n")
                self.db = self.root / f"lipid-{index}.sqlite"
                report = self.ingest(granularity="representative")
                self.assertEqual(validate_run(self.db, report.run_id).valid, expected)
                with closing(sqlite3.connect(self.db)) as c:
                    self.assertEqual(c.execute("SELECT name, name_at_supported_level, is_lipid_chains_match FROM msdial_annotation_candidate").fetchone(),
                                     (name, supported or None, 1))

    def test_optional_members_and_mztab(self):
        report = self.ingest(granularity="representative", members=False, mztab="none")
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.alignment_members, 0)
        validation = validate_run(self.db, report.run_id)
        self.assertTrue(validation.valid, validation.errors)
        self.assertEqual(validation.counts["mztab_small_molecule_features"], 0)

    def test_full_dual_name_is_accepted_as_published(self):
        name = "Invented alpha|Invented beta"
        for filename in ("AlignResult.mdcandidate.tsv", "AlignResult.mdalign"):
            path = self.unit / "output" / filename
            path.write_text(path.read_text().replace("Invented alpha", name))
        report = self.ingest(granularity="representative")
        self.assertTrue(validate_run(self.db, report.run_id).valid)

    def test_blank_candidate_evidence_is_null_in_inventory(self):
        from msdial_spectrum_catalog.catalog_build import build_report
        path = self.unit / "output/AlignResult.mdcandidate.tsv"
        path.write_text(path.read_text().replace("ReferenceSpectrum", "").replace(" mz; msms ", "   "))
        self.ingest(granularity="representative")
        with closing(sqlite3.connect(self.db)) as c:
            self.assertEqual(c.execute("SELECT evidence_source, measured_terms FROM msdial_annotation_candidate").fetchone(),
                             (None, None))
        self.assertEqual(build_report(self.db)["candidate_evidence_null"], 1)

    def test_missing_representative_is_invalid(self):
        report = self.ingest(granularity="representative")
        with closing(sqlite3.connect(self.db)) as c, c:
            c.execute("UPDATE alignment_feature SET representative_feature_id=NULL")
        self.assertFalse(validate_run(self.db, report.run_id).valid)

    def test_candidate_evidence_verbatim_and_missing_columns(self):
        report = self.ingest(granularity="representative")
        with closing(sqlite3.connect(self.db)) as c, c:
            self.assertEqual(c.execute("SELECT evidence_source, measured_terms FROM msdial_annotation_candidate").fetchone(),
                             ("ReferenceSpectrum", " mz; msms "))
        sidecar = self.unit / "output/AlignResult.mdcandidate.tsv"
        lines = ["\t".join(line.split("\t")[:-2]) for line in sidecar.read_text().splitlines()]
        sidecar.write_text("\n".join(lines) + "\n")
        second = self.ingest(granularity="full")
        with closing(sqlite3.connect(self.db)) as c, c:
            self.assertEqual(c.execute("SELECT evidence_source, measured_terms FROM msdial_annotation_candidate WHERE run_id=?", (second.run_id,)).fetchone(), (None, None))

    def test_invalid_options_write_nothing(self):
        with self.assertRaises(ValueError):
            self.ingest(granularity="unsupported")
        self.assertFalse(self.db.exists())

    def test_feature_candidates_full_only_with_evidence(self):
        path = self.unit / "output/sample_a.mdcandidate.tsv"
        path.write_text("subject_kind\tfile_name\tpeak_id\tcandidate_rank\tcandidate_count\tis_representative\tname\tevidence_source\tmeasured_terms\n"
                        "feature\tsample_a\t7\t1\t1\ttrue\tInvented alpha\tInSilico\tformula\n")
        rep = self.ingest(granularity="representative")
        self.assertTrue(rep.valid, rep.errors)
        full = self.ingest(granularity="full")
        self.assertTrue(full.valid, full.errors)
        with closing(sqlite3.connect(self.db)) as c:
            self.assertEqual(c.execute("SELECT evidence_source, measured_terms FROM msdial_annotation_candidate WHERE subject_kind='feature'").fetchone(), ("InSilico", "formula"))
