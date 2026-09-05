import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.ingest import ingest_run
from msdial_spectrum_catalog.validate import validate_run


CANDIDATE_COLUMNS = (
    "alignment_master_id", "alignment_local_id", "parent_alignment_id", "candidate_rank",
    "candidate_count", "is_representative", "annotator_id", "database_id", "source", "priority",
    "library_id", "name", "formula", "ontology", "inchikey", "smiles", "reference_mz",
    "reference_rt_min", "reference_adduct", "annotation_tag_vs1", "is_reference_matched",
    "is_annotation_suggested", "is_precursor_mz_match", "is_spectrum_match",
    "is_spectrum_comparison_performed", "total_score", "mz_similarity", "rt_similarity",
    "ri_similarity", "ccs_similarity", "isotope_similarity", "simple_dot_product",
    "weighted_dot_product", "reverse_dot_product", "matched_peaks_count", "matched_peaks_percentage",
)
CANDIDATE_HEADER = "\t".join(CANDIDATE_COLUMNS)


def _candidate(**overrides) -> str:
    """One row of the sidecar, written from named fields so a call reads as what it asserts."""
    values = {name: "" for name in CANDIDATE_COLUMNS}
    values.update({
        "alignment_master_id": 3,
        "alignment_local_id": 3,
        "parent_alignment_id": -1,
        "candidate_rank": 1,
        "candidate_count": 1,
        "is_representative": "true",
        "annotator_id": "msp",
        "database_id": "MspDB",
        "source": "MspDB",
        "priority": 1,
        "library_id": 7,
        "name": "Feature A",
        "is_reference_matched": "true",
        "is_annotation_suggested": "false",
        "is_precursor_mz_match": "true",
        "is_spectrum_match": "true",
        "is_spectrum_comparison_performed": "true",
        "total_score": 0.75,
        "simple_dot_product": 0.5,
        "weighted_dot_product": 0.25,
        "reverse_dot_product": 0.75,
        "matched_peaks_count": 12,
        "matched_peaks_percentage": 0.5,
    })
    values.update(overrides)
    return "\t".join(str(values[name]) for name in CANDIDATE_COLUMNS) + "\n"


def _build_run(root: Path, candidate_rows: str, *, alignment_name: str = "Feature A") -> None:
    (root / "sample_a.mdpeak").write_text(
        "Peak ID\tName\tScan\tRT (min)\tPrecursor m/z\tHeight\tArea\tAdduct\n"
        "7\tFeature A\t100\t2.5\t300.1\t1000\t4000\t[M-H]-\n",
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
        f"3\t2.5\t300.1\t{alignment_name}\tsample_a\t1000\n",
        encoding="utf-8",
    )
    (root / "AlignResult.mdmsp").write_text(
        "NAME: Feature A\nPRECURSORMZ: 300.1\nIONMODE: Negative\nRETENTIONTIME: 2.5\n"
        "COMMENT: |PEAKID=3\nNum Peaks: 2\n100 10\n150 20\n\n",
        encoding="utf-8",
    )
    (root / "AlignResult.mdpeakid.tsv").write_text(
        "alignment_master_id\tsample_a\n3\t7\n",
        encoding="utf-8",
    )
    (root / "AlignResult.mdprovenance.tsv").write_text(
        "alignment_master_id\talignment_local_id\tparent_alignment_id\tfile_id\tfile_name\t"
        "is_representative\thas_source_peak\tpeak_origin\tsource_master_peak_id\tsource_peak_id\t"
        "source_parent_peak_id\tms1_raw_spectrum_id\tms1_raw_spectrum_id_top\t"
        "ms2_raw_spectrum_id\tms2_raw_spectrum_ids\tms2_collision_energies\trt_min\tmz\t"
        "height\tarea_above_zero\tarea_above_baseline\n"
        "3\t3\t-1\t0\tsample_a\ttrue\ttrue\tdetected\t7\t7\t-1\t100\t100\t101\t101\t101:20\t"
        "2.5\t300.1\t1000\t4000\t3500\n",
        encoding="utf-8",
    )
    (root / "AlignResult.mdcandidate.tsv").write_text(
        CANDIDATE_HEADER + "\n" + candidate_rows,
        encoding="utf-8",
    )


def _ingest(root: Path):
    return ingest_run(root / "catalog.sqlite", root, "test", "S1", "unit")


def _rows(root: Path, sql: str = "SELECT * FROM msdial_annotation_candidate ORDER BY candidate_rank"):
    connection = sqlite3.connect(root / "catalog.sqlite")
    connection.row_factory = sqlite3.Row
    try:
        return [dict(row) for row in connection.execute(sql)]
    finally:
        connection.close()


class CandidateIngestTests(unittest.TestCase):
    def test_every_candidate_is_stored_with_its_rank(self):
        # The whole point: an alignment feature whose search kept three references it could not separate
        # is three rows, not one, and the two it used to drop are the ones that make the claim honest.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, "".join([
                _candidate(candidate_rank=1, candidate_count=3, is_representative="true",
                           name="Feature A", library_id=7, total_score=0.75),
                _candidate(candidate_rank=2, candidate_count=3, is_representative="false",
                           name="Feature B", library_id=8, total_score=0.74),
                _candidate(candidate_rank=3, candidate_count=3, is_representative="false",
                           name="Feature C", library_id=9, total_score=0.73),
            ]))
            report = _ingest(root)
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.msdial_annotation_candidates, 3)

            rows = _rows(root)
            self.assertEqual([row["name"] for row in rows], ["Feature A", "Feature B", "Feature C"])
            self.assertEqual([row["is_representative"] for row in rows], [1, 0, 0])
            self.assertEqual({row["candidate_count"] for row in rows}, {3})
            self.assertEqual(rows[0]["subject_kind"], "alignment_feature")
            self.assertEqual(rows[0]["score_convention"], "cosine")

            validation = validate_run(root / "catalog.sqlite", report.run_id)
            self.assertTrue(validation.valid, validation.errors)
            self.assertEqual(validation.counts["annotation_candidates"], 3)
            self.assertEqual(validation.counts["ambiguous_alignments"], 1)

    def test_a_single_candidate_is_not_an_ambiguity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, _candidate())
            report = _ingest(root)
            self.assertTrue(report.valid, report.errors)
            validation = validate_run(root / "catalog.sqlite", report.run_id)
            self.assertEqual(validation.counts["annotated_alignments"], 1)
            self.assertEqual(validation.counts["ambiguous_alignments"], 0)

    def test_an_uncompared_spectrum_stores_no_spectral_score(self):
        # The exporter writes those cells empty. Storing 0 instead would recreate exactly the defect the
        # .mdpeak score columns used to carry: a comparison that never happened, rendered as a number.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, _candidate(
                is_reference_matched="false", is_annotation_suggested="true",
                is_spectrum_match="false", is_spectrum_comparison_performed="false",
                simple_dot_product="", weighted_dot_product="", reverse_dot_product="",
                matched_peaks_count="", matched_peaks_percentage="",
            ))
            report = _ingest(root)
            self.assertTrue(report.valid, report.errors)
            row = _rows(root)[0]
            for column in ("simple_dot_product", "weighted_dot_product", "reverse_dot_product",
                           "matched_peaks_count", "matched_peaks_percentage"):
                self.assertIsNone(row[column], column)
            # The aggregate score and the precursor term are still measurements.
            self.assertAlmostEqual(row["total_score"], 0.75)
            self.assertIsNone(row["score_convention"])
            self.assertEqual(row["is_annotation_suggested"], 1)

    def test_a_score_without_a_comparison_is_refused_rather_than_stored(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, _candidate(is_spectrum_comparison_performed="false"))
            report = _ingest(root)
            self.assertFalse(report.valid)
            self.assertTrue(
                any("no spectrum comparison was performed" in error for error in report.errors),
                report.errors,
            )
            # Refused as evidence, not merely reported: the contradicted columns are still empty.
            self.assertIsNone(_rows(root)[0]["simple_dot_product"])

    def test_a_rank_gap_fails_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, "".join([
                _candidate(candidate_rank=1, candidate_count=3),
                _candidate(candidate_rank=3, candidate_count=3, is_representative="false"),
            ]))
            report = _ingest(root)
            validation = validate_run(root / "catalog.sqlite", report.run_id)
            self.assertFalse(validation.valid)
            self.assertTrue(
                any("ranked 1..candidate_count" in error for error in validation.errors),
                validation.errors,
            )

    def test_a_representative_below_rank_one_fails_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, "".join([
                _candidate(candidate_rank=1, candidate_count=2, is_representative="false"),
                _candidate(candidate_rank=2, candidate_count=2, is_representative="true"),
            ]))
            report = _ingest(root)
            validation = validate_run(root / "catalog.sqlite", report.run_id)
            self.assertFalse(validation.valid)
            self.assertTrue(
                any("rank the representative first" in error for error in validation.errors),
                validation.errors,
            )

    def test_the_two_artifacts_must_name_the_same_winner(self):
        # .mdalign and the sidecar describe one decision from different sides. If they disagree about
        # which candidate won, one of them is describing a different run.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, _candidate(name="Feature Z"), alignment_name="Feature A")
            report = _ingest(root)
            validation = validate_run(root / "catalog.sqlite", report.run_id)
            self.assertFalse(validation.valid)
            self.assertTrue(
                any("name a different winner" in error for error in validation.errors),
                validation.errors,
            )

    def test_the_mdalign_suggestion_prefix_is_not_a_disagreement(self):
        # .mdalign renders a precursor-only suggestion as "no MS2: <name>"; the sidecar publishes the
        # raw reference name. Those are the same winner.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(
                root,
                _candidate(
                    is_reference_matched="false", is_annotation_suggested="true",
                    is_spectrum_match="false", is_spectrum_comparison_performed="false",
                    simple_dot_product="", weighted_dot_product="", reverse_dot_product="",
                    matched_peaks_count="", matched_peaks_percentage="",
                ),
                alignment_name="no MS2: Feature A",
            )
            report = _ingest(root)
            self.assertTrue(report.valid, report.errors)
            validation = validate_run(root / "catalog.sqlite", report.run_id)
            self.assertTrue(validation.valid, validation.errors)

    def test_an_in_house_identifier_is_not_a_compound_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(
                root,
                _candidate(name="RIKEN N-VS1 ID-97 from Mouse_Muscle_WT_CTX0_Ctr"),
                alignment_name="RIKEN N-VS1 ID-97 from Mouse_Muscle_WT_CTX0_Ctr",
            )
            report = _ingest(root)
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(_rows(root)[0]["candidate_is_named"], 0)

    def test_a_run_without_the_sidecar_still_ingests(self):
        # The export is off by default, so its absence is the normal case rather than a fault.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, _candidate())
            (root / "AlignResult.mdcandidate.tsv").unlink()
            report = _ingest(root)
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.msdial_annotation_candidates, 0)
            validation = validate_run(root / "catalog.sqlite", report.run_id)
            self.assertTrue(validation.valid, validation.errors)


if __name__ == "__main__":
    unittest.main()
