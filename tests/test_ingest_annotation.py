import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.ingest import (
    classify_metabolite_name,
    ingest_run,
    is_annotated_metabolite_name,
)


PEAK_COLUMNS = [
    "Peak ID", "Name", "Scan", "RT left(min)", "RT (min)", "RT right (min)", "Precursor m/z",
    "Height", "Area", "Model masses", "Adduct", "Isotope", "Comment", "Reference RT",
    "Reference m/z", "Formula", "Ontology", "InChIKey", "SMILES", "Annotation tag (VS1.0)",
    "RT matched", "m/z matched", "MS/MS matched", "RT similarity", "m/z similarity",
    "Simple dot product", "Weighted dot product", "Reverse dot product", "Matched peaks count",
    "Matched peaks percentage", "Total score", "S/N", "MS1 isotopes", "MSMS spectrum",
]

ALIGN_COLUMNS = [
    "Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type",
    "Post curation result", "Fill %", "MS/MS assigned", "Reference RT", "Reference m/z",
    "Formula", "Ontology", "INCHIKEY", "SMILES", "Annotation tag (VS1.0)", "RT matched",
    "m/z matched", "MS/MS matched", "Comment", "Manually modified for quantification",
    "Manually modified for annotation", "Isotope tracking parent ID",
    "Isotope tracking weight number", "RT similarity", "m/z similarity", "Simple dot product",
    "Weighted dot product", "Reverse dot product", "Matched peaks count",
    "Matched peaks percentage", "Total score", "S/N average", "Spectrum reference file name",
    "MS1 isotopic spectrum", "MS/MS spectrum",
]


def _row(columns: list[str], values: dict[str, str]) -> str:
    return "\t".join(values.get(column, "null") for column in columns)


def _peak_row(peak_id: int, name: str, **overrides: str) -> str:
    values = {
        "Peak ID": str(peak_id), "Name": name, "Scan": str(100 + peak_id),
        "RT (min)": "2.5", "Precursor m/z": "300.1", "Height": "1000", "Area": "4000",
        "Adduct": "[M-H]-", "Comment": "", "Formula": "C18H30O4", "Ontology": "Fatty acid",
        "InChIKey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "SMILES": "O=C(O)CCCC",
        "Annotation tag (VS1.0)": "410", "RT matched": "False", "m/z matched": "True",
        "MS/MS matched": "True", "RT similarity": "0.00", "m/z similarity": "1.00",
        "Simple dot product": "0.125", "Weighted dot product": "0.153",
        "Reverse dot product": "0.866", "Matched peaks count": "1.00",
        "Matched peaks percentage": "2.00", "Total score": "2.216", "S/N": "12.7",
    }
    values.update(overrides)
    return _row(PEAK_COLUMNS, values)


def _align_row(alignment_id: int, name: str, **overrides: str) -> str:
    values = {
        "Alignment ID": str(alignment_id), "Average Rt(min)": "2.5", "Average Mz": "300.1",
        "Metabolite name": name, "Adduct type": "[M-H]-", "Fill %": "1.00",
        "MS/MS assigned": "True", "Formula": "C18H30O4", "Ontology": "Fatty acid",
        "INCHIKEY": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "SMILES": "O=C(O)CCCC",
        "Annotation tag (VS1.0)": "430", "RT matched": "False", "m/z matched": "True",
        "MS/MS matched": "True", "Comment": "Normalized unit: Intensity",
        "RT similarity": "0.00", "m/z similarity": "1.00", "Simple dot product": "0.616",
        "Weighted dot product": "0.612", "Reverse dot product": "0.839",
        "Matched peaks count": "2.00", "Matched peaks percentage": "1.00",
        "Total score": "1.763", "S/N average": "12.70",
        "Spectrum reference file name": "sample_a",
    }
    values.update(overrides)
    return _row(ALIGN_COLUMNS, values) + "\t1000"


def _build_run(root: Path, *, sme_refs: tuple[str, str] = ("10|11", "12")) -> None:
    (root / "sample_a.mdpeak").write_text(
        "\t".join(PEAK_COLUMNS) + "\n"
        + _peak_row(1, "FA 18:3;O2") + "\n"
        + _peak_row(2, "Unknown", **{"Formula": "null", "Annotation tag (VS1.0)": "999"}) + "\n"
        + _peak_row(3, "w/o MS2:FA 5:0") + "\n"
        + _peak_row(4, "RIKEN N-VS1 ID-1988 from Mouse_Feces_WT_N_Ctr") + "\n"
        + _peak_row(5, "") + "\n",
        encoding="utf-8",
    )
    (root / "sample_a.mdmsp").write_text(
        "NAME: FA 18:3;O2\nPRECURSORMZ: 300.1\nIONMODE: Negative\nRETENTIONTIME: 2.5\n"
        "COMMENT: |PEAKID=1|MS1SCAN=101|MS2SCAN=102\nNum Peaks: 2\n100 10\n150 20\n\n",
        encoding="utf-8",
    )
    preamble = "".join(
        "\t" * (len(ALIGN_COLUMNS) - 1) + f"{label}\tSample\n"
        for label in ("Class", "File type", "Injection order", "Batch ID")
    )
    (root / "AlignResult-1.mdalign").write_text(
        preamble
        + "\t".join(ALIGN_COLUMNS) + "\tsample_a\n"
        + _align_row(0, "Unknown", **{"Formula": "null", "Annotation tag (VS1.0)": "999"}) + "\n"
        + _align_row(1, "ST 24:2;O4") + "\n",
        encoding="utf-8",
    )
    (root / "AlignResult-1.mdmsp").write_text(
        "NAME: ST 24:2;O4\nPRECURSORMZ: 300.1\nIONMODE: Negative\nRETENTIONTIME: 2.5\n"
        "COMMENT: |PEAKID=1\nNum Peaks: 1\n100 10\n\n",
        encoding="utf-8",
    )
    (root / "AlignResult-1.mdpeakid.tsv").write_text(
        "alignment_master_id\tsample_a\n0\t-1\n1\t1\n",
        encoding="utf-8",
    )
    (root / "AlignResult-1.mzTab").write_text(
        "SMH\tSML_ID\tSMF_ID_REFS\tchemical_name\n"
        "SML\t0\t0\tUnknown\n"
        "SML\t1\t1\tST 24:2;O4\n"
        "SFH\tSMF_ID\tSME_ID_REFS\tSME_ID_REF_ambiguity_code\tadduct_ion\texp_mass_to_charge\n"
        f"SMF\t0\t{sme_refs[0]}\tnull\t[M-H]1-\t300.1\n"
        f"SMF\t1\t{sme_refs[1]}\tnull\t[M-H]1-\t300.1\n"
        "SEH\tSME_ID\tdatabase_identifier\tchemical_formula\tchemical_name\tadduct_ion\t"
        "theoretical_mass_to_charge\tspectra_ref\tid_confidence_measure[1]\n"
        "SME\t10\tLbmDB:ST 24:2;O4\tC24H40O4\tST 24:2;O4\t[M-H]1-\t300.101\tms_run[1]:ms1scanID=1\t1.7\n"
        "SME\t11\tLbmDB:ST 24:2;O5\tC24H40O5\tST 24:2;O5\t[M-H]1-\t300.102\tms_run[1]:ms1scanID=1\t1.6\n"
        "SME\t12\tLbmDB:FA 5:0\tC5H10O2\tFA 5:0\t[M-H]1-\t300.103\tms_run[1]:ms1scanID=1\t2.5\n",
        encoding="utf-8",
    )


class MetaboliteNamePredicateTests(unittest.TestCase):
    def test_msdial_placeholder_names_are_not_annotations(self):
        for name in (None, "", "   ", "null", "Unknown", "unknown", "Unknown metabolite"):
            self.assertFalse(is_annotated_metabolite_name(name), name)

    def test_real_names_are_annotations(self):
        for name in ("FA 18:3;O2", "ST 24:2;O4", "low score: Cer 18:1;2O/24:0", "HexCer 18:1;O2/16:0"):
            self.assertTrue(is_annotated_metabolite_name(name), name)

    def test_suggestion_prefixes_are_classified_rather_than_discarded(self):
        # MsdialCore/Utility/DataAccess.cs SetMoleculeMsPropertyAsSuggested writes "no MS2: " when no
        # product-ion spectrum was acquired and "low score: " when one was acquired but failed the search.
        # Both are annotations; only their evidence strength differs, so they must survive ingestion with
        # that difference recorded.
        for name in ("no MS2: FA 5:0", "low score: NAGly 11:0", "RIKEN N-VS1 ID-1988"):
            self.assertTrue(is_annotated_metabolite_name(name), name)

    def test_classify_metabolite_name(self):
        self.assertEqual(classify_metabolite_name("FA 5:0"), ("msms_matched", "FA 5:0", True))
        self.assertEqual(classify_metabolite_name("no MS2: FA 5:0"), ("precursor_only", "FA 5:0", True))
        # MS-DIAL 4 omits the space after the colon; MS-DIAL 5 includes it.
        self.assertEqual(classify_metabolite_name("w/o MS2:FA 5:0"), ("precursor_only", "FA 5:0", True))
        self.assertEqual(
            classify_metabolite_name("low score: NAGly 11:0"), ("low_score", "NAGly 11:0", True)
        )
        self.assertEqual(
            classify_metabolite_name("no MS2: RIKEN N-VS1 ID-1988"),
            ("precursor_only", "RIKEN N-VS1 ID-1988", False),
        )
        self.assertEqual(
            classify_metabolite_name("RIKEN N-VS1 ID-1988"),
            ("msms_matched", "RIKEN N-VS1 ID-1988", False),
        )


class MsdialAnnotationResultTests(unittest.TestCase):
    def test_annotated_rows_are_promoted_with_squared_cosine_convention(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root)
            database = root / "catalog.sqlite"
            report = ingest_run(database, root, "test", "S1", "unit")
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.msdial_annotation_results, 4)
            connection = sqlite3.connect(database)
            connection.row_factory = sqlite3.Row
            try:
                feature_row = connection.execute(
                    "SELECT * FROM msdial_annotation_result WHERE subject_kind = 'feature'"
                ).fetchone()
                self.assertEqual(feature_row["metabolite_name"], "FA 18:3;O2")
                self.assertEqual(feature_row["formula"], "C18H30O4")
                self.assertEqual(feature_row["ontology"], "Fatty acid")
                self.assertEqual(feature_row["inchikey"], "AAAAAAAAAAAAAA-BBBBBBBBBB-C")
                self.assertEqual(feature_row["smiles"], "O=C(O)CCCC")
                self.assertEqual(feature_row["adduct"], "[M-H]-")
                self.assertEqual(feature_row["annotation_tag"], "410")
                self.assertEqual(feature_row["is_rt_matched"], 0)
                self.assertEqual(feature_row["is_mz_matched"], 1)
                self.assertEqual(feature_row["is_msms_matched"], 1)
                self.assertAlmostEqual(feature_row["simple_dot_product"], 0.125)
                self.assertAlmostEqual(feature_row["weighted_dot_product"], 0.153)
                self.assertAlmostEqual(feature_row["reverse_dot_product"], 0.866)
                self.assertAlmostEqual(feature_row["matched_peaks_count"], 1.0)
                self.assertAlmostEqual(feature_row["matched_peaks_percentage"], 2.0)
                self.assertAlmostEqual(feature_row["total_score"], 2.216)
                self.assertIsNone(feature_row["ccs_similarity"])
                self.assertEqual(feature_row["score_convention"], "squared_cosine")
                self.assertEqual(feature_row["annotation_kind"], "msms_matched")
                self.assertEqual(feature_row["candidate_name"], "FA 18:3;O2")
                self.assertEqual(feature_row["candidate_is_named"], 1)
                self.assertEqual(feature_row["rank"], 1)
                self.assertEqual(feature_row["annotator_id"], "")
                self.assertIsNotNone(feature_row["source_artifact_id"])
                self.assertEqual(feature_row["source_row"], 2)
                self.assertEqual(feature_row["subject_id"], feature_row["feature_id"])
                self.assertIsNone(feature_row["alignment_feature_id"])

                align_row = connection.execute(
                    "SELECT * FROM msdial_annotation_result WHERE subject_kind = 'alignment_feature'"
                ).fetchone()
                self.assertEqual(align_row["metabolite_name"], "ST 24:2;O4")
                self.assertEqual(align_row["annotation_tag"], "430")
                self.assertAlmostEqual(align_row["reverse_dot_product"], 0.839)
                self.assertEqual(align_row["score_convention"], "squared_cosine")
                self.assertEqual(align_row["comment"], "Normalized unit: Intensity")
                self.assertEqual(align_row["subject_id"], align_row["alignment_feature_id"])
                self.assertIsNone(align_row["feature_id"])
            finally:
                connection.close()

    def test_unknown_produces_no_row_while_suggestions_are_classified(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root)
            database = root / "catalog.sqlite"
            ingest_run(database, root, "test", "S1", "unit")
            connection = sqlite3.connect(database)
            try:
                rows = {
                    name: (kind, candidate, is_named)
                    for name, kind, candidate, is_named in connection.execute(
                        "SELECT metabolite_name, annotation_kind, candidate_name, candidate_is_named "
                        "FROM msdial_annotation_result ORDER BY metabolite_name"
                    )
                }
                # 'Unknown' asserts nothing and must not become a row.
                self.assertNotIn("Unknown", rows)
                self.assertEqual(rows["FA 18:3;O2"], ("msms_matched", "FA 18:3;O2", 1))
                self.assertEqual(rows["ST 24:2;O4"], ("msms_matched", "ST 24:2;O4", 1))
                # A precursor-only suggestion is retained, but marked so it can never be read as an
                # MS/MS library match.
                self.assertEqual(rows["w/o MS2:FA 5:0"], ("precursor_only", "FA 5:0", 1))
                # An in-house record whose NAME is an identifier is retained and flagged as unnamed.
                self.assertEqual(
                    rows["RIKEN N-VS1 ID-1988 from Mouse_Feces_WT_N_Ctr"],
                    ("msms_matched", "RIKEN N-VS1 ID-1988 from Mouse_Feces_WT_N_Ctr", 0),
                )
            finally:
                connection.close()

    def test_not_applicable_negative_sentinels_become_null(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            header = [
                "Peak ID", "Name", "Scan", "RT (min)", "Precursor m/z", "Height", "Area", "Adduct",
                "Simple dot product", "Weighted dot product", "Reverse dot product",
                "Matched peaks count", "Matched peaks percentage", "Total score",
            ]
            values = [
                "7", "no MS2: FA 18:3;O2", "100", "2.5", "300.1", "1000", "4000", "[M-H]-",
                "0.000", "0.000", "0.000",
                "-1.00", "-1.00", "-0.146",
            ]
            (root / "sample_a.mdpeak").write_text(
                "\t".join(header) + "\n" + "\t".join(values) + "\n",
                encoding="utf-8",
            )
            database = root / "catalog.sqlite"
            report = ingest_run(database, root, "test", "S4", "unit")
            self.assertTrue(report.valid, report.errors)
            connection = sqlite3.connect(database)
            connection.row_factory = sqlite3.Row
            try:
                row = connection.execute("SELECT * FROM msdial_annotation_result").fetchone()
                self.assertEqual(row["annotation_kind"], "precursor_only")
                # -1 is MS-DIAL's not-applicable marker; a count and a percentage cannot be negative.
                self.assertIsNone(row["matched_peaks_count"])
                self.assertIsNone(row["matched_peaks_percentage"])
                # .mdpeak writes 0.000 dot products for a precursor-only row while .mdalign writes null.
                # 0.000 is the unset MsScanMatchResult default, not a comparison that scored zero.
                self.assertIsNone(row["simple_dot_product"])
                self.assertIsNone(row["weighted_dot_product"])
                self.assertIsNone(row["reverse_dot_product"])
                # Total score is an unnormalized composite, so its negative value is kept verbatim and
                # no cosine convention is claimed for the row.
                self.assertAlmostEqual(row["total_score"], -0.146)
                self.assertIsNone(row["score_convention"])
            finally:
                connection.close()

    def test_peak_table_without_score_columns_still_ingests(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample_a.mdpeak").write_text(
                "Peak ID\tName\tScan\tRT (min)\tPrecursor m/z\tHeight\tArea\tAdduct\tFormula\n"
                "7\tFA 18:3;O2\t100\t2.5\t300.1\t1000\t4000\t[M-H]-\tC18H30O4\n",
                encoding="utf-8",
            )
            database = root / "catalog.sqlite"
            report = ingest_run(database, root, "test", "S3", "unit")
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.msdial_annotation_results, 1)
            connection = sqlite3.connect(database)
            connection.row_factory = sqlite3.Row
            try:
                row = connection.execute("SELECT * FROM msdial_annotation_result").fetchone()
                self.assertEqual(row["metabolite_name"], "FA 18:3;O2")
                self.assertEqual(row["formula"], "C18H30O4")
                self.assertIsNone(row["simple_dot_product"])
                self.assertIsNone(row["total_score"])
                self.assertIsNone(row["is_msms_matched"])
                self.assertIsNone(row["score_convention"])
            finally:
                connection.close()

    def test_reingesting_the_same_directory_does_not_duplicate_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root)
            database = root / "catalog.sqlite"
            first = ingest_run(database, root, "test", "S1", "unit")
            second = ingest_run(database, root, "test", "S1", "unit")
            self.assertTrue(second.valid, second.errors)
            self.assertEqual(first.run_id, second.run_id)
            self.assertEqual(second.msdial_annotation_results, first.msdial_annotation_results)
            self.assertEqual(second.mztab_sme_linked, first.mztab_sme_linked)
            connection = sqlite3.connect(database)
            try:
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM msdial_annotation_result").fetchone()[0], 4
                )
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM mztab_record").fetchone()[0], 7
                )
            finally:
                connection.close()


class MztabEvidenceLinkTests(unittest.TestCase):
    def test_every_sme_row_is_linked_to_an_alignment_feature(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root)
            database = root / "catalog.sqlite"
            report = ingest_run(database, root, "test", "S1", "unit")
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.mztab_sme_linked, 3)
            connection = sqlite3.connect(database)
            try:
                unlinked = connection.execute(
                    "SELECT COUNT(*) FROM mztab_record WHERE section = 'SME' AND alignment_feature_id IS NULL"
                ).fetchone()[0]
                self.assertEqual(unlinked, 0)
                grouped = dict(
                    connection.execute(
                        "SELECT record_id, alignment_feature_id FROM mztab_record WHERE section = 'SME'"
                    )
                )
                self.assertEqual(grouped["10"], grouped["11"])
                self.assertNotEqual(grouped["10"], grouped["12"])
            finally:
                connection.close()
            self.assertEqual(
                [warning for warning in report.warnings if "SME" in warning], []
            )

    def test_unlinked_sme_rows_are_reported_as_a_warning(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, sme_refs=("10|11", "null"))
            report = ingest_run(root / "catalog.sqlite", root, "test", "S1", "unit")
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.mztab_sme_linked, 2)
            self.assertTrue(
                any("1 SME rows" in warning and "could not be linked" in warning for warning in report.warnings),
                report.warnings,
            )

    def test_sme_shared_by_two_smf_rows_warns_and_keeps_the_first_link(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _build_run(root, sme_refs=("10,12", "10"))
            database = root / "catalog.sqlite"
            report = ingest_run(database, root, "test", "S1", "unit")
            self.assertTrue(report.valid, report.errors)
            self.assertTrue(
                any(
                    "SME_ID 10 is referenced by SMF_ID 0 and 1" in warning
                    for warning in report.warnings
                ),
                report.warnings,
            )
            connection = sqlite3.connect(database)
            try:
                grouped = dict(
                    connection.execute(
                        "SELECT record_id, alignment_feature_id FROM mztab_record WHERE section = 'SME'"
                    )
                )
                self.assertEqual(grouped["10"], grouped["12"])
                self.assertIsNone(grouped["11"])
            finally:
                connection.close()


if __name__ == "__main__":
    unittest.main()
