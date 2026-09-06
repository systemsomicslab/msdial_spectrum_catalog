import hashlib
import json
import sqlite3
import tempfile
import unittest
import zlib
from pathlib import Path

from msdial_spectrum_catalog.annotation import (
    AnnotationReport,
    CandidateInput,
    EvidenceInput,
    ToolRunInput,
    criteria_rule_id_for,
    list_assertions,
    load_assertion,
    notation_for,
    record_assertion,
    record_criteria_set,
    record_tool_run,
    validate_annotations,
)
from msdial_spectrum_catalog.identifiers import make_id
from msdial_spectrum_catalog.storage import connect, initialize

REPOSITORY = "mb_post"
ACCESSION = "MPST000007"
UNIT = "lcms_neg_dda"

STUDY_ID = make_id("study", REPOSITORY, ACCESSION)
UNIT_ID = make_id("unit", REPOSITORY, ACCESSION, UNIT)
RUN_ID = make_id("run", REPOSITORY, ACCESSION, UNIT, "fingerprint")
SAMPLE_ID = make_id("sample", REPOSITORY, ACCESSION, UNIT, "sample_a")
ALIGNMENT_ID = make_id("alignment", RUN_ID, 3)
SPECTRUM_ID = make_id("spectrum", RUN_ID, "alignment_consensus", 3)


def _put_blob(connection, payload):
    raw = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode("utf-8")
    payload_sha = hashlib.sha256(raw).hexdigest()
    connection.execute(
        "INSERT OR IGNORE INTO spectrum_blob(payload_sha256, compression, uncompressed_bytes, payload) "
        "VALUES (?, 'zlib-json', ?, ?)",
        (payload_sha, len(raw), zlib.compress(raw, level=9)),
    )
    return payload_sha


def build_catalog(root: Path) -> Path:
    """Create the smallest catalog that can anchor an annotation claim."""
    database = root / "catalog.sqlite"
    initialize(database)
    connection = connect(database)
    try:
        connection.execute(
            "INSERT INTO study(study_id, repository, accession, title) VALUES (?, ?, ?, ?)",
            (STUDY_ID, REPOSITORY, ACCESSION, "Lenticin worked example"),
        )
        connection.execute(
            "INSERT INTO analysis_unit(analysis_unit_id, study_id, external_unit_id, separation_type, "
            "ion_mode, acquisition_type) VALUES (?, ?, ?, ?, ?, ?)",
            (UNIT_ID, STUDY_ID, UNIT, "lc", "negative", "dda"),
        )
        connection.execute(
            "INSERT INTO analysis_run(run_id, analysis_unit_id, run_fingerprint, output_directory) "
            "VALUES (?, ?, ?, ?)",
            (RUN_ID, UNIT_ID, "fingerprint", str(root)),
        )
        connection.execute(
            "INSERT INTO sample(sample_id, analysis_unit_id, sample_name) VALUES (?, ?, ?)",
            (SAMPLE_ID, UNIT_ID, "sample_a"),
        )
        connection.execute(
            "INSERT INTO alignment_feature(alignment_feature_id, run_id, alignment_master_id, "
            "average_rt_min, average_mz, name) VALUES (?, ?, ?, ?, ?, ?)",
            (ALIGNMENT_ID, RUN_ID, 3, 2.5, 300.1, "Unknown"),
        )
        payload_sha = _put_blob(connection, {"peaks": [[100.0, 10.0], [150.0, 20.0]]})
        connection.execute(
            "INSERT INTO spectrum(spectrum_id, run_id, sample_id, alignment_feature_id, spectrum_kind, "
            "precursor_mz, rt_min, ion_mode, peak_count, payload_sha256) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (SPECTRUM_ID, RUN_ID, SAMPLE_ID, ALIGNMENT_ID, "alignment_consensus", 300.1, 2.5, "Negative", 2, payload_sha),
        )
        connection.commit()
    finally:
        connection.close()
    return database


def add_reference_library(database: Path, library_kind: str) -> str:
    library_id = make_id("library", "demo", library_kind)
    reference_spectrum_id = make_id("reference-spectrum", library_id, 1)
    connection = connect(database)
    try:
        payload_sha = _put_blob(connection, {"peaks": [[100.0, 10.0]]})
        connection.execute(
            "INSERT INTO reference_library(library_id, library_name, library_version, library_kind, sha256) "
            "VALUES (?, ?, ?, ?, ?)",
            (library_id, f"demo {library_kind}", "1.0", library_kind, "0" * 64),
        )
        connection.execute(
            "INSERT INTO reference_spectrum(reference_spectrum_id, library_id, library_record_index, "
            "record_name, peak_count, payload_sha256) VALUES (?, ?, ?, ?, ?, ?)",
            (reference_spectrum_id, library_id, 1, "Lenticin analogue", 1, payload_sha),
        )
        connection.commit()
    finally:
        connection.close()
    return reference_spectrum_id


class AnnotationWriteTests(unittest.TestCase):
    def test_report_has_the_house_shape(self):
        report = AnnotationReport()
        self.assertTrue(report.valid)
        self.assertEqual((report.assertions, report.candidates, report.evidence), (0, 0, 0))
        report.errors.append("boom")
        self.assertFalse(report.valid)

    def test_records_a_minimal_assertion_against_a_minimal_catalog(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L4",
                evidence=[EvidenceInput(tag="FM", metric="mass_error_ppm", measured_value=1.2, measured_unit="ppm")],
                alignment_feature_id=ALIGNMENT_ID,
                formula="C20H32O2",
            )
            stored = load_assertion(database, assertion_id)
            self.assertIsNotNone(stored)
            self.assertEqual(stored["annotation_level"], "L4")
            self.assertEqual(stored["notation_verbatim"], "L4[FM]")
            self.assertEqual(stored["subject_kind"], "alignment_feature")
            self.assertEqual(stored["alignment_feature_id"], ALIGNMENT_ID)
            self.assertIsNone(stored["claim_concept_id"])
            self.assertEqual(stored["claim_unresolved"], 0)
            self.assertEqual(len(stored["evidence"]), 1)
            self.assertEqual(stored["evidence"][0]["evidence_tag"], "FM")
            self.assertEqual(stored["evidence"][0]["evidence_concept_id"], "smb:evidence/formula_mass")
            self.assertEqual(notation_for(database, assertion_id), "L4[FM]")

            listed = list_assertions(database, run_id=RUN_ID)
            self.assertEqual([item["assertion_id"] for item in listed], [assertion_id])
            self.assertEqual(list_assertions(database, spectrum_id=SPECTRUM_ID), listed)
            self.assertEqual(list_assertions(database, alignment_feature_id=ALIGNMENT_ID), listed)

            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertEqual((report.assertions, report.evidence, report.candidates), (1, 1, 0))

    def test_lenticin_case_three_is_a_substructure_complete_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("SC",),
                evidence=[
                    EvidenceInput(tag="FM"),
                    EvidenceInput(tag="DF"),
                    EvidenceInput(tag="MN"),
                    EvidenceInput(tag="CO"),
                ],
                compound_name="lenticin",
                curation_comment="Case 3: components resolved, connectivity unresolved",
            )
            stored = load_assertion(database, assertion_id)
            self.assertEqual(stored["notation_verbatim"], "L3-SC[FM,DF,MN,CO]")
            self.assertEqual(stored["claim_concept_id"], "smb:claim/substructure_complete")
            self.assertEqual(stored["annotation_claim"], "SC")
            self.assertEqual(stored["claim_unresolved"], 0)
            self.assertEqual(
                [item["claim_concept_id"] for item in stored["claim_components"]],
                ["smb:claim/substructure_complete"],
            )
            self.assertEqual(
                [item["evidence_concept_id"] for item in stored["evidence"]],
                [
                    "smb:evidence/formula_mass",
                    "smb:evidence/diagnostic_fragment",
                    "smb:evidence/molecular_network",
                    "smb:evidence/contextual",
                ],
            )
            self.assertNotIn("smb:evidence/reference_standard", [item["evidence_concept_id"] for item in stored["evidence"]])
            self.assertEqual(notation_for(database, assertion_id), "L3-SC[FM,DF,MN,CO]")

    def test_evidence_order_does_not_change_the_canonical_notation(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("SC",),
                evidence=[EvidenceInput(tag="CO"), EvidenceInput(tag="MN"), EvidenceInput(tag="DF"), EvidenceInput(tag="FM")],
            )
            self.assertEqual(notation_for(database, assertion_id), "L3-SC[FM,DF,MN,CO]")

    def test_supplied_notation_is_stored_verbatim_and_a_mismatch_raises(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("SC",),
                evidence=[EvidenceInput(tag="CO"), EvidenceInput(tag="FM")],
                notation="L3-SC[CO,FM]",
            )
            stored = load_assertion(database, assertion_id)
            self.assertEqual(stored["notation_verbatim"], "L3-SC[CO,FM]")
            self.assertEqual(stored["notation"], "L3-SC[FM,CO]")
            with self.assertRaises(ValueError):
                record_assertion(
                    database,
                    spectrum_id=SPECTRUM_ID,
                    level="L3",
                    claim_tokens=("SC",),
                    evidence=[EvidenceInput(tag="FM")],
                    notation="L3-SC[FM,DF]",
                )

    def test_claim_token_on_level_two_raises(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            with self.assertRaises(ValueError):
                record_assertion(
                    database,
                    spectrum_id=SPECTRUM_ID,
                    level="L2",
                    claim_tokens=("SP",),
                    evidence=[EvidenceInput(tag="SL"), EvidenceInput(tag="FM")],
                )

    def test_unknown_level_raises_and_is_never_inferred(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            with self.assertRaises(ValueError):
                record_assertion(
                    database,
                    spectrum_id=SPECTRUM_ID,
                    level="L1.5",
                    evidence=[EvidenceInput(tag="RS")],
                )

    def test_level_three_without_a_claim_is_incomplete_not_invalid(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                evidence=[EvidenceInput(tag="FM"), EvidenceInput(tag="DF")],
            )
            stored = load_assertion(database, assertion_id)
            self.assertEqual(stored["claim_unresolved"], 1)
            self.assertEqual(stored["claim_components"], [])
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertTrue(
                any("no claim component" in warning for warning in report.warnings), report.warnings
            )

    def test_the_rebound_cp_token_is_stored_as_two_different_concepts(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            draft_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("CP",),
                evidence=[EvidenceInput(tag="DF")],
                vocab_version="smb-v1-draft",
            )
            consensus_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("CP",),
                evidence=[EvidenceInput(tag="DF")],
                vocab_version="smb-v2-consensus",
            )
            self.assertNotEqual(draft_id, consensus_id)
            draft = load_assertion(database, draft_id)
            consensus = load_assertion(database, consensus_id)
            self.assertEqual(draft["claim_concept_id"], "smb:claim/substructure_complete")
            self.assertEqual(consensus["claim_concept_id"], "smb:claim/class")
            self.assertEqual(draft["annotation_claim"], "CP")
            self.assertEqual(consensus["annotation_claim"], "CP")
            self.assertEqual(notation_for(database, draft_id, vocab_version="smb-v2-consensus"), "L3-SC[DF]")
            self.assertEqual(notation_for(database, consensus_id), "L3-CP[DF]")
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.assertions, 2)

    def test_unknown_claim_token_is_strict_error_and_quarantines_under_quarantine_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            with self.assertRaises(ValueError):
                record_assertion(
                    database,
                    spectrum_id=SPECTRUM_ID,
                    level="L3",
                    claim_tokens=("MO",),
                    evidence=[EvidenceInput(tag="DF")],
                )
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("MO",),
                evidence=[EvidenceInput(tag="DF"), EvidenceInput(tag="CX")],
                mode="quarantine",
            )
            stored = load_assertion(database, assertion_id)
            self.assertEqual(stored["claim_unresolved"], 1)
            self.assertEqual(stored["annotation_claim"], "MO")
            self.assertEqual(stored["claim_components"], [])
            tags = {item["evidence_tag"]: item["evidence_concept_id"] for item in stored["evidence"]}
            self.assertEqual(tags["DF"], "smb:evidence/diagnostic_fragment")
            self.assertIsNone(tags["CX"])
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertTrue(any("no resolved concept" in warning for warning in report.warnings))


class AnnotationCandidateTests(unittest.TestCase):
    def _candidates(self):
        return [
            CandidateInput(
                rank=1,
                compound_name="isomer A",
                formula="C20H32O2",
                inchikey="AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                smiles="CCO",
                score=0.81,
                score_type="squared_cosine",
                score_gap_to_next=0.01,
            ),
            CandidateInput(
                rank=2,
                compound_name="isomer B",
                formula="C20H32O2",
                inchikey="DDDDDDDDDDDDDD-EEEEEEEEEE-F",
                score=0.80,
                score_type="squared_cosine",
            ),
        ]

    def test_a_or_b_without_an_ambiguity_class_warns(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("SP",),
                evidence=[EvidenceInput(tag="SL"), EvidenceInput(tag="FM")],
                candidates=self._candidates(),
            )
            stored = load_assertion(database, assertion_id)
            self.assertEqual(stored["candidate_count"], 2)
            self.assertEqual(stored["claim_concept_id"], "smb:claim/structure")
            self.assertEqual([item["rank"] for item in stored["candidates"]], [1, 2])
            self.assertEqual(stored["candidates"][0]["inchikey_skeleton"], "AAAAAAAAAAAAAA")
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertTrue(
                any("without an ambiguity_class_id" in warning for warning in report.warnings),
                report.warnings,
            )

    def test_a_or_b_with_an_ambiguity_class_does_not_warn(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            ambiguity_class_id = make_id("ambiguity-class", "demo", 1)
            connection = connect(database)
            try:
                connection.execute(
                    "INSERT INTO ambiguity_class(ambiguity_class_id, definition_id, library_scope_json, "
                    "member_count, linkage_rule, min_pairwise_score, score_convention) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (ambiguity_class_id, "demo-definition", json.dumps(["demo"]), 2, "complete_linkage", 0.9, "squared_cosine"),
                )
                connection.commit()
            finally:
                connection.close()
            record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("SP",),
                evidence=[EvidenceInput(tag="SL"), EvidenceInput(tag="FM")],
                candidates=self._candidates(),
                ambiguity_class_id=ambiguity_class_id,
            )
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertFalse(
                [warning for warning in report.warnings if "ambiguity_class_id" in warning],
                report.warnings,
            )

    def test_non_contiguous_candidate_ranks_raise(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            with self.assertRaises(ValueError):
                record_assertion(
                    database,
                    spectrum_id=SPECTRUM_ID,
                    level="L3",
                    claim_tokens=("SP",),
                    candidates=[CandidateInput(rank=1), CandidateInput(rank=3)],
                )
            with self.assertRaises(ValueError):
                record_assertion(
                    database,
                    spectrum_id=SPECTRUM_ID,
                    level="L3",
                    claim_tokens=("SP",),
                    candidates=[CandidateInput(rank=1), CandidateInput(rank=1)],
                )

    def test_score_convention_round_trips_without_conversion(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L2",
                evidence=[
                    EvidenceInput(
                        tag="SL",
                        metric="weighted_dot_product",
                        measured_value=0.64,
                        comparison=">=",
                        threshold_value=0.36,
                        passed=True,
                        value={"score_convention": "squared_cosine", "msdial_cutoff": "0.6F*0.6F"},
                    )
                ],
                candidates=[
                    CandidateInput(rank=1, compound_name="hit", score=0.64, score_type="squared_cosine")
                ],
            )
            stored = load_assertion(database, assertion_id)
            evidence_row = stored["evidence"][0]
            self.assertEqual(evidence_row["measured_value"], 0.64)
            self.assertEqual(evidence_row["threshold_value"], 0.36)
            self.assertEqual(json.loads(evidence_row["evidence_value_json"])["score_convention"], "squared_cosine")
            self.assertEqual(stored["candidates"][0]["score"], 0.64)
            self.assertEqual(stored["candidates"][0]["score_type"], "squared_cosine")


class AnnotationIdempotenceTests(unittest.TestCase):
    def test_recording_the_identical_assertion_twice_updates_in_place(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            arguments = dict(
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("SC",),
                evidence=[EvidenceInput(tag="FM"), EvidenceInput(tag="DF")],
                candidates=[CandidateInput(rank=1, compound_name="lenticin")],
                compound_name="lenticin",
            )
            first = record_assertion(database, **arguments)
            second = record_assertion(database, **arguments)
            self.assertEqual(first, second)
            connection = sqlite3.connect(database)
            try:
                counts = {
                    table: connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                    for table in (
                        "annotation_assertion",
                        "annotation_claim_component",
                        "annotation_evidence",
                        "annotation_candidate",
                    )
                }
            finally:
                connection.close()
            self.assertEqual(counts["annotation_assertion"], 1)
            self.assertEqual(counts["annotation_claim_component"], 1)
            self.assertEqual(counts["annotation_evidence"], 2)
            self.assertEqual(counts["annotation_candidate"], 1)

    def test_tool_run_is_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            tool_run = ToolRunInput(
                tool_name="MS-DIAL",
                tool_version="5.5.250101",
                parameters={"ms2_tolerance_da": 0.05},
                provenance={"console": "MsdialConsoleApp"},
                input_fingerprint="deadbeef",
            )
            first = record_tool_run(database, tool_run, run_id=RUN_ID)
            second = record_tool_run(database, tool_run, run_id=RUN_ID)
            self.assertEqual(first, second)
            connection = sqlite3.connect(database)
            try:
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM annotation_tool_run").fetchone()[0], 1
                )
                row = connection.execute(
                    "SELECT tool_name, tool_version, status FROM annotation_tool_run"
                ).fetchone()
            finally:
                connection.close()
            self.assertEqual(row, ("MS-DIAL", "5.5.250101", "completed"))


class CriteriaSetTests(unittest.TestCase):
    def test_criteria_rules_are_queryable_and_link_to_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            criteria_set_id = make_id("criteria-set", "pilot", "v1")
            record_criteria_set(
                database,
                criteria_set_id,
                "LC-MS annotation pilot",
                version="1.0",
                description="Study-level operational criteria",
                rules=[
                    {
                        "tag": "SL",
                        "metric": "weighted_dot_product",
                        "comparison": ">=",
                        "threshold_value": 0.36,
                        "operational_criterion": "MS-DIAL squared cosine above the 0.6F*0.6F cutoff",
                        "notes": "squared_cosine convention",
                    },
                    {
                        "tag": "FM",
                        "metric": "mass_error_ppm",
                        "comparison": "<=",
                        "threshold_value": 5.0,
                        "threshold_unit": "ppm",
                    },
                ],
            )
            connection = connect(database)
            try:
                rows = connection.execute(
                    "SELECT evidence_token, evidence_concept_id, metric, comparison, threshold_value "
                    "FROM criteria_rule WHERE criteria_set_id = ? ORDER BY evidence_token",
                    (criteria_set_id,),
                ).fetchall()
                stored_version = connection.execute(
                    "SELECT vocab_version FROM criteria_set WHERE criteria_set_id = ?", (criteria_set_id,)
                ).fetchone()[0]
            finally:
                connection.close()
            self.assertEqual(
                [tuple(row) for row in rows],
                [
                    ("FM", "smb:evidence/formula_mass", "mass_error_ppm", "<=", 5.0),
                    ("SL", "smb:evidence/spectral_library", "weighted_dot_product", ">=", 0.36),
                ],
            )
            self.assertEqual(stored_version, "smb-v2-consensus")

            rule_id = criteria_rule_id_for(criteria_set_id, "smb:evidence/spectral_library", "weighted_dot_product")
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L2",
                evidence=[
                    EvidenceInput(
                        tag="SL",
                        metric="weighted_dot_product",
                        measured_value=0.72,
                        comparison=">=",
                        threshold_value=0.36,
                        passed=True,
                        criteria_rule_id=rule_id,
                    )
                ],
                criteria_set_id=criteria_set_id,
            )
            stored = load_assertion(database, assertion_id)
            evidence_row = stored["evidence"][0]
            self.assertEqual(evidence_row["criteria_rule_id"], rule_id)
            self.assertEqual(evidence_row["measured_value"], 0.72)
            self.assertEqual(evidence_row["threshold_value"], 0.36)
            self.assertEqual(evidence_row["passed"], 1)
            self.assertEqual(stored["criteria_set_id"], criteria_set_id)

            record_criteria_set(
                database,
                criteria_set_id,
                "LC-MS annotation pilot",
                rules=[{"tag": "FM", "metric": "mass_error_ppm", "threshold_value": 3.0}],
            )
            connection = connect(database)
            try:
                remaining = connection.execute(
                    "SELECT COUNT(*) FROM criteria_rule WHERE criteria_set_id = ?", (criteria_set_id,)
                ).fetchone()[0]
            finally:
                connection.close()
            self.assertEqual(remaining, 1)

    def test_unknown_criteria_tag_raises(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            with self.assertRaises(ValueError):
                record_criteria_set(
                    database,
                    make_id("criteria-set", "pilot", "bad"),
                    "bad",
                    rules=[{"tag": "CX"}],
                )


class AnnotationValidationTests(unittest.TestCase):
    def test_in_silico_library_cited_as_spectral_library_evidence_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            predicted_id = add_reference_library(database, "in_silico_predicted")
            record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L2",
                evidence=[
                    EvidenceInput(tag="SL", source_reference_spectrum_id=predicted_id),
                    EvidenceInput(tag="FM"),
                ],
            )
            report = validate_annotations(database, RUN_ID)
            self.assertFalse(report.valid)
            self.assertTrue(
                any("in_silico_predicted" in error for error in report.errors), report.errors
            )

    def test_experimental_library_cited_as_spectral_library_evidence_is_fine(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            measured_id = add_reference_library(database, "experimental_reference")
            record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L2",
                evidence=[
                    EvidenceInput(tag="SL", source_reference_spectrum_id=measured_id),
                    EvidenceInput(tag="FM"),
                ],
            )
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)

    def test_failing_evidence_row_warns(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L2",
                evidence=[
                    EvidenceInput(tag="SL", metric="weighted_dot_product", measured_value=0.2, passed=False),
                ],
            )
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertTrue(any("did not pass" in warning for warning in report.warnings), report.warnings)

    def test_missing_run_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            report = validate_annotations(database, "urn:msdial:run:absent")
            self.assertFalse(report.valid)
            self.assertIn("Run does not exist", report.errors)

    def test_dangling_spectrum_reference_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            assertion_id = record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L4",
                evidence=[EvidenceInput(tag="FM")],
                alignment_feature_id=ALIGNMENT_ID,
            )
            # Foreign keys are enforced by the write API, so the only way to reach this check is to
            # simulate a catalog whose spectrum row was removed outside the API.
            connection = sqlite3.connect(database)
            try:
                connection.execute("PRAGMA foreign_keys = OFF")
                connection.execute(
                    "UPDATE annotation_assertion SET spectrum_id = ? WHERE assertion_id = ?",
                    ("urn:msdial:spectrum:gone", assertion_id),
                )
                connection.commit()
            finally:
                connection.close()
            report = validate_annotations(database, RUN_ID)
            self.assertFalse(report.valid)
            self.assertTrue(any("is missing from spectrum" in error for error in report.errors), report.errors)

    def test_discouraged_claim_combination_warns(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_catalog(Path(directory))
            record_assertion(
                database,
                spectrum_id=SPECTRUM_ID,
                level="L3",
                claim_tokens=("SP", "CP"),
                evidence=[EvidenceInput(tag="SL"), EvidenceInput(tag="FM")],
            )
            report = validate_annotations(database, RUN_ID)
            self.assertTrue(report.valid, report.errors)
            self.assertTrue(report.warnings)


if __name__ == "__main__":
    unittest.main()
