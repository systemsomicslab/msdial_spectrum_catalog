import hashlib
import json
import sqlite3
import tempfile
import unittest
import zlib
from pathlib import Path

from msdial_spectrum_catalog.ambiguity import (
    AmbiguityReport,
    ClassDefinition,
    ambiguity_class_for,
    blocking_key,
    collision_energy_bin,
    compute_ambiguity_classes,
    discriminating_evidence_needed,
    discriminating_ions,
    is_clique,
    iter_blocks,
)
from msdial_spectrum_catalog.identifiers import make_id
from msdial_spectrum_catalog.storage import connect, initialize

LIBRARY_ID = make_id("library", "public-demo", "experimental")

MZ = [100.1, 120.1, 140.1, 160.1, 180.1, 200.1, 220.1]

# Three spectra with A ~ B and B ~ C above the default thresholds while A ~ C is below both of them.
# Measured with similarity.compare at the default tolerance: A-B 0.950/0.941, B-C 0.962/0.938,
# A-C 0.828/0.771. This is the non-transitivity the anchored design exists to preserve.
INTENSITIES_A = [1.00, 0.90, 0.80, 0.45, 0.25, 0.15, 0.08]
INTENSITIES_B = [0.75, 0.72, 0.70, 0.62, 0.55, 0.50, 0.45]
INTENSITIES_C = [0.30, 0.35, 0.50, 0.70, 0.85, 0.95, 1.00]
INTENSITIES_A_NEAR = [0.98, 0.91, 0.78, 0.46, 0.24, 0.16, 0.08]


def spectrum(intensities: list[float]) -> list[list[float]]:
    return [[mz, value] for mz, value in zip(MZ, intensities) if value > 0]


def _put_blob(connection, payload) -> str:
    raw = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode("utf-8")
    payload_sha = hashlib.sha256(raw).hexdigest()
    connection.execute(
        "INSERT OR IGNORE INTO spectrum_blob(payload_sha256, compression, uncompressed_bytes, payload) "
        "VALUES (?, 'zlib-json', ?, ?)",
        (payload_sha, len(raw), zlib.compress(raw, level=9)),
    )
    return payload_sha


def reference_id(index: int) -> str:
    return make_id("reference-spectrum", LIBRARY_ID, index)


def build_library(root: Path, records: list[dict]) -> Path:
    """Insert reference_library and reference_spectrum rows directly, with no dependency on the loader."""
    database = root / "catalog.sqlite"
    initialize(database)
    connection = connect(database)
    try:
        connection.execute(
            "INSERT OR IGNORE INTO reference_library(library_id, library_name, library_version, "
            "library_kind, source_uri, license, record_count) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (
                LIBRARY_ID, "public demo positive", "vDemo", "experimental",
                "https://example.invalid/demo.msp", "CC-BY", len(records),
            ),
        )
        for record in records:
            payload_sha = _put_blob(connection, {"peaks": record["peaks"]})
            inchikey = record.get("inchikey")
            connection.execute(
                """INSERT INTO reference_spectrum(
                       reference_spectrum_id, library_id, library_record_index, record_name, inchikey,
                       inchikey_skeleton, smiles, formula, ontology, precursor_mz, precursor_type,
                       ion_mode, instrument_type, instrument_class, collision_energy_raw,
                       collision_energy_value, rt_min, ccs, peak_count, payload_sha256)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    reference_id(record["index"]),
                    LIBRARY_ID,
                    record["index"],
                    record.get("name", f"record {record['index']}"),
                    inchikey,
                    inchikey.split("-")[0] if inchikey else None,
                    record.get("smiles"),
                    record.get("formula", "C15H10O6"),
                    record.get("ontology", "Flavonoid"),
                    record.get("precursor_mz", 287.0550),
                    record.get("precursor_type", "[M+H]+"),
                    record.get("ion_mode", "Positive"),
                    record.get("instrument_type", "LC-ESI-QTOF"),
                    record.get("instrument_class", "TOF"),
                    record.get("collision_energy_raw", "20"),
                    record.get("collision_energy_value", 20.0),
                    record.get("rt_min"),
                    record.get("ccs"),
                    len(record["peaks"]),
                    payload_sha,
                ),
            )
        connection.commit()
    finally:
        connection.close()
    return database


def near_identical_pair() -> list[dict]:
    """Two entries that no spectrum can separate: same formula, same condition, different skeleton."""
    return [
        {
            "index": 1, "name": "Kaempferol", "inchikey": "IYRMWMYZSQPJKC-UHFFFAOYSA-N",
            "formula": "C15H10O6", "peaks": spectrum(INTENSITIES_A), "rt_min": 5.10, "ccs": 168.2,
        },
        {
            "index": 2, "name": "Luteolin", "inchikey": "IQPNAANSBPBGFQ-UHFFFAOYSA-N",
            "formula": "C15H10O6", "peaks": spectrum(INTENSITIES_A_NEAR), "rt_min": 5.90, "ccs": 171.9,
        },
    ]


def counts(database) -> dict[str, int]:
    connection = sqlite3.connect(database)
    try:
        return {
            "similarity": connection.execute("SELECT COUNT(*) FROM spectrum_similarity").fetchone()[0],
            "classes": connection.execute("SELECT COUNT(*) FROM ambiguity_class").fetchone()[0],
            "members": connection.execute("SELECT COUNT(*) FROM ambiguity_class_member").fetchone()[0],
        }
    finally:
        connection.close()


class DefinitionTests(unittest.TestCase):
    def test_report_has_the_house_shape(self):
        report = AmbiguityReport()
        self.assertTrue(report.valid)
        self.assertEqual((report.blocks, report.edges, report.classes, report.singletons), (0, 0, 0, 0))
        report.errors.append("boom")
        self.assertFalse(report.valid)

    def test_as_rules_round_trips_every_parameter(self):
        # The threshold is a convention with no attached error rate, so a class is only recalibratable if
        # every parameter that produced it survives the round trip.
        definition = ClassDefinition(
            definition_id="ambiguity-v9",
            weighted_cosine_threshold=0.80,
            entropy_similarity_threshold=0.75,
            mz_tolerance_da=0.05,
            minimum_informative_peaks=8,
            minimum_matched_peaks=5,
            relative_floor=0.02,
            precursor_tolerance_da=0.02,
            precursor_tolerance_ppm=5.0,
            require_formula_agreement=False,
            require_condition_match=False,
            linkage="anchored_neighbourhood",
        )
        rules = definition.as_rules()
        self.assertEqual(ClassDefinition(**rules), definition)
        self.assertEqual(
            sorted(rules),
            [
                "definition_id", "entropy_similarity_threshold", "linkage", "minimum_informative_peaks",
                "minimum_matched_peaks", "mz_tolerance_da", "precursor_tolerance_da",
                "precursor_tolerance_ppm", "relative_floor", "require_condition_match",
                "require_formula_agreement", "weighted_cosine_threshold",
            ],
        )

    def test_default_definition_states_the_provisional_convention(self):
        definition = ClassDefinition()
        self.assertEqual(definition.definition_id, "ambiguity-v1")
        self.assertAlmostEqual(definition.weighted_cosine_threshold, 0.90)
        self.assertAlmostEqual(definition.entropy_similarity_threshold, 0.85)
        self.assertAlmostEqual(definition.mz_tolerance_da, 0.025)
        self.assertEqual(definition.linkage, "anchored_neighbourhood")

    def test_collision_energy_is_binned_and_missing_values_stay_unknown(self):
        self.assertEqual(collision_energy_bin(20.0), "20-30")
        self.assertEqual(collision_energy_bin(29.9), "20-30")
        self.assertEqual(collision_energy_bin(45.0), "40-50")
        # A quarter of the public positive records carry an unusable collision energy; the condition is
        # then unestablished rather than equal to anything.
        self.assertEqual(collision_energy_bin(None), "unknown")

    def test_blocking_key_names_mode_precursor_type_and_window(self):
        row = {"ion_mode": "Positive", "precursor_type": "[M+H]+", "precursor_mz": 287.0550}
        self.assertEqual(blocking_key(row), "positive|[M+H]+|28705")
        self.assertEqual(
            blocking_key({"ion_mode": None, "precursor_type": None, "precursor_mz": None}),
            "unknown|unknown|unknown",
        )


class TargetCaseTests(unittest.TestCase):
    def test_indistinguishable_isomers_become_a_class(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), near_identical_pair())
            report = compute_ambiguity_classes(database)
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.blocks, 1)
            self.assertEqual(report.pairs_compared, 1)
            self.assertEqual(report.edges, 1)
            self.assertEqual(report.pairs_insufficient_evidence, 0)
            self.assertEqual(report.pairs_condition_mismatch, 0)
            self.assertEqual(report.pairs_isobaric_not_isomeric, 0)
            self.assertEqual(report.singletons, 0)
            # One class per anchor: the neighbourhood is defined relative to an entry, not carved out of
            # a partition, so the symmetric pair yields N(1) and N(2) with the same two members.
            self.assertEqual(report.classes, 2)

            first = ambiguity_class_for(database, reference_id(1))
            self.assertIsNotNone(first)
            self.assertEqual(first["member_count"], 2)
            self.assertEqual(first["formula_agreement"], "same_formula")
            self.assertEqual(first["skeleton_agreement"], "different_skeleton")
            self.assertEqual(first["linkage_rule"], "anchored_neighbourhood")
            self.assertEqual(first["score_convention"], "cosine")
            self.assertTrue(first["clique"])
            self.assertGreaterEqual(first["min_pairwise_score"], 0.90)
            self.assertEqual(first["blocking_key"], "positive|[M+H]+|28705")
            self.assertEqual(first["condition_scope"]["instrument_classes"], ["tof"])
            self.assertEqual(first["condition_scope"]["collision_energy_bins"], ["20-30"])
            self.assertTrue(first["condition_scope"]["condition_established"])
            # The rule set travels with the class because the DDL has no column for it and a threshold
            # that cannot be recovered cannot be recalibrated.
            self.assertEqual(
                first["condition_scope"]["definition_rules"], ClassDefinition().as_rules()
            )
            self.assertEqual(
                first["library_scope"]["completeness"], "systematically_incomplete"
            )
            self.assertEqual(
                [library["library_id"] for library in first["library_scope"]["compared_libraries"]],
                [LIBRARY_ID],
            )
            self.assertEqual(first["library_scope"]["member_libraries"], [LIBRARY_ID])
            # Both members carry a retention time and a CCS that differ, so both would break the tie.
            self.assertEqual(first["discriminating_evidence_needed"], "RT,IM,RS")
            self.assertEqual(
                [member["reference_spectrum_id"] for member in first["members"]],
                [reference_id(1), reference_id(2)],
            )
            self.assertEqual(
                [member["compound_name"] for member in first["members"]], ["Kaempferol", "Luteolin"]
            )

            second = ambiguity_class_for(database, reference_id(2))
            self.assertEqual(second["members"][0]["reference_spectrum_id"], reference_id(2))
            self.assertNotEqual(second["ambiguity_class_id"], first["ambiguity_class_id"])

    def test_same_skeleton_members_are_reported_as_a_stereochemical_question(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["inchikey"] = "IYRMWMYZSQPJKC-ZZZZZZZZZZ-N"
            database = build_library(Path(directory), records)
            compute_ambiguity_classes(database)
            found = ambiguity_class_for(database, reference_id(1))
            # The 2D structure is determined; only stereochemistry or labelling is open, which is a
            # different outcome from a genuine constitutional-isomer class.
            self.assertEqual(found["skeleton_agreement"], "same_skeleton")
            self.assertEqual(found["formula_agreement"], "same_formula")


class RefusalTests(unittest.TestCase):
    def test_differing_formula_is_isobaric_not_ambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["formula"] = "C11H14N4O4"
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(database)
            self.assertEqual(report.pairs_compared, 1)
            # Resolvable by formula and mass evidence alone, so an artifact of the window width.
            self.assertEqual(report.pairs_isobaric_not_isomeric, 1)
            self.assertEqual(report.edges, 0)
            self.assertEqual(report.classes, 0)
            self.assertEqual(report.singletons, 2)
            self.assertEqual(counts(database), {"similarity": 0, "classes": 0, "members": 0})
            self.assertIsNone(ambiguity_class_for(database, reference_id(1)))

    def test_two_peak_spectra_are_insufficient_evidence_not_indistinguishable(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[0]["peaks"] = [[100.1, 1.0], [200.1, 0.5]]
            records[1]["peaks"] = [[100.1, 1.0], [200.1, 0.5]]
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(database)
            self.assertEqual(report.pairs_compared, 1)
            self.assertEqual(report.pairs_insufficient_evidence, 1)
            # Two spectra with two peaks each look identical without that meaning anything, so the pair
            # must never be counted, stored or reported as indistinguishable.
            self.assertEqual(report.edges, 0)
            self.assertEqual(report.classes, 0)
            self.assertEqual(report.pairs_condition_mismatch, 0)
            self.assertEqual(report.pairs_isobaric_not_isomeric, 0)
            self.assertEqual(counts(database), {"similarity": 0, "classes": 0, "members": 0})

    def test_different_collision_energy_bins_are_a_condition_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["collision_energy_value"] = 45.0
            records[1]["collision_energy_raw"] = "45HCD"
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(database)
            self.assertEqual(report.pairs_compared, 1)
            # A third outcome: neither distinguishable nor indistinguishable, because the two records
            # were not measured under the same condition.
            self.assertEqual(report.pairs_condition_mismatch, 1)
            self.assertEqual(report.pairs_insufficient_evidence, 0)
            self.assertEqual(report.edges, 0)
            self.assertEqual(report.classes, 0)

    def test_different_instrument_class_is_a_condition_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["instrument_class"] = "FT"
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(database)
            self.assertEqual(report.pairs_condition_mismatch, 1)
            self.assertEqual(report.edges, 0)

    def test_condition_mismatch_can_be_relaxed_by_the_definition(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["collision_energy_value"] = 45.0
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(
                database, definition=ClassDefinition(require_condition_match=False)
            )
            self.assertEqual(report.pairs_condition_mismatch, 0)
            self.assertEqual(report.edges, 1)
            found = ambiguity_class_for(database, reference_id(1))
            self.assertEqual(found["condition_scope"]["collision_energy_bins"], ["20-30", "40-50"])
            self.assertFalse(found["condition_scope"]["require_condition_match"])


class BlockingTests(unittest.TestCase):
    def test_a_precursor_outside_the_window_is_never_compared(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["precursor_mz"] = 287.4550
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(database)
            # Separable by MS1 mass alone, so the pair never reaches the similarity step at all.
            self.assertEqual(report.blocks, 0)
            self.assertEqual(report.pairs_compared, 0)
            self.assertEqual(report.edges, 0)
            self.assertEqual(report.singletons, 2)

    def test_different_precursor_type_or_ion_mode_never_share_a_block(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["precursor_type"] = "[M+Na]+"
            database = build_library(Path(directory), records)
            self.assertEqual(report_pairs(database), 0)
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["ion_mode"] = "Negative"
            database = build_library(Path(directory), records)
            self.assertEqual(report_pairs(database), 0)

    def test_blocks_are_reported_with_their_size(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), non_transitive_triple())
            observed: list[dict] = []
            compute_ambiguity_classes(database, progress=observed.append)
            self.assertEqual([event["block_size"] for event in observed], [3])
            self.assertEqual(observed[0]["blocking_key"], "positive|[M+H]+|28705")

    def test_iter_blocks_yields_each_membership_once(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), non_transitive_triple())
            blocks = list(iter_blocks(database, ClassDefinition()))
            self.assertEqual(len(blocks), 1)
            self.assertEqual(
                sorted(row["reference_spectrum_id"] for row in blocks[0]),
                [reference_id(1), reference_id(2), reference_id(3)],
            )
            self.assertEqual(len(blocks[0][0]["peaks"]), len(MZ))


def report_pairs(database) -> int:
    return compute_ambiguity_classes(database).pairs_compared


def non_transitive_triple() -> list[dict]:
    return [
        {"index": 1, "name": "A", "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N", "peaks": spectrum(INTENSITIES_A)},
        {"index": 2, "name": "B", "inchikey": "BBBBBBBBBBBBBB-UHFFFAOYSA-N", "peaks": spectrum(INTENSITIES_B)},
        {"index": 3, "name": "C", "inchikey": "CCCCCCCCCCCCCC-UHFFFAOYSA-N", "peaks": spectrum(INTENSITIES_C)},
    ]


class NonTransitivityTests(unittest.TestCase):
    def test_neighbourhoods_are_anchored_and_no_closure_is_taken(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), non_transitive_triple())
            report = compute_ambiguity_classes(database)
            self.assertTrue(report.valid, report.errors)
            self.assertEqual(report.pairs_compared, 3)
            self.assertEqual(report.edges, 2)
            self.assertEqual(report.classes, 3)

            middle = ambiguity_class_for(database, reference_id(2))
            self.assertEqual(
                [member["reference_spectrum_id"] for member in middle["members"]],
                [reference_id(2), reference_id(1), reference_id(3)],
            )
            # A and C are both indistinguishable from B while being distinguishable from each other, so
            # the honest statement is weaker than a clique and is recorded as such.
            self.assertFalse(middle["clique"])
            self.assertEqual(middle["member_count"], 3)

            first = ambiguity_class_for(database, reference_id(1))
            self.assertEqual(
                [member["reference_spectrum_id"] for member in first["members"]],
                [reference_id(1), reference_id(2)],
            )
            # The transitive closure would have pulled C into N(A). Similarity is not transitive, so it
            # must not, and A's class stays a statement about A.
            self.assertNotIn(
                reference_id(3), [member["reference_spectrum_id"] for member in first["members"]]
            )
            self.assertTrue(first["clique"])

            last = ambiguity_class_for(database, reference_id(3))
            self.assertEqual(
                [member["reference_spectrum_id"] for member in last["members"]],
                [reference_id(3), reference_id(2)],
            )
            self.assertTrue(last["clique"])

    def test_is_clique_reads_an_edge_list_or_an_edge_map(self):
        members = ["a", "b", "c"]
        edges = {("a", "b"): {}, ("b", "c"): {}}
        self.assertFalse(is_clique(members, edges))
        self.assertTrue(is_clique(members, {**edges, ("a", "c"): {}}))
        self.assertTrue(is_clique(["a", "b"], [("b", "a")]))
        self.assertTrue(is_clique(["a"], []))


class StorageTests(unittest.TestCase):
    def test_edges_are_stored_once_in_canonical_order(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), near_identical_pair())
            compute_ambiguity_classes(database)
            connection = sqlite3.connect(database)
            connection.row_factory = sqlite3.Row
            try:
                rows = connection.execute("SELECT * FROM spectrum_similarity").fetchall()
                self.assertEqual(len(rows), 1)
                row = rows[0]
                self.assertEqual(row["subject_kind_a"], "reference_spectrum")
                self.assertEqual(row["subject_kind_b"], "reference_spectrum")
                # The relation is symmetric, so the pair is stored once with its ids sorted as strings.
                self.assertLess(row["subject_id_a"], row["subject_id_b"])
                self.assertEqual(row["method"], "ambiguity_weighted_cosine")
                # The label plus a digest of the rules it stood for. The label alone was the stored
                # version, and it does not change when a threshold does.
                self.assertEqual(
                    row["method_version"], f"ambiguity-v1+{ClassDefinition().rules_sha256}"
                )
                self.assertEqual(row["score_convention"], "cosine")
                self.assertEqual(row["secondary_method"], "entropy_similarity")
                self.assertGreaterEqual(row["score"], 0.90)
                self.assertGreaterEqual(row["secondary_score"], 0.85)
                self.assertEqual(row["matched_peak_count"], len(MZ))
                self.assertAlmostEqual(row["mz_tolerance_da"], 0.025)
            finally:
                connection.close()

    def test_recomputing_is_idempotent(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), non_transitive_triple())
            first = compute_ambiguity_classes(database)
            before = counts(database)
            second = compute_ambiguity_classes(database)
            self.assertEqual(
                (second.blocks, second.pairs_compared, second.edges, second.classes),
                (first.blocks, first.pairs_compared, first.edges, first.classes),
            )
            self.assertEqual(counts(database), before)
            self.assertEqual(before, {"similarity": 2, "classes": 3, "members": 7})

    def test_a_stricter_definition_does_not_overwrite_the_looser_one(self):
        # The identity of a class includes the rules it was computed under, so a second run at a
        # different threshold is a different answer to a different question rather than a replacement
        # of the first. It used to be a silent replacement: the identifier encoded only definition_id,
        # which defaults to "ambiguity-v1" whatever thresholds the command line carried, so the same
        # rows came back holding a different meaning with nothing stored that said so. An annotation
        # citing a class would have become a claim about rules it never saw.
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), near_identical_pair())
            compute_ambiguity_classes(database)
            before = counts(database)
            self.assertEqual(before, {"similarity": 1, "classes": 2, "members": 4})

            report = compute_ambiguity_classes(
                database, definition=ClassDefinition(weighted_cosine_threshold=0.999999)
            )

            self.assertEqual(report.edges, 0)
            self.assertEqual(report.classes, 0)
            # The stricter run admitted nothing, and the looser run's rows are still exactly as they
            # were, still resolvable by anything that cited them.
            self.assertEqual(counts(database), before)

    def test_the_rules_digest_changes_with_any_threshold(self):
        base = ClassDefinition()
        self.assertNotEqual(base.rules_sha256, ClassDefinition(weighted_cosine_threshold=0.95).rules_sha256)
        self.assertNotEqual(base.rules_sha256, ClassDefinition(entropy_similarity_threshold=0.5).rules_sha256)
        self.assertNotEqual(base.rules_sha256, ClassDefinition(mz_tolerance_da=0.05).rules_sha256)
        # The same rules under a different label are the same rules.
        self.assertEqual(base.rules_sha256, ClassDefinition().rules_sha256)

    def test_the_tool_run_is_recorded_on_every_row_it_produced(self):
        with tempfile.TemporaryDirectory() as directory:
            database = build_library(Path(directory), near_identical_pair())
            connection = connect(database)
            try:
                connection.execute(
                    "INSERT INTO annotation_tool_run(tool_run_id, tool_name, tool_version) "
                    "VALUES (?, ?, ?)",
                    ("urn:msdial:tool-run:ambiguity:test", "msdial_spectrum_catalog.ambiguity", "0.1.0"),
                )
                connection.commit()
            finally:
                connection.close()
            compute_ambiguity_classes(database, tool_run_id="urn:msdial:tool-run:ambiguity:test")
            connection = sqlite3.connect(database)
            try:
                self.assertEqual(
                    connection.execute(
                        "SELECT COUNT(*) FROM spectrum_similarity WHERE tool_run_id = ?",
                        ("urn:msdial:tool-run:ambiguity:test",),
                    ).fetchone()[0],
                    1,
                )
                self.assertEqual(
                    connection.execute(
                        "SELECT COUNT(*) FROM ambiguity_class WHERE tool_run_id = ?",
                        ("urn:msdial:tool-run:ambiguity:test",),
                    ).fetchone()[0],
                    2,
                )
            finally:
                connection.close()

    def test_lookup_returns_none_for_an_entry_with_no_class(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["precursor_mz"] = 287.4550
            database = build_library(Path(directory), records)
            compute_ambiguity_classes(database)
            self.assertIsNone(ambiguity_class_for(database, reference_id(1)))
            self.assertIsNone(ambiguity_class_for(database, "urn:msdial:reference-spectrum:absent"))


class DiscriminationTests(unittest.TestCase):
    def test_discriminating_ions_report_the_ion_one_member_lacks(self):
        members = [
            {
                "reference_spectrum_id": "a",
                "peaks": [(100.1, 1.0), (150.1, 0.5), (275.1, 0.4)],
            },
            {
                "reference_spectrum_id": "b",
                "peaks": [(100.1, 1.0), (150.1, 0.5)],
            },
        ]
        ions = discriminating_ions(members, tolerance=0.025)
        self.assertEqual(len(ions), 1)
        self.assertAlmostEqual(ions[0]["mz"], 275.1)
        self.assertEqual(ions[0]["present_in"], ["a"])
        self.assertEqual(ions[0]["absent_from"], ["b"])
        self.assertAlmostEqual(ions[0]["mean_relative_intensity"], 0.4)
        self.assertAlmostEqual(ions[0]["split_fraction"], 0.5)

    def test_the_most_evenly_splitting_ion_comes_first(self):
        members = [
            {"reference_spectrum_id": "a", "peaks": [(100.1, 1.0), (200.1, 0.6), (300.1, 0.5)]},
            {"reference_spectrum_id": "b", "peaks": [(100.1, 1.0), (200.1, 0.6)]},
            {"reference_spectrum_id": "c", "peaks": [(100.1, 1.0), (200.1, 0.6)]},
            {"reference_spectrum_id": "d", "peaks": [(100.1, 1.0), (400.1, 0.5)]},
        ]
        ions = discriminating_ions(members, tolerance=0.025)
        # 200.1 splits three against one; 300.1 and 400.1 split one against three. The most even split
        # is the most informative single measurement, so it leads.
        self.assertAlmostEqual(ions[0]["mz"], 200.1)
        self.assertEqual(ions[0]["present_in"], ["a", "b", "c"])
        self.assertNotIn(100.1, [ion["mz"] for ion in ions])

    def test_evidence_needed_is_derived_from_what_the_rows_carry(self):
        # RS always: an authentic standard resolves any class.
        self.assertEqual(discriminating_evidence_needed([{"rt_min": None}, {"rt_min": None}]), "RS")
        self.assertEqual(
            discriminating_evidence_needed([{"rt_min": 5.1, "ccs": 168.0}, {"rt_min": 5.9, "ccs": 168.2}]),
            "RT,RS",
        )
        self.assertEqual(
            discriminating_evidence_needed([{"rt_min": 5.1, "ccs": 168.0}, {"rt_min": 5.11, "ccs": 180.0}]),
            "IM,RS",
        )
        # One member with a retention time cannot separate anything, so RT is not claimed.
        self.assertEqual(discriminating_evidence_needed([{"rt_min": 5.1}, {"rt_min": None}]), "RS")

    def test_a_class_carries_the_ions_that_would_separate_it(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["peaks"] = records[1]["peaks"] + [[275.1, 0.05]]
            database = build_library(Path(directory), records)
            compute_ambiguity_classes(database)
            found = ambiguity_class_for(database, reference_id(1))
            self.assertEqual(
                [ion["mz"] for ion in found["discriminating_mz"]], [275.1]
            )
            self.assertEqual(found["discriminating_mz"][0]["present_in"], [reference_id(2)])


class UnestablishedConditionTests(unittest.TestCase):
    def test_an_unusable_collision_energy_weakens_the_class_rather_than_blocking_it(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            for record in records:
                record["collision_energy_raw"] = "nan"
                record["collision_energy_value"] = None
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(database)
            self.assertEqual(report.edges, 1)
            self.assertTrue(
                any("unestablished collision-energy bin" in warning for warning in report.warnings),
                report.warnings,
            )
            found = ambiguity_class_for(database, reference_id(1))
            self.assertEqual(found["condition_scope"]["collision_energy_bins"], ["unknown"])
            self.assertFalse(found["condition_scope"]["condition_established"])

    def test_a_missing_formula_is_refused_and_reported_separately(self):
        with tempfile.TemporaryDirectory() as directory:
            records = near_identical_pair()
            records[1]["formula"] = None
            database = build_library(Path(directory), records)
            report = compute_ambiguity_classes(database)
            self.assertEqual(report.edges, 0)
            self.assertEqual(report.pairs_isobaric_not_isomeric, 1)
            self.assertTrue(
                any("a FORMULA is missing" in warning for warning in report.warnings), report.warnings
            )


if __name__ == "__main__":
    unittest.main()
