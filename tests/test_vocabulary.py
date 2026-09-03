import json
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.vocabulary import (
    AXES,
    CANONICAL_EVIDENCE_ORDER,
    CLAIM_AXIS,
    DEFAULT_VOCABULARY,
    EVIDENCE_AXIS,
    LEVEL_AXIS,
    REGISTRY_DIR,
    AmbiguousClaimTokenError,
    ClaimReading,
    NotationError,
    available_versions,
    emit_notation,
    find_migration,
    load_use_cases,
    load_vocabulary,
    migrate_reading,
    parse_notation,
    parse_notation_any,
    resolve_claim_token,
    validate_combination,
)

V1 = "smb-v1-draft"
V2 = "smb-v2-consensus"

STRUCTURE = "smb:claim/structure"
SUBSTRUCTURE_COMPLETE = "smb:claim/substructure_complete"
SUBSTRUCTURE_INCOMPLETE = "smb:claim/substructure_incomplete"
CLASS = "smb:claim/class"
CONTEXTUAL = "smb:evidence/contextual"
FORMULA_MASS = "smb:evidence/formula_mass"
SPECTRAL_SIMILARITY = "smb:evidence/spectral_similarity"


class RegistryTests(unittest.TestCase):
    def test_both_versions_are_registered_and_default_is_the_consensus(self):
        self.assertEqual(available_versions(), (V1, V2))
        self.assertEqual(DEFAULT_VOCABULARY, V2)
        self.assertEqual(load_vocabulary(V2).status, "accepted")
        self.assertEqual(load_vocabulary(V1).status, "superseded")

    def test_registry_files_are_data_and_carry_the_required_keys(self):
        for version in available_versions():
            payload = json.loads((REGISTRY_DIR / f"{version}.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["version"], version)
            self.assertIn("schema_version", payload)
            self.assertIn("status", payload)
            self.assertEqual(set(payload["axes"]), set(AXES))
            for key in ("source_file", "source_date", "sections", "notes"):
                self.assertIn(key, payload["provenance"])

    def test_every_term_round_trips_token_to_concept_to_token(self):
        for version in available_versions():
            vocabulary = load_vocabulary(version)
            for axis in AXES:
                self.assertTrue(vocabulary.axes[axis])
                for term in vocabulary.axes[axis]:
                    by_token = vocabulary.term(axis, term.token)
                    by_concept = vocabulary.by_concept(axis, by_token.concept_id)
                    self.assertEqual(by_concept.token, term.token)
                    self.assertTrue(term.definition)
                    self.assertTrue(term.definition_source)

    def test_accepted_token_sets_match_the_confirmed_consensus(self):
        v2 = load_vocabulary(V2)
        self.assertEqual([term.token for term in v2.axes[LEVEL_AXIS]], ["L1", "L2", "L3", "L4", "L5"])
        self.assertEqual([term.token for term in v2.axes[CLAIM_AXIS]], ["SP", "SC", "SI", "CP"])
        self.assertEqual(
            [term.token for term in v2.axes[EVIDENCE_AXIS]],
            ["RS", "FM", "SL", "DF", "RT", "IM", "IS", "MN", "HO", "CO", "UN", "OS"],
        )
        self.assertTrue(all(term.status == "accepted" for term in v2.axes[EVIDENCE_AXIS]))
        v1 = load_vocabulary(V1)
        self.assertEqual([term.token for term in v1.axes[CLAIM_AXIS]], ["SP", "CP", "MO", "CL"])
        self.assertEqual(
            {term.token for term in v1.axes[EVIDENCE_AXIS] if term.status == "proposed"}, {"SM", "HO", "UN"}
        )

    def test_canonical_evidence_order_is_the_section_four_table_order(self):
        self.assertEqual(
            CANONICAL_EVIDENCE_ORDER,
            (
                "smb:evidence/reference_standard",
                "smb:evidence/formula_mass",
                "smb:evidence/spectral_library",
                "smb:evidence/diagnostic_fragment",
                "smb:evidence/retention",
                "smb:evidence/ion_mobility",
                "smb:evidence/in_silico",
                "smb:evidence/molecular_network",
                "smb:evidence/homologue",
                "smb:evidence/contextual",
                "smb:evidence/unclassified",
                "smb:evidence/other_spectroscopy",
            ),
        )

    def test_specificity_rank_is_declared_and_never_a_confidence_ranking(self):
        v2 = load_vocabulary(V2)
        self.assertEqual([term.specificity_rank for term in v2.axes[CLAIM_AXIS]], [1, 2, 3, 4])
        self.assertTrue(all(term.specificity_rank is None for term in v2.axes[LEVEL_AXIS]))
        self.assertTrue(any("not a confidence ranking" in rule["rule"] for rule in v2.decision_rules))

    def test_open_issues_are_machine_visible(self):
        issues = {issue["issue_id"] for issue in load_vocabulary(V2).open_issues}
        self.assertEqual(
            issues,
            {
                "ISSUE-SECTIONS-5-8-STALE",
                "ISSUE-CP-COLLISION",
                "ISSUE-SM-DROPPED",
                "ISSUE-LEVEL-NOT-DERIVABLE",
                "ISSUE-EVIDENCE-DELIMITER",
            },
        )
        self.assertTrue(all(issue["status"] == "open" for issue in load_vocabulary(V2).open_issues))

    def test_unknown_version_is_refused(self):
        with self.assertRaises(ValueError):
            load_vocabulary("smb-v3-imaginary")

    def test_registry_is_loaded_relative_to_the_package_not_the_working_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertTrue(Path(directory).is_dir())
            self.assertTrue((REGISTRY_DIR / f"{V2}.json").is_file())
            self.assertEqual(load_vocabulary(V2).version, V2)


class ParseTests(unittest.TestCase):
    def test_consensus_example_parses_and_re_emits_byte_identically(self):
        reading = parse_notation("L3-SC[FM,DF,MN,CO]", V2)
        self.assertEqual(reading.level, "L3")
        self.assertEqual(reading.claim_concept_ids, (SUBSTRUCTURE_COMPLETE,))
        self.assertEqual(reading.unknown_tokens, ())
        self.assertFalse(reading.unresolved)
        self.assertEqual(emit_notation(reading), "L3-SC[FM,DF,MN,CO]")

    def test_documented_grammar_examples_all_parse(self):
        for notation in ("L3-SC[FM,DF,MN,CO]", "L1[RS,RT,SL,FM]", "L4[FM]", "L5[]", "L5", "L3-SC+CP[DF,FM]"):
            self.assertEqual(parse_notation(notation, V2).notation_verbatim, notation)

    def test_level_five_parses_with_and_without_the_bracket_and_always_emits_it(self):
        for notation in ("L5", "L5[]"):
            reading = parse_notation(notation, V2)
            self.assertEqual(reading.evidence_concept_ids, ())
            self.assertEqual(emit_notation(reading), "L5[]")

    def test_claim_tag_is_only_legal_on_level_three(self):
        with self.assertRaises(NotationError):
            parse_notation("L2-SP[SL]", V2)
        with self.assertRaises(NotationError):
            parse_notation("L1-SP[RS]", V2)

    def test_slash_alternation_is_rejected(self):
        with self.assertRaises(NotationError) as caught:
            parse_notation("L3-SC/SI[FM,DF]", V2)
        self.assertIn("slash", str(caught.exception))

    def test_or_alternation_is_rejected(self):
        with self.assertRaises(NotationError) as caught:
            parse_notation("L3-SP[FM,DF] or L3-SC[FM,DF]", V2)
        self.assertIn("'or'", str(caught.exception))

    def test_prose_conditional_is_rejected(self):
        with self.assertRaises(NotationError) as caught:
            parse_notation("L2[SL,FM] if unambiguous; otherwise L3-SP[SL,FM]", V2)
        self.assertIn("prose conditional", str(caught.exception))

    def test_malformed_strings_are_rejected(self):
        for notation in ("", "   ", "L6[FM]", "L3-[FM]", "L3-SC[FM,]", "L3-SC[FM,FM]", "L3-SC+SC[FM]", "L3-sc[FM]"):
            with self.assertRaises(NotationError):
                parse_notation(notation, V2)

    def test_vocab_version_is_required(self):
        with self.assertRaises(NotationError):
            parse_notation("L3-CP[FM]", "")
        with self.assertRaises(ValueError):
            parse_notation("L3-CP[FM]", "not-a-version")

    def test_unknown_mode_is_refused(self):
        with self.assertRaises(ValueError):
            parse_notation("L4[FM]", V2, mode="lenient")


class CollisionTests(unittest.TestCase):
    def test_the_same_string_means_different_things_in_each_version(self):
        as_draft = parse_notation("L3-CP[FM]", V1)
        as_consensus = parse_notation("L3-CP[FM]", V2)
        self.assertEqual(as_draft.claim_concept_ids, (SUBSTRUCTURE_COMPLETE,))
        self.assertEqual(as_consensus.claim_concept_ids, (CLASS,))
        self.assertNotEqual(as_draft.claim_concept_ids, as_consensus.claim_concept_ids)
        self.assertEqual(emit_notation(as_draft), "L3-CP[FM]")
        self.assertEqual(emit_notation(as_consensus), "L3-CP[FM]")

    def test_parse_notation_any_returns_one_reading_per_version(self):
        readings = parse_notation_any("L3-CP[FM]")
        self.assertEqual(len(readings), 2)
        self.assertEqual({reading.vocab_version for reading in readings}, {V1, V2})
        self.assertEqual(
            {reading.vocab_version: reading.claim_concept_ids for reading in readings},
            {V1: (SUBSTRUCTURE_COMPLETE,), V2: (CLASS,)},
        )

    def test_resolving_a_rebound_token_without_a_version_is_an_error(self):
        with self.assertRaises(AmbiguousClaimTokenError) as caught:
            resolve_claim_token("CP")
        self.assertEqual(caught.exception.token, "CP")
        self.assertEqual(
            set(caught.exception.candidates), {(V1, SUBSTRUCTURE_COMPLETE), (V2, CLASS)}
        )

    def test_a_stable_token_resolves_without_a_version(self):
        self.assertEqual(resolve_claim_token("SP"), (V1, STRUCTURE))

    def test_parse_notation_any_returns_nothing_when_generations_are_mixed(self):
        self.assertEqual(parse_notation_any("L3-SC[FM,DF,MN,CX]"), ())


class EmitTests(unittest.TestCase):
    def test_evidence_is_emitted_in_canonical_order_regardless_of_input_order(self):
        scrambled = ClaimReading(
            level="L3",
            claim_concept_ids=(SUBSTRUCTURE_COMPLETE,),
            evidence_concept_ids=(
                CONTEXTUAL,
                "smb:evidence/molecular_network",
                FORMULA_MASS,
                "smb:evidence/diagnostic_fragment",
            ),
            vocab_version=V2,
            notation_verbatim="",
        )
        self.assertEqual(emit_notation(scrambled), "L3-SC[FM,DF,MN,CO]")

    def test_claims_are_emitted_in_specificity_order(self):
        reading = parse_notation("L3-CP+SC[DF,FM]", V2)
        self.assertEqual(emit_notation(reading), "L3-SC+CP[FM,DF]")

    def test_delimiter_is_configurable_but_never_the_mztab_separator(self):
        reading = parse_notation("L1[RS,RT,SL,FM]", V2)
        self.assertEqual(emit_notation(reading, delimiter=";"), "L1[RS,FM,SL,RT]".replace(",", ";"))
        with self.assertRaises(ValueError):
            emit_notation(reading, delimiter="|")

    def test_emitting_a_concept_the_target_version_does_not_know_is_refused(self):
        reading = ClaimReading(
            level="L3",
            claim_concept_ids=(SUBSTRUCTURE_COMPLETE,),
            evidence_concept_ids=(SPECTRAL_SIMILARITY,),
            vocab_version=V2,
            notation_verbatim="",
        )
        with self.assertRaises(NotationError):
            emit_notation(reading)


class MigrationTests(unittest.TestCase):
    def test_migration_entry_is_version_pair_scoped_and_declares_the_hazard(self):
        migration = find_migration(V1, V2)
        self.assertEqual(migration["map"][CLAIM_AXIS], {"SP": "SP", "CP": "SC", "MO": "SI", "CL": "CP"})
        self.assertEqual(migration["map"][EVIDENCE_AXIS], {"CX": "CO"})
        self.assertEqual(migration["reused_tokens"], ["CP"])
        self.assertEqual(migration["dropped_terms"], ["SM"])
        self.assertEqual(migration["confidence"], "high")
        self.assertFalse(migration["reviewed"])
        self.assertTrue(migration["evidence"])
        with self.assertRaises(ValueError):
            find_migration(V2, V1)

    def test_component_becomes_substructure_complete_and_contextual_is_renamed(self):
        drafted = parse_notation("L3-CP[FM,CX]", V1)
        migrated = migrate_reading(drafted, V2)
        self.assertEqual(migrated.vocab_version, V2)
        self.assertEqual(migrated.claim_concept_ids, (SUBSTRUCTURE_COMPLETE,))
        self.assertEqual(migrated.evidence_concept_ids, (FORMULA_MASS, CONTEXTUAL))
        self.assertFalse(migrated.unresolved)
        self.assertEqual(emit_notation(migrated), "L3-SC[FM,CO]")

    def test_motif_and_class_migrate_to_the_consensus_tokens(self):
        self.assertEqual(emit_notation(migrate_reading(parse_notation("L3-MO[DF,FM]", V1), V2)), "L3-SI[FM,DF]")
        self.assertEqual(emit_notation(migrate_reading(parse_notation("L3-CL[DF,FM]", V1), V2)), "L3-CP[FM,DF]")

    def test_spectral_similarity_is_dropped_and_marks_the_reading_unresolved(self):
        drafted = parse_notation("L3-SP[SM,FM]", V1, mode="permissive")
        self.assertEqual(drafted.evidence_concept_ids, (SPECTRAL_SIMILARITY, FORMULA_MASS))
        migrated = migrate_reading(drafted, V2)
        self.assertEqual(migrated.evidence_concept_ids, (FORMULA_MASS,))
        self.assertIn("SM", migrated.unknown_tokens)
        self.assertTrue(migrated.unresolved)
        self.assertEqual(emit_notation(migrated), "L3-SP[FM]")

    def test_migrating_to_the_same_version_is_an_identity(self):
        reading = parse_notation("L4[FM]", V2)
        self.assertIs(migrate_reading(reading, V2), reading)

    def test_emit_to_another_version_migrates_first(self):
        drafted = parse_notation("L3-CP[FM,CX]", V1)
        self.assertEqual(emit_notation(drafted, V2), "L3-SC[FM,CO]")


class ModeTests(unittest.TestCase):
    def test_strict_rejects_a_proposed_status_token(self):
        with self.assertRaises(NotationError) as caught:
            parse_notation("L3-SP[SM,FM]", V1)
        self.assertIn("proposed", str(caught.exception))

    def test_permissive_accepts_a_proposed_status_token(self):
        reading = parse_notation("L3-SP[SM,FM]", V1, mode="permissive")
        self.assertEqual(reading.evidence_concept_ids, (SPECTRAL_SIMILARITY, FORMULA_MASS))
        self.assertEqual(reading.unknown_tokens, ())
        self.assertFalse(reading.unresolved)

    def test_unknown_token_is_refused_recorded_then_quarantined(self):
        with self.assertRaises(NotationError):
            parse_notation("L3-SC[ZZ,FM]", V2)
        permissive = parse_notation("L3-SC[ZZ,FM]", V2, mode="permissive")
        self.assertEqual(permissive.unknown_tokens, ("ZZ",))
        self.assertEqual(permissive.evidence_concept_ids, (FORMULA_MASS,))
        self.assertFalse(permissive.unresolved)
        quarantined = parse_notation("L3-SC[ZZ,FM]", V2, mode="quarantine")
        self.assertEqual(quarantined.unknown_tokens, ("ZZ",))
        self.assertTrue(quarantined.unresolved)

    def test_quarantine_marks_even_a_clean_reading_as_unresolved(self):
        self.assertTrue(parse_notation("L4[FM]", V2, mode="quarantine").unresolved)

    def test_quarantine_accepts_any_status(self):
        reading = parse_notation("L3-SP[SM,FM]", V1, mode="quarantine")
        self.assertEqual(reading.evidence_concept_ids, (SPECTRAL_SIMILARITY, FORMULA_MASS))
        self.assertEqual(reading.unknown_tokens, ())


class CombinationTests(unittest.TestCase):
    def test_discouraged_and_avoided_combinations_warn(self):
        self.assertTrue(validate_combination(parse_notation("L3-SP+CP[FM]", V2)))
        self.assertIn("discouraged", validate_combination(parse_notation("L3-SP+CP[FM]", V2))[0])
        self.assertTrue(validate_combination(parse_notation("L3-SC+SI[FM]", V2)))
        self.assertIn("avoid", validate_combination(parse_notation("L3-SC+SI[FM]", V2))[0])

    def test_allowed_combinations_and_single_claims_are_silent(self):
        for notation in ("L3-SC+CP[DF,FM]", "L3-SI+CP[DF,FM]", "L3-SP[FM]", "L3-SC[FM]", "L4[FM]", "L5[]"):
            self.assertEqual(validate_combination(parse_notation(notation, V2)), ())

    def test_a_combination_the_proposal_does_not_describe_is_flagged(self):
        warnings = validate_combination(parse_notation("L3-SP+SI[FM]", V2))
        self.assertEqual(len(warnings), 1)
        self.assertIn("not described", warnings[0])

    def test_validate_combination_never_raises(self):
        broken = ClaimReading(
            level="L3",
            claim_concept_ids=("smb:claim/nonsense", "smb:claim/other"),
            evidence_concept_ids=(),
            vocab_version="smb-v9-missing",
            notation_verbatim="",
        )
        self.assertEqual(validate_combination(broken), ())


class UseCaseTests(unittest.TestCase):
    def setUp(self):
        self.use_cases = load_use_cases()

    def test_all_twelve_consensus_use_cases_are_recorded(self):
        cases = self.use_cases["cases"]
        self.assertEqual(len(cases), 12)
        self.assertEqual([case["case_id"] for case in cases], [f"case-{number}" for number in range(1, 13)])
        for case in cases:
            self.assertIn("v2-docx", [notation["source"] for notation in case["notations"]])
            self.assertIn("v2-docx", case["sources_in_agreement"])

    def test_recorded_notations_behave_exactly_as_the_flags_claim(self):
        recorded = [notation for case in self.use_cases["cases"] for notation in case["notations"]]
        recorded.extend(self.use_cases["unmatched_source_cases"])
        self.assertTrue(recorded)
        for notation in recorded:
            text = notation["notation_verbatim"]
            if notation["parseable"]:
                self.assertIsNone(notation["reason_outside_grammar"])
                readings = parse_notation_any(text)
                self.assertEqual([reading.vocab_version for reading in readings], notation["resolvable_under"])
            else:
                self.assertTrue(notation["reason_outside_grammar"])
                self.assertEqual(notation["resolvable_under"], [])
                with self.assertRaises(NotationError):
                    parse_notation_any(text)

    def test_the_out_of_grammar_forms_observed_in_the_sources_are_all_rejected(self):
        rejected = {
            notation["reason_outside_grammar"].split(":")[0]
            for case in self.use_cases["cases"]
            for notation in case["notations"]
            if not notation["parseable"]
        }
        self.assertEqual(rejected, {"or_alternation", "prose_conditional", "slash_alternation"})

    def test_level_cannot_be_derived_from_the_evidence_bracket(self):
        case = next(case for case in self.use_cases["cases"] if case["case_id"] == "case-7")
        notation = next(entry for entry in case["notations"] if entry["source"] == "v2-docx")
        self.assertEqual(notation["notation_verbatim"], "L2[SL,FM] or L3-SP[SL,FM]")
        first = parse_notation("L2[SL,FM]", V2)
        second = parse_notation("L3-SP[SL,FM]", V2)
        self.assertEqual(first.evidence_concept_ids, second.evidence_concept_ids)
        self.assertNotEqual(first.level, second.level)

    def test_divergent_cases_keep_every_source_string_verbatim(self):
        case = next(case for case in self.use_cases["cases"] if case["case_id"] == "case-3")
        strings = {notation["source"]: notation["notation_verbatim"] for notation in case["notations"]}
        self.assertEqual(strings["v2-docx"], "L3-SC[FM,DF,MN,CX]")
        self.assertEqual(strings["table-s9"], "L3-MO[FM,DF,MN,CX]")
        self.assertTrue(case["divergent"])


if __name__ == "__main__":
    unittest.main()
