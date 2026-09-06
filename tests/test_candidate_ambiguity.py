"""Why a candidate set holds more than one entry, and what each reason permits.

Several candidates is not one situation. The cascade separates the reasons because they have opposite
consequences: a rule-based lipid annotation is already at the resolution its fragments support, a
mass-window artifact is separable by accurate mass, a set with no product-ion comparison establishes
nothing either way, and only the last two -- the library cannot separate these, or the library can and
this run did not -- say anything about a claim ceiling.

The order of the cascade is itself a decision under test. A fact the run establishes on its own is
never displaced by the library being unavailable, and anything undecidable is not_assessed with the
reason rather than one of the decided states.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.candidate_ambiguity import (
    STATE_LIBRARY_INDISTINGUISHABLE,
    STATE_MASS_SEPARABLE,
    STATE_NO_MS2_EVIDENCE,
    STATE_NOT_ASSESSED,
    STATE_QUERY_NON_DISCRIMINATION,
    STATE_RULE_BASED_LIBRARY,
    STATE_SINGLE_CANDIDATE,
    AssessmentDefinition,
    assess_candidate_set,
    assess_run_candidates,
)


def _candidate(**overrides) -> dict:
    candidate = {
        "candidate_rank": 1,
        "is_spectrum_match": 1,
        "is_spectrum_comparison_performed": 1,
        "annotation_tag": "430",
        "database_id": "MspDB",
        "formula": "C15H10O6",
        "inchikey": "IYRMWMYZSQPJKC-UHFFFAOYSA-N",
        "reference_mz": 287.0550,
        "weighted_dot_product": 0.90,
        "name": "Kaempferol",
    }
    candidate.update(overrides)
    return candidate


def _tied_pair() -> list[dict]:
    """Two MS/MS matches the product-ion evidence does not separate."""
    return [
        _candidate(candidate_rank=1, weighted_dot_product=0.9000),
        _candidate(
            candidate_rank=2,
            weighted_dot_product=0.8995,
            inchikey="IQPNAANSBPBGFQ-UHFFFAOYSA-N",
            name="Luteolin",
        ),
    ]


class CascadeOrderTests(unittest.TestCase):
    def test_no_product_ion_comparison_decides_nothing(self):
        # Just over half the annotated features of the reference run land here: precursor-window
        # suggestions whose alignment consensus spectrum has no peaks at all.
        candidates = [
            _candidate(is_spectrum_match=0, is_spectrum_comparison_performed=0, annotation_tag="530"),
            _candidate(
                candidate_rank=2, is_spectrum_match=0, is_spectrum_comparison_performed=0,
                annotation_tag="530", inchikey="IQPNAANSBPBGFQ-UHFFFAOYSA-N",
            ),
        ]

        verdict = assess_candidate_set(candidates)

        self.assertEqual(STATE_NO_MS2_EVIDENCE, verdict["state"])
        self.assertEqual("unknown", verdict["ceiling_effect"])

    def test_a_rule_based_lipid_annotation_is_out_of_scope(self):
        # The LBM library's registered spectra gate the match; the name and its resolution come from
        # diagnostic-ion logic in MS-DIAL's source. Asking whether two library entries are spectrally
        # indistinguishable is not the question that applies.
        candidates = [
            _candidate(annotation_tag="410", database_id="LbmDB", name="PC 16:0_18:1"),
            _candidate(
                candidate_rank=2, annotation_tag="410", database_id="LbmDB",
                name="PC 18:1_16:0", inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            ),
        ]

        verdict = assess_candidate_set(candidates)

        self.assertEqual(STATE_RULE_BASED_LIBRARY, verdict["state"])
        self.assertEqual("no_effect", verdict["ceiling_effect"])

    def test_the_lipid_tag_is_enough_without_the_database_name(self):
        # The tag comes from flags only the rule-based path sets, so it is evidence from the code
        # rather than from a naming convention.
        verdict = assess_candidate_set(
            [
                _candidate(annotation_tag="420", database_id="SomeOtherDB"),
                _candidate(candidate_rank=2, annotation_tag="420", database_id="SomeOtherDB"),
            ]
        )
        self.assertEqual(STATE_RULE_BASED_LIBRARY, verdict["state"])

    def test_a_suggestion_beside_a_match_is_not_a_two_way_ambiguity(self):
        # Rank is a lexicographic order on (evidence tier, score): a real MS/MS match outranks a
        # higher-scoring precursor-only suggestion. Counting the suggestion as a contender would
        # overstate the ambiguity roughly fourfold on the reference run.
        candidates = [
            _candidate(candidate_rank=1),
            _candidate(
                candidate_rank=2, is_spectrum_match=0, annotation_tag="530",
                weighted_dot_product=0.0, inchikey="IQPNAANSBPBGFQ-UHFFFAOYSA-N",
            ),
        ]

        verdict = assess_candidate_set(candidates)

        self.assertEqual(STATE_SINGLE_CANDIDATE, verdict["state"])
        self.assertEqual(1, verdict["contender_count"])
        self.assertEqual("spectrum_matched", verdict["contender_selection_rule"])

    def test_contenders_separable_by_accurate_mass_are_a_window_artifact(self):
        candidates = [
            _candidate(reference_mz=287.0550),
            _candidate(
                candidate_rank=2, reference_mz=287.0610, formula="C16H14O5",
                inchikey="IQPNAANSBPBGFQ-UHFFFAOYSA-N",
            ),
        ]

        verdict = assess_candidate_set(candidates)

        self.assertEqual(STATE_MASS_SEPARABLE, verdict["state"])
        self.assertAlmostEqual(6.0, verdict["max_reference_mz_spread_mda"], places=3)

    def test_a_clear_spectral_winner_is_not_an_ambiguity(self):
        # Every contender was scored against the same query spectrum by the same annotator on the same
        # convention, so this within-set comparison is legitimate even though cross-feature ones are
        # not. On the reference run only 76 of 349 multi-match sets are genuinely tied.
        candidates = [
            _candidate(weighted_dot_product=0.95),
            _candidate(
                candidate_rank=2, weighted_dot_product=0.60,
                inchikey="IQPNAANSBPBGFQ-UHFFFAOYSA-N",
            ),
        ]

        verdict = assess_candidate_set(candidates)

        self.assertEqual(STATE_SINGLE_CANDIDATE, verdict["state"])
        self.assertAlmostEqual(0.35, verdict["discrimination_margin"], places=6)


class LibraryStageTests(unittest.TestCase):
    def test_without_a_library_the_verdict_is_not_assessed(self):
        # The state that must never be merged with "the library can separate them": one tests nothing.
        verdict = assess_candidate_set(_tied_pair(), library_available=False)

        self.assertEqual(STATE_NOT_ASSESSED, verdict["state"])
        self.assertEqual("library_unavailable", verdict["not_assessed_reason"])
        self.assertEqual("unknown", verdict["ceiling_effect"])

    def test_no_covering_class_means_this_run_was_the_limitation(self):
        verdict = assess_candidate_set(_tied_pair(), library_available=True, covering_class=None)

        self.assertEqual(STATE_QUERY_NON_DISCRIMINATION, verdict["state"])
        self.assertEqual("run_limited", verdict["ceiling_effect"])

    def test_a_covering_class_without_discriminating_ions_lowers_the_ceiling(self):
        verdict = assess_candidate_set(
            _tied_pair(),
            library_available=True,
            covering_class={
                "ambiguity_class_id": "urn:class:1",
                "members": [{"reference_spectrum_id": "a"}, {"reference_spectrum_id": "b"}],
                "discriminating_mz": [],
            },
        )

        self.assertEqual(STATE_LIBRARY_INDISTINGUISHABLE, verdict["state"])
        self.assertEqual("lowers_to_substructure_or_class", verdict["ceiling_effect"])
        self.assertEqual("urn:class:1", verdict["ambiguity_class_id"])
        self.assertIn("no discriminating product ion", verdict["state_reason"])

    def test_a_covering_class_with_discriminating_ions_says_so(self):
        verdict = assess_candidate_set(
            _tied_pair(),
            library_available=True,
            covering_class={
                "ambiguity_class_id": "urn:class:2",
                "members": [{"reference_spectrum_id": "a"}, {"reference_spectrum_id": "b"}],
                "discriminating_mz": [{"mz": 151.0}],
            },
        )

        self.assertEqual(STATE_LIBRARY_INDISTINGUISHABLE, verdict["state"])
        self.assertIn("threshold effect", verdict["state_reason"])


class DefinitionTests(unittest.TestCase):
    def test_the_rules_digest_changes_with_any_threshold(self):
        base = AssessmentDefinition()
        self.assertNotEqual(base.rules_sha256, AssessmentDefinition(mass_separation_mda=5.0).rules_sha256)
        self.assertNotEqual(
            base.rules_sha256, AssessmentDefinition(discrimination_margin=0.2).rules_sha256
        )
        # A different label over the same rules is the same rules.
        self.assertEqual(
            base.rules_sha256, AssessmentDefinition(label="something-else").rules_sha256
        )

    def test_as_rules_round_trips(self):
        definition = AssessmentDefinition(mass_separation_mda=2.5, discrimination_margin=0.05)
        self.assertEqual(AssessmentDefinition(**definition.as_rules()), definition)


class RunLevelTests(unittest.TestCase):
    def test_an_unknown_run_is_an_error_not_an_empty_success(self):
        with tempfile.TemporaryDirectory() as directory:
            report = assess_run_candidates(Path(directory) / "catalog.sqlite", "urn:msdial:run:none")

            self.assertFalse(report.valid)
            self.assertEqual(0, report.subjects)


if __name__ == "__main__":
    unittest.main()
