"""Tests for the spectral similarity kernel, with every convention pinned to a number."""

from __future__ import annotations

import math
import random
import unittest

from msdial_spectrum_catalog.similarity import (
    Peaks,
    compare,
    entropy_similarity,
    informative_peak_count,
    msdial_weighted_dot_product,
    normalize_to_base_peak,
    simple_cosine,
    spectral_entropy,
    weighted_cosine,
)


TOLERANCE = 0.05

# One double-precision step at 1.0 is 2**-52. The symmetric measures accumulate identical terms in
# identical order, so any residue below a couple of steps is rounding, not an ordering difference.
COSINE_ULP = 4e-16
ENTROPY_ULP = 4e-15

THREE_PEAKS: Peaks = [(100.0, 100.0), (150.0, 50.0), (200.0, 10.0)]
FIVE_PEAKS: Peaks = [(100.0, 100.0), (150.0, 80.0), (200.0, 60.0), (250.0, 40.0), (300.0, 20.0)]
FIVE_PEAKS_DISJOINT: Peaks = [(120.0, 100.0), (170.0, 80.0), (220.0, 60.0), (270.0, 40.0), (320.0, 20.0)]

# (label, peaks_a, peaks_b) covering the shapes an ambiguity class actually meets.
SYMMETRY_CASES: list[tuple[str, Peaks, Peaks]] = [
    ("identical", THREE_PEAKS, list(THREE_PEAKS)),
    ("partial_overlap", THREE_PEAKS, [(100.0, 80.0), (175.0, 60.0), (200.0, 30.0)]),
    ("disjoint", THREE_PEAKS, [(120.0, 100.0), (170.0, 50.0), (220.0, 10.0)]),
    ("one_against_five", [(100.0, 100.0)], FIVE_PEAKS),
    ("wildly_different_scales", [(100.0, 1e-3), (150.0, 5e-4)], [(100.0, 1e9), (150.0, 5e8)]),
    ("straddling_the_tolerance", [(100.0, 100.0), (200.0, 50.0)], [(100.049, 100.0), (200.051, 50.0)]),
    ("duplicate_mz", [(100.0, 60.0), (100.0, 40.0), (150.0, 50.0)], [(100.0, 100.0), (150.0, 50.0)]),
    ("unsorted_input", [(200.0, 10.0), (100.0, 100.0), (150.0, 50.0)], [(150.0, 50.0), (100.0, 100.0)]),
    ("huge_dynamic_range", [(100.0, 1e6), (500.0, 1.0), (900.0, 3.0)], [(100.0, 1e6), (500.0, 2.0), (900.0, 1.0)]),
    ("single_peaks_inside_tolerance", [(100.0, 1.0)], [(100.04, 1.0)]),
    ("single_peaks_outside_tolerance", [(100.0, 1.0)], [(100.06, 1.0)]),
    ("very_large_mz", [(1.0e6, 1.0), (2.0e6, 0.5)], [(1.0e6, 1.0), (2.0e6, 0.25)]),
]

# The one-ulp residue below is not shape-dependent noise; it is the division order in
# `covariance * covariance / scalar_a / scalar_b`. Swapping the arguments swaps the two divisions, and
# (x / a) / b is not always the same double as (x / b) / a. Dividing by the product would be exact.
DIVISION_ORDER_CASES: list[tuple[str, Peaks, Peaks]] = [
    ("partial_overlap", THREE_PEAKS, [(100.0, 80.0), (175.0, 60.0), (200.0, 30.0)]),
    (
        "one_against_a_decaying_five",
        [(100.0, 100.0)],
        [(100.0, 100.0), (150.0, 50.0), (200.0, 10.0), (250.0, 5.0), (300.0, 1.0)],
    ),
    ("unsorted_input", [(200.0, 10.0), (100.0, 100.0), (150.0, 50.0)], [(150.0, 50.0), (100.0, 100.0)]),
]

SIMPLE_COSINE_DIVISION_ORDER: tuple[Peaks, Peaks] = (
    [(238.0, 164.0), (311.0, 661.0)],
    [(154.0, 320.0), (311.0, 785.0)],
)


def _random_spectra(seed: int, count: int) -> list[tuple[Peaks, Peaks]]:
    """Build a reproducible set of random spectrum pairs for a symmetry sweep."""
    generator = random.Random(seed)
    pairs = []
    for _ in range(count):
        spectra = []
        for _ in range(2):
            size = generator.randint(1, 8)
            spectra.append(
                sorted(
                    (round(generator.uniform(50.0, 400.0), 3), generator.uniform(1e-3, 1e3))
                    for _ in range(size)
                )
            )
        pairs.append((spectra[0], spectra[1]))
    return pairs


class SymmetryTests(unittest.TestCase):
    """The class relation must not depend on argument order."""

    def test_weighted_cosine_is_symmetric_for_every_shape(self):
        for label, peaks_a, peaks_b in SYMMETRY_CASES:
            with self.subTest(label):
                forward, matched_forward = weighted_cosine(peaks_a, peaks_b, TOLERANCE)
                reverse, matched_reverse = weighted_cosine(peaks_b, peaks_a, TOLERANCE)
                self.assertEqual(matched_forward, matched_reverse)
                self.assertIsNotNone(forward)
                self.assertAlmostEqual(forward, reverse, delta=COSINE_ULP)

    def test_simple_cosine_is_symmetric_for_every_shape(self):
        for label, peaks_a, peaks_b in SYMMETRY_CASES:
            with self.subTest(label):
                forward = simple_cosine(peaks_a, peaks_b, TOLERANCE)
                reverse = simple_cosine(peaks_b, peaks_a, TOLERANCE)
                self.assertIsNotNone(forward)
                self.assertAlmostEqual(forward, reverse, delta=COSINE_ULP)

    def test_entropy_similarity_is_symmetric_for_every_shape(self):
        for label, peaks_a, peaks_b in SYMMETRY_CASES:
            with self.subTest(label):
                forward = entropy_similarity(peaks_a, peaks_b)
                reverse = entropy_similarity(peaks_b, peaks_a)
                self.assertIsNotNone(forward)
                self.assertAlmostEqual(forward, reverse, delta=ENTROPY_ULP)

    def test_symmetry_survives_a_randomized_sweep(self):
        for peaks_a, peaks_b in _random_spectra(20260904, 1500):
            forward, matched_forward = weighted_cosine(peaks_a, peaks_b, TOLERANCE)
            reverse, matched_reverse = weighted_cosine(peaks_b, peaks_a, TOLERANCE)
            self.assertEqual(matched_forward, matched_reverse, (peaks_a, peaks_b))
            self.assertAlmostEqual(forward, reverse, delta=COSINE_ULP, msg=(peaks_a, peaks_b))
            self.assertAlmostEqual(
                simple_cosine(peaks_a, peaks_b, TOLERANCE),
                simple_cosine(peaks_b, peaks_a, TOLERANCE),
                delta=COSINE_ULP,
                msg=(peaks_a, peaks_b),
            )
            self.assertAlmostEqual(
                entropy_similarity(peaks_a, peaks_b),
                entropy_similarity(peaks_b, peaks_a),
                delta=ENTROPY_ULP,
                msg=(peaks_a, peaks_b),
            )

    def test_disjoint_spectra_are_symmetric_and_score_zero(self):
        forward, matched = weighted_cosine(FIVE_PEAKS, FIVE_PEAKS_DISJOINT, TOLERANCE)
        reverse, _ = weighted_cosine(FIVE_PEAKS_DISJOINT, FIVE_PEAKS, TOLERANCE)
        self.assertEqual(forward, 0.0)
        self.assertEqual(reverse, 0.0)
        self.assertEqual(matched, 0)
        self.assertEqual(entropy_similarity(FIVE_PEAKS, FIVE_PEAKS_DISJOINT), 0.0)
        self.assertEqual(entropy_similarity(FIVE_PEAKS_DISJOINT, FIVE_PEAKS), 0.0)


class BitIdentityTests(unittest.TestCase):
    """Pin the one place where symmetry is only rounding-close, so a fix is noticed."""

    def test_weighted_cosine_is_bit_identical_under_argument_swap(self):
        """The cases that used to differ by one double step now agree exactly.

        `weighted_cosine` divides by the product of the two scalars rather than dividing twice, because
        `(x / a) / b` and `(x / b) / a` are not always the same double. IEEE-754 multiplication is
        commutative, so the product makes the relation exactly symmetric. That matters at a threshold
        boundary, where a one-ulp difference could make `sim(A, B)` pass while `sim(B, A)` fails.
        """
        for label, peaks_a, peaks_b in DIVISION_ORDER_CASES:
            with self.subTest(label):
                self.assertEqual(
                    weighted_cosine(peaks_a, peaks_b, TOLERANCE)[0],
                    weighted_cosine(peaks_b, peaks_a, TOLERANCE)[0],
                )

    def test_simple_cosine_is_bit_identical_under_argument_swap(self):
        peaks_a, peaks_b = SIMPLE_COSINE_DIVISION_ORDER
        self.assertEqual(
            simple_cosine(peaks_a, peaks_b, TOLERANCE),
            simple_cosine(peaks_b, peaks_a, TOLERANCE),
        )

    def test_entropy_similarity_is_bit_identical_under_argument_swap(self):
        """Both order dependencies are gone.

        The two single-spectrum entropies are summed before being subtracted, and the binning helpers
        sort each frame before reducing it with `math.fsum`, so neither the subtraction order nor the
        concatenation order inside `_combined` can move the last bit.
        """
        peaks_a, peaks_b = THREE_PEAKS, [(100.0, 80.0), (175.0, 60.0), (200.0, 30.0)]
        self.assertEqual(entropy_similarity(peaks_a, peaks_b), entropy_similarity(peaks_b, peaks_a))

    def test_no_residue_remains_on_the_cases_that_had_one(self):
        for label, peaks_a, peaks_b in DIVISION_ORDER_CASES:
            with self.subTest(label):
                forward = weighted_cosine(peaks_a, peaks_b, TOLERANCE)[0]
                reverse = weighted_cosine(peaks_b, peaks_a, TOLERANCE)[0]
                self.assertEqual(forward, reverse)

class MsdialAsymmetryTests(unittest.TestCase):
    """The faithful port is asymmetric on purpose; these numbers say so out loud."""

    def test_the_faithful_port_is_asymmetric_by_design(self):
        """Pin both directions of GetWeightedDotProduct so a symmetry 'fix' fails loudly.

        MS-DIAL derives peakCountPenalty from the reference side only and applies its 0.01 intensity
        cutoff to the measured side only. Three measured peaks against one reference peak therefore
        score 0.3846 with a 0.75 penalty, while the same pair with the sides exchanged scores 0.88 with
        a 0.88 penalty. Anyone tempted to make these agree should read the module docstring first: the
        symmetric relation is `weighted_cosine`, and this function exists only so a class edge can be
        held against an MS-DIAL annotation score.
        """
        forward = msdial_weighted_dot_product(THREE_PEAKS, [(100.0, 100.0)], TOLERANCE)
        reverse = msdial_weighted_dot_product([(100.0, 100.0)], THREE_PEAKS, TOLERANCE)
        self.assertEqual(forward, (0.3846153846153847, 0.75))
        self.assertEqual(reverse, (0.88, 0.88))
        self.assertNotEqual(forward[0], reverse[0])
        self.assertNotEqual(forward[1], reverse[1])

    def test_the_penalty_ladder_matches_the_reference_peak_count(self):
        expected = {1: 0.75, 2: 0.88, 3: 0.94, 4: 0.97, 5: 1.0, 6: 1.0}
        for count, penalty in expected.items():
            with self.subTest(count):
                peaks = [(100.0 + 50.0 * index, 100.0 - 5.0 * index) for index in range(count)]
                value, applied = msdial_weighted_dot_product(peaks, peaks, TOLERANCE)
                self.assertEqual(applied, penalty)
                self.assertEqual(value, penalty)
                # The symmetric relation carries no penalty, so a self-comparison stays exactly 1.
                self.assertEqual(weighted_cosine(peaks, peaks, TOLERANCE)[0], 1.0)

    def test_the_squared_value_is_never_read_as_a_cosine(self):
        two_peaks: Peaks = [(100.0, 100.0), (150.0, 50.0)]
        self.assertEqual(weighted_cosine(two_peaks, two_peaks, TOLERANCE)[0], 1.0)
        self.assertEqual(msdial_weighted_dot_product(two_peaks, two_peaks, TOLERANCE), (0.88, 0.88))

    def test_the_squared_value_is_the_square_of_the_cosine_when_no_penalty_applies(self):
        peaks_b: Peaks = [(100.0, 90.0), (150.0, 20.0), (200.0, 70.0), (250.0, 55.0), (300.0, 30.0)]
        cosine, _ = weighted_cosine(FIVE_PEAKS, peaks_b, TOLERANCE)
        squared, penalty = msdial_weighted_dot_product(FIVE_PEAKS, peaks_b, TOLERANCE)
        self.assertEqual(penalty, 1.0)
        self.assertEqual(cosine, 0.9625646738468938)
        self.assertEqual(squared, 0.926530751337977)
        self.assertEqual(cosine * cosine, squared)


class SelfComparisonTests(unittest.TestCase):
    """A spectrum against itself is exactly 1.0, and never more."""

    SHAPES: list[tuple[str, Peaks]] = [
        ("single_peak", [(180.0876, 4210.0)]),
        ("three_peaks", THREE_PEAKS),
        ("five_peaks", FIVE_PEAKS),
        ("huge_dynamic_range", [(100.0, 1.0), (301.1408, 5.0e8), (700.5, 3.0)]),
        ("tiny_intensities", [(100.0, 1e-9), (200.0, 4e-9)]),
        ("duplicate_mz", [(100.0, 60.0), (100.0, 40.0), (150.0, 50.0)]),
        ("very_large_mz", [(1.0e6, 1.0), (2.0e6, 0.5)]),
    ]

    def test_every_measure_returns_exactly_one_for_a_spectrum_against_itself(self):
        for label, peaks in self.SHAPES:
            with self.subTest(label):
                self.assertEqual(weighted_cosine(peaks, list(peaks), TOLERANCE)[0], 1.0)
                self.assertEqual(simple_cosine(peaks, list(peaks), TOLERANCE), 1.0)
                self.assertEqual(entropy_similarity(peaks, list(peaks)), 1.0)

    def test_the_clamp_holds_on_a_case_that_overshoots_one(self):
        """Pin a self-comparison whose unclamped square is 1.0000000000000002.

        `covariance` accumulates sqrt(i * i) * mz while the scalars accumulate i * mz, and sqrt(i * i)
        is not always the same double as i, so a perfect match can exceed 1 before the clamp.
        """
        peaks: Peaks = [(189.0, 274.0), (820.0, 204.0)]
        normalized = normalize_to_base_peak(peaks)
        covariance = scalar_a = scalar_b = 0.0
        for mz, intensity in normalized:
            covariance += math.sqrt(intensity * intensity) * mz
            scalar_a += intensity * mz
            scalar_b += intensity * mz
        self.assertGreater(covariance * covariance / scalar_a / scalar_b, 1.0)
        self.assertEqual(weighted_cosine(peaks, list(peaks), TOLERANCE)[0], 1.0)

    def test_no_random_pair_ever_exceeds_one(self):
        for peaks_a, peaks_b in _random_spectra(11, 1500):
            cosine, _ = weighted_cosine(peaks_a, peaks_b, TOLERANCE)
            self.assertLessEqual(cosine, 1.0, (peaks_a, peaks_b))
            self.assertGreaterEqual(cosine, 0.0, (peaks_a, peaks_b))
            self.assertLessEqual(simple_cosine(peaks_a, peaks_b, TOLERANCE), 1.0)


class AdmissibilityGateTests(unittest.TestCase):
    """A gate failure is 'not comparable', which is not the same statement as 'not similar'."""

    GATE = {
        "tolerance": TOLERANCE,
        "minimum_informative_peaks": 3,
        "minimum_matched_peaks": 2,
        "relative_floor": 0.05,
    }

    THIN: Peaks = [(100.0, 100.0), (150.0, 90.0)]
    RICH: Peaks = [(100.0, 100.0), (150.0, 90.0), (200.0, 80.0), (250.0, 70.0), (300.0, 60.0)]
    DISJOINT: Peaks = [(120.0, 100.0), (170.0, 90.0), (220.0, 80.0), (270.0, 70.0), (320.0, 60.0)]
    ONE_MATCH: Peaks = [(100.0, 100.0), (400.0, 90.0), (450.0, 80.0), (500.0, 70.0)]

    def test_too_few_informative_peaks_is_not_comparable(self):
        result = compare(self.THIN, self.RICH, **self.GATE)
        self.assertFalse(result.comparable)
        self.assertIn("insufficient_informative_peaks", result.reason)
        self.assertIsNone(result.weighted_cosine)
        self.assertIsNone(result.entropy_similarity)
        self.assertIsNone(result.msdial_weighted_squared)
        self.assertEqual(result.informative_peak_count_a, 2)
        self.assertEqual(result.informative_peak_count_b, 5)

    def test_too_few_matched_peaks_is_not_comparable(self):
        result = compare(self.RICH, self.DISJOINT, **self.GATE)
        self.assertFalse(result.comparable)
        self.assertIn("insufficient_matched_peaks", result.reason)
        self.assertIsNone(result.weighted_cosine)
        self.assertIsNone(result.entropy_similarity)
        self.assertIsNone(result.msdial_weighted_squared)
        self.assertEqual(result.matched_peak_count, 0)
        self.assertEqual(result.informative_peak_count_a, 5)
        self.assertEqual(result.informative_peak_count_b, 5)

    def test_one_matched_peak_below_the_minimum_is_not_comparable(self):
        result = compare(self.RICH, self.ONE_MATCH, **self.GATE)
        self.assertFalse(result.comparable)
        self.assertIn("insufficient_matched_peaks", result.reason)
        self.assertEqual(result.matched_peak_count, 1)
        self.assertEqual(result.informative_peak_count_a, 5)
        self.assertEqual(result.informative_peak_count_b, 4)

    def test_the_two_failure_reasons_are_distinguishable(self):
        informative = compare(self.THIN, self.RICH, **self.GATE).reason
        matched = compare(self.RICH, self.DISJOINT, **self.GATE).reason
        self.assertNotEqual(informative, matched)
        self.assertNotIn("insufficient_matched_peaks", informative)
        self.assertNotIn("insufficient_informative_peaks", matched)
        self.assertIn("minimum 3", informative)
        self.assertIn("minimum 2", matched)

    def test_both_failure_modes_report_informative_counts(self):
        for label, result in (
            ("informative", compare(self.THIN, self.RICH, **self.GATE)),
            ("matched", compare(self.RICH, self.DISJOINT, **self.GATE)),
        ):
            with self.subTest(label):
                self.assertGreater(result.informative_peak_count_a, 0)
                self.assertGreater(result.informative_peak_count_b, 0)

    def test_a_gate_failure_is_not_reported_as_a_low_similarity(self):
        """A pair the gate rejects must carry no number at all, not a small one."""
        direct, matched = weighted_cosine(self.RICH, self.DISJOINT, TOLERANCE)
        self.assertEqual(direct, 0.0)
        self.assertEqual(matched, 0)
        result = compare(self.RICH, self.DISJOINT, **self.GATE)
        self.assertFalse(result.comparable)
        self.assertIsNone(result.weighted_cosine)

    def test_an_admissible_pair_carries_every_number(self):
        result = compare(self.RICH, list(self.RICH), **self.GATE)
        self.assertTrue(result.comparable)
        self.assertIsNone(result.reason)
        self.assertEqual(result.weighted_cosine, 1.0)
        self.assertEqual(result.entropy_similarity, 1.0)
        self.assertEqual(result.matched_peak_count, 5)
        self.assertEqual(result.msdial_weighted_squared, 1.0)
        self.assertEqual(result.msdial_peak_count_penalty, 1.0)

    def test_the_gate_is_evaluated_inside_the_mass_range(self):
        result = compare(self.RICH, list(self.RICH), mass_begin=90.0, mass_end=210.0, **self.GATE)
        self.assertTrue(result.comparable)
        self.assertEqual(result.informative_peak_count_a, 3)
        self.assertEqual(result.matched_peak_count, 3)

    def test_the_relative_floor_selects_which_peaks_count(self):
        peaks: Peaks = [(100.0, 100.0), (150.0, 3.0), (200.0, 2.0)]
        self.assertEqual(informative_peak_count(peaks, 0.05), 1)
        self.assertEqual(informative_peak_count(peaks, 0.01), 3)
        self.assertEqual(informative_peak_count([(100.0, 100.0), (150.0, 5.0)], 0.05), 2)


class EntropyParityTests(unittest.TestCase):
    """Shannon entropy in bits, and the Li-Fiehn combination MS-DIAL uses."""

    def test_entropy_of_hand_computable_spectra(self):
        self.assertEqual(spectral_entropy([(100.0, 5.0)]), 0.0)
        self.assertEqual(spectral_entropy([(100.0, 1.0), (200.0, 1.0)]), 1.0)
        self.assertEqual(
            spectral_entropy([(100.0, 7.0), (200.0, 7.0), (300.0, 7.0), (400.0, 7.0)]), 2.0
        )
        self.assertEqual(spectral_entropy([(100.0, 1.0), (200.0, 3.0)]), 0.8112781244591328)

    def test_entropy_of_a_degenerate_spectrum_is_zero(self):
        self.assertEqual(spectral_entropy([]), 0.0)
        self.assertEqual(spectral_entropy([(100.0, 0.0), (200.0, 0.0)]), 0.0)

    def test_entropy_similarity_on_a_hand_computed_case(self):
        peaks_a: Peaks = [(100.0, 1.0), (200.0, 1.0)]
        peaks_b: Peaks = [(100.0, 1.0), (200.0, 3.0)]
        # Unit-normalized, the combined spectrum is (0.375, 0.625) after the 0.5 factor.
        merged = spectral_entropy([(100.0, 0.375), (200.0, 0.625)])
        # 0.9512050593046014 is the correctly rounded double: the exact value evaluated to 60 decimal
        # digits is 0.95120505930460146741..., which rounds up. The naive float expression below lands
        # one ulp low, which is why the assertion on it is approximate and the assertion on the function
        # is exact -- entropy_similarity reduces each frame with math.fsum.
        expected = 1.0 - (2.0 * merged - 1.0 - 0.8112781244591328) * 0.5
        self.assertAlmostEqual(expected, 0.9512050593046014, places=15)
        self.assertEqual(entropy_similarity(peaks_a, peaks_b), 0.9512050593046014)

    def test_entropy_similarity_of_disjoint_spectra_is_zero(self):
        self.assertEqual(entropy_similarity([(100.0, 1.0)], [(300.0, 1.0)]), 0.0)

    def test_entropy_similarity_ignores_a_shared_scale(self):
        peaks: Peaks = [(100.0, 3.0), (200.0, 1.0)]
        scaled: Peaks = [(100.0, 3.0e6), (200.0, 1.0e6)]
        self.assertEqual(entropy_similarity(peaks, scaled), 1.0)


class FixedGridBinningTests(unittest.TestCase):
    """MS-DIAL bins entropy on a fixed int(mass / bin) grid; that is reproduced, not repaired."""

    def test_peaks_two_millidalton_apart_can_land_in_different_frames(self):
        """Pin the fixed-grid artefact: 0.002 apart scores 0, 0.048 apart scores 1.

        `int(99.999 / 0.05)` is 1999 and `int(100.001 / 0.05)` is 2000, so two peaks 0.002 apart that
        straddle the boundary are treated as sharing nothing, while 100.001 and 100.049 -- 24 times
        further apart -- share a frame and are treated as identical. This is
        SpectrumHandler.GetBinnedSpectrum's behaviour reproduced deliberately, so entropy values stay
        comparable with MS-DIAL's own. It is not a defect in this port, and it is the reason
        entropy_similarity is reported alongside the tolerance-window cosine rather than instead of it.
        """
        self.assertEqual(int(99.999 / 0.05), 1999)
        self.assertEqual(int(100.001 / 0.05), 2000)
        self.assertEqual(int(100.049 / 0.05), 2000)
        self.assertEqual(entropy_similarity([(100.001, 1.0)], [(99.999, 1.0)]), 0.0)
        self.assertEqual(entropy_similarity([(100.001, 1.0)], [(100.049, 1.0)]), 1.0)

    def test_the_tolerance_window_cosine_does_not_share_the_artefact(self):
        straddling = weighted_cosine([(100.001, 1.0)], [(99.999, 1.0)], TOLERANCE)[0]
        self.assertEqual(straddling, 1.0)

    def test_which_pairs_share_a_frame_depends_on_the_grid_not_the_gap(self):
        # bin_width 0.3 puts both peaks in frame 333; bin_width 2.0 splits them across 49 and 50. The
        # gap never changed, only where the fixed grid happens to fall.
        self.assertEqual(entropy_similarity([(100.001, 1.0)], [(99.999, 1.0)], bin_width=0.3), 1.0)
        self.assertEqual(entropy_similarity([(100.001, 1.0)], [(99.999, 1.0)], bin_width=2.0), 0.0)


class EdgeCaseTests(unittest.TestCase):
    """Establish what each degenerate input returns, so callers can rely on it."""

    PEAKS: Peaks = [(100.0, 100.0), (150.0, 50.0)]

    def test_missing_or_unusable_peaks_give_none_not_zero(self):
        cases: list[tuple[str, Peaks, Peaks]] = [
            ("both_empty", [], []),
            ("left_empty", [], self.PEAKS),
            ("right_empty", self.PEAKS, []),
            ("all_zero_intensities", [(100.0, 0.0), (150.0, 0.0)], self.PEAKS),
            ("negative_intensities", [(100.0, -5.0), (150.0, -1.0)], self.PEAKS),
            ("zero_on_both_sides", [(100.0, 0.0)], [(100.0, 0.0)]),
        ]
        for label, peaks_a, peaks_b in cases:
            with self.subTest(label):
                self.assertEqual(weighted_cosine(peaks_a, peaks_b, TOLERANCE), (None, 0))
                self.assertIsNone(simple_cosine(peaks_a, peaks_b, TOLERANCE))

    def test_nothing_left_inside_the_mass_range_gives_none(self):
        self.assertEqual(
            weighted_cosine(self.PEAKS, self.PEAKS, TOLERANCE, mass_begin=500.0, mass_end=600.0),
            (None, 0),
        )
        self.assertIsNone(
            simple_cosine(self.PEAKS, self.PEAKS, TOLERANCE, mass_begin=500.0, mass_end=600.0)
        )
        self.assertEqual(
            msdial_weighted_dot_product(self.PEAKS, self.PEAKS, TOLERANCE, mass_begin=500.0),
            (None, 1.0),
        )

    def test_entropy_similarity_gives_none_for_unusable_input(self):
        self.assertIsNone(entropy_similarity([], self.PEAKS))
        self.assertIsNone(entropy_similarity(self.PEAKS, []))
        self.assertIsNone(entropy_similarity([], []))
        self.assertIsNone(entropy_similarity([(100.0, 0.0)], self.PEAKS))
        self.assertIsNone(entropy_similarity([(100.0, -5.0), (150.0, 1.0)], self.PEAKS))

    def test_a_single_peak_each_is_handled(self):
        self.assertEqual(weighted_cosine([(100.0, 1.0)], [(100.0, 7.0)], TOLERANCE), (1.0, 1))
        self.assertEqual(simple_cosine([(100.0, 1.0)], [(100.0, 7.0)], TOLERANCE), 1.0)
        self.assertEqual(entropy_similarity([(100.0, 1.0)], [(100.0, 7.0)]), 1.0)
        self.assertEqual(
            msdial_weighted_dot_product([(100.0, 1.0)], [(100.0, 7.0)], TOLERANCE), (0.75, 0.75)
        )

    def test_unsorted_input_is_sorted_before_alignment(self):
        unsorted_peaks: Peaks = [(200.0, 10.0), (100.0, 100.0), (150.0, 50.0)]
        self.assertEqual(
            weighted_cosine(unsorted_peaks, THREE_PEAKS, TOLERANCE),
            weighted_cosine(THREE_PEAKS, THREE_PEAKS, TOLERANCE),
        )
        self.assertEqual(simple_cosine(unsorted_peaks, THREE_PEAKS, TOLERANCE), 1.0)

    def test_duplicate_mz_values_are_summed_into_one_frame(self):
        split: Peaks = [(100.0, 60.0), (100.0, 40.0), (150.0, 50.0)]
        merged: Peaks = [(100.0, 100.0), (150.0, 50.0)]
        self.assertEqual(weighted_cosine(split, merged, TOLERANCE), (1.0, 2))
        self.assertEqual(simple_cosine(split, merged, TOLERANCE), 1.0)

    def test_very_large_mz_values_do_not_break_the_weighting(self):
        large: Peaks = [(1.0e6, 1.0), (2.0e6, 0.5)]
        self.assertEqual(weighted_cosine(large, list(large), TOLERANCE), (1.0, 2))
        self.assertLess(weighted_cosine(large, [(1.0e6, 1.0), (2.0e6, 5.0)], TOLERANCE)[0], 1.0)

    def test_non_positive_peaks_are_dropped_by_normalization(self):
        self.assertEqual(normalize_to_base_peak([]), [])
        self.assertEqual(normalize_to_base_peak([(100.0, 0.0), (150.0, -3.0)]), [])
        self.assertEqual(
            normalize_to_base_peak([(100.0, 50.0), (150.0, 25.0), (200.0, 0.0), (250.0, -3.0)]),
            [(100.0, 1.0), (150.0, 0.5)],
        )
        self.assertEqual(informative_peak_count([], 0.05), 0)

    def test_the_intensity_cutoff_hides_minor_peaks_under_a_huge_base_peak(self):
        """Pin that the inherited 0.01 cutoff can score two different spectra as identical.

        With a base peak six orders of magnitude above the rest, every minor frame falls below the 0.01
        relative cutoff and is dropped from both sides, so spectra that differ only in those peaks score
        exactly 1.0. That is MS-DIAL's cutoff, and it is why the admissibility gate counts informative
        peaks with its own relative floor rather than trusting the cosine alone.
        """
        peaks_a: Peaks = [(100.0, 1.0e6), (500.0, 1.0), (900.0, 3.0)]
        peaks_b: Peaks = [(100.0, 1.0e6), (500.0, 2.0), (900.0, 1.0)]
        self.assertEqual(weighted_cosine(peaks_a, peaks_b, TOLERANCE), (1.0, 1))
        self.assertEqual(informative_peak_count(peaks_a, 0.01), 1)

    def test_spectral_entropy_ignores_negative_intensities(self):
        """Non-positive peaks are dropped before the total is taken, so entropy is never negative.

        Dividing by a signed total would let a share exceed 1 and return a negative "entropy", which is
        a wrong number rather than a degenerate one. These two cases returned -2.0 and about -2.3e8
        before the guard was added.
        """
        # One negative and one positive peak leaves a single positive peak: a degenerate distribution,
        # entropy 0, not a negative number.
        self.assertEqual(spectral_entropy([(100.0, -5.0), (150.0, 10.0)]), 0.0)
        self.assertEqual(spectral_entropy([(100.0, -9.999999), (150.0, 10.0)]), 0.0)
        # An all-negative spectrum has nothing to measure.
        self.assertEqual(spectral_entropy([(100.0, -1.0), (150.0, -2.0)]), 0.0)
        # The guard must not disturb a normal spectrum.
        self.assertEqual(spectral_entropy([(100.0, 1.0), (150.0, 1.0)]), 1.0)

    def test_comparable_always_means_a_number_is_present(self):
        """Zero minimums no longer let an empty pair through as comparable.

        `compare` gates on the minimums it is given, so minimums of 0 once admitted two empty spectra
        and reported comparable=True with every similarity None. A caller reading `comparable` as "there
        is a number here" would have been wrong, so `compare` now refuses when no cosine could be
        computed at all.
        """
        result = compare(
            [],
            [],
            tolerance=TOLERANCE,
            minimum_informative_peaks=0,
            minimum_matched_peaks=0,
            relative_floor=0.05,
        )
        self.assertFalse(result.comparable)
        self.assertEqual(result.reason, "no_usable_peaks_in_range")
        self.assertIsNone(result.weighted_cosine)
        # And whenever comparable is True, every headline number is populated.
        good = compare(
            THREE_PEAKS,
            THREE_PEAKS,
            tolerance=TOLERANCE,
            minimum_informative_peaks=1,
            minimum_matched_peaks=1,
            relative_floor=0.01,
        )
        self.assertTrue(good.comparable)
        self.assertIsNotNone(good.weighted_cosine)
        self.assertIsNotNone(good.entropy_similarity)

    def test_msdial_port_returns_zero_when_no_measured_peak_survives_the_cutoff(self):
        value, penalty = msdial_weighted_dot_product(
            [(100.0, 1.0e6), (500.0, 1.0)], [(700.0, 5.0)], TOLERANCE
        )
        self.assertEqual(value, 0.0)
        self.assertEqual(penalty, 0.75)


class MonotonicityTests(unittest.TestCase):
    """Moving one spectrum away from another must never raise the similarity."""

    def test_mixing_in_a_disjoint_spectrum_is_non_increasing(self):
        sequence = []
        for step in range(11):
            fraction = step / 10.0
            mixed = sorted(
                [(mz, intensity * (1.0 - fraction)) for mz, intensity in FIVE_PEAKS]
                + [(mz, intensity * fraction) for mz, intensity in FIVE_PEAKS_DISJOINT]
            )
            cosine, _ = weighted_cosine(FIVE_PEAKS, mixed, TOLERANCE)
            sequence.append(cosine)
        self.assertEqual(sequence[0], 1.0)
        self.assertEqual(sequence[-1], 0.0)
        for earlier, later in zip(sequence, sequence[1:]):
            self.assertGreaterEqual(earlier, later, sequence)

    def test_walking_one_peak_out_of_the_tolerance_window_is_non_increasing(self):
        sequence = []
        for shift in (0.0, 0.01, 0.02, 0.03, 0.04, 0.049, 0.06, 0.5, 5.0):
            shifted = FIVE_PEAKS[:-1] + [(FIVE_PEAKS[-1][0] + shift, FIVE_PEAKS[-1][1])]
            cosine, _ = weighted_cosine(FIVE_PEAKS, shifted, TOLERANCE)
            sequence.append(cosine)
        self.assertEqual(sequence[0], 1.0)
        self.assertLess(sequence[-1], 1.0)
        for earlier, later in zip(sequence, sequence[1:]):
            self.assertGreaterEqual(earlier, later, sequence)

    def test_narrowing_the_mass_range_cannot_invent_a_match(self):
        wide, wide_matched = weighted_cosine(FIVE_PEAKS, FIVE_PEAKS_DISJOINT, TOLERANCE)
        narrow, narrow_matched = weighted_cosine(
            FIVE_PEAKS, FIVE_PEAKS_DISJOINT, TOLERANCE, mass_begin=90.0, mass_end=210.0
        )
        self.assertEqual(wide, 0.0)
        self.assertEqual(narrow, 0.0)
        self.assertEqual(wide_matched, 0)
        self.assertEqual(narrow_matched, 0)


if __name__ == "__main__":
    unittest.main()
