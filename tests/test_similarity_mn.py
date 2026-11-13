"""Hand anchors, deterministic properties and identity grammar contracts."""
import ast
from dataclasses import asdict, replace
import itertools
import math
from pathlib import Path
import random
import unittest

from msdial_spectrum_catalog import similarity_mn as mn
from msdial_spectrum_catalog import similarity_mn_tools as mt


class KernelAnchors(unittest.TestCase):
    def setUp(self):
        self.a = mn.Spectrum(250, ((100, 3), (150, 4)))

    def test_single_shared_fragment(self):
        b = mn.Spectrum(264, ((100, 3), (180, 12)))
        mod = mn.msdial_moddot_legacy(self.a, b)
        self.assertEqual((mod.score, mod.matched_peak_count), (1, 1))
        self.assertEqual(mn.msdial_bonanza_v1(self.a, b).score, 9 / 169)
        self.assertAlmostEqual(mn.modified_cosine(self.a, b, tol=.025, assign="greedy").score,
                               0.1455213750217998, places=15)

    def test_shift_identity_and_entropy_has_no_shift(self):
        b = mn.Spectrum(264, ((100, 3), (164, 4)))
        results = [mn.msdial_moddot_legacy(self.a, b), mn.msdial_bonanza_v1(self.a, b)]
        results += [mn.modified_cosine(self.a, b, tol=.025, assign=x) for x in ("greedy", "optimal")]
        for result in results:
            self.assertEqual((result.score, result.matched_peak_count), (1, 2))
        self.assertLess(mn.entropy_li2021(self.a, b).score, 1)
        self.assertEqual(mn.msdial_all_v1(self.a, b), (1, 2, 1, 1))

    def test_equal_precursor_role_asymmetry(self):
        a = mn.Spectrum(200, ((100, 10), (100.02, 25)))
        b = mn.Spectrum(200, ((100.01, 20),))
        self.assertEqual(mn.msdial_match(a, b).e, "b")
        self.assertEqual(mn.msdial_bonanza_v1(a, b).score, 200 / 825)
        self.assertEqual(mn.msdial_bonanza_v1(b, a).score, 500 / 600)

    def test_strict_and_inclusive_tolerance(self):
        a, b = mn.Spectrum(200, ((100, 10),)), mn.Spectrum(200, ((100.25, 10),))
        self.assertEqual(mn.msdial_match(a, b, .25).pairs, ())
        for assign in ("greedy", "optimal"):
            result = mn.modified_cosine(a, b, tol=.25, assign=assign)
            self.assertEqual((result.score, result.matched_peak_count), (1, 1))

    def test_modified_cosine_boundary_asymmetry_matches_matchms(self):
        # Inclusive floating-point endpoints need not survive argument reversal.
        # The same candidate-window arithmetic is used by matchms 0.34.0.
        a = mn.Spectrum(200, ((0.01, 1),))
        b = mn.Spectrum(200, ((0.035, 1),))
        self.assertEqual(0.01 + 0.025, 0.035)
        self.assertGreater(0.035 - 0.025, 0.01)
        for assign in ("greedy", "optimal"):
            with self.subTest(assign=assign):
                forward = mn.modified_cosine(a, b, tol=.025, assign=assign)
                reverse = mn.modified_cosine(b, a, tol=.025, assign=assign)
                self.assertEqual((forward.score, forward.matched_peak_count), (1, 1))
                self.assertEqual((reverse.score, reverse.matched_peak_count), (0, 0))

    def test_small_shift_does_not_add_shifted_candidates(self):
        a, b = mn.Spectrum(200, ((100, 10),)), mn.Spectrum(200.02, ((100.04, 10),))
        self.assertEqual(mn.modified_cosine(a, b, tol=.025, assign="greedy").matched_peak_count, 0)

    def test_intensity_selection_and_first_tie(self):
        a, b = mn.Spectrum(200, ((90, 9), (100, 3))), mn.Spectrum(210, ((100, 10),))
        self.assertEqual(mn.msdial_match(a, b).pairs, ((0, 0, "neutral_loss"),))
        a = mn.Spectrum(200, ((99.99, 8), (100.01, 12)))
        self.assertEqual(mn.msdial_match(a, b).pairs, ((0, 0, "product"),))

    def test_low_level_preserves_input_order(self):
        a = mn.Spectrum(200, ((100.01, 12), (99.99, 8)))
        b = mn.Spectrum(210, ((100, 10),))
        self.assertEqual(mn.msdial_match(a, b).pairs, ((0, 0, "product"),))
        result = mn.score("msdial.bonanza.v1?prep=none&tol=0.025", a, b)
        self.assertEqual(result.score, 80 / (80 + 144))

    def test_nan_and_all_mode_zero_guard(self):
        a, b = mn.Spectrum(200, ((100, 0),)), mn.Spectrum(214, ((100, 0),))
        self.assertTrue(math.isnan(mn.msdial_moddot_legacy(a, b).score))
        self.assertTrue(math.isnan(mn.msdial_bonanza_v1(a, b).score))
        all_scores = mn.msdial_all_v1(a, b)
        self.assertTrue(math.isnan(all_scores[0]))
        self.assertEqual(all_scores[1:], (1, 0, 0))
        self.assertEqual(mn.modified_cosine(a, b, tol=.025, assign="greedy").score, 0)

    def test_empty_and_entropy(self):
        empty = mn.Spectrum(200, ())
        self.assertEqual(mn.msdial_all_v1(empty, self.a), (0, 0, 0, 0))
        self.assertEqual(mn.msdial_moddot_legacy(empty, self.a).score, 0)
        self.assertIsNone(mn.entropy_li2021(empty, self.a).score)
        a, b = mn.Spectrum(200, ((90, 1),)), mn.Spectrum(200, ((110, 1),))
        self.assertEqual(mn.entropy_li2021(a, b).score, 0)
        a, b = mn.Spectrum(200, ((90, 1), (110, 1))), mn.Spectrum(200, ((90, 1), (130, 1)))
        self.assertEqual(mn.entropy_li2021(a, b).score, .5)

    def test_greedy_can_lose_to_optimal(self):
        a = mn.Spectrum(200, ((100, 9), (100.04, 10)))
        b = mn.Spectrum(200, ((100.02, 10), (100.06, 9)))
        greedy = mn.modified_cosine(a, b, tol=.025, assign="greedy")
        optimal = mn.modified_cosine(a, b, tol=.025, assign="optimal")
        self.assertAlmostEqual(greedy.score, 100 / 181)
        self.assertAlmostEqual(optimal.score, 180 / 181)
        self.assertGreater(optimal.score - greedy.score, .05)

    def test_entropy_zero_and_nonfinite_intensities_are_not_comparable(self):
        invalid = [mn.Spectrum(200, ((100, 0),)),
                   mn.preprocess(mn.Spectrum(200, ((100, 0),)))]
        invalid += [mn.Spectrum(200, ((100, 1), (150, value)))
                    for value in (math.nan, math.inf, -math.inf)]
        for spectrum in invalid:
            for a, b in ((spectrum, spectrum), (spectrum, self.a), (self.a, spectrum)):
                with self.subTest(peaks=a.peaks, other=b.peaks):
                    self.assertIsNone(mn.entropy_li2021(a, b).score)

    def test_greedy_ties_reverse_stable_candidate_order(self):
        a = mn.Spectrum(200, ((100, 1), (100.01, 1)))
        b = mn.Spectrum(200, ((100, 1), (100.01, 1)))
        self.assertEqual(mn.modified_cosine(a, b, tol=.025, assign="greedy").details["pairs"], [[1, 1], [0, 0]])

    def test_nonzero_mz_and_intensity_powers_full_norm(self):
        a = mn.Spectrum(250, ((100, 9), (150, 16)))
        b = mn.Spectrum(264, ((100, 9), (180, 144)))
        result = mn.modified_cosine(a, b, tol=.025, assign="greedy", mz_power=1, intensity_power=.5)
        self.assertAlmostEqual(result.score, 90000 / (math.sqrt(450000) * math.sqrt(4755600)))

    def test_edge_gate_boundaries_and_nan(self):
        good = mn.KernelResult("", .5, "bonanza_ratio", 2, {})
        self.assertTrue(mn.msdial_edge_passes(good, -400))
        self.assertFalse(mn.msdial_edge_passes(good, 400.00001))
        self.assertFalse(mn.msdial_edge_passes(replace(good, matched_peak_count=1), 0))
        self.assertFalse(mn.msdial_edge_passes(replace(good, score=.49999), 0))
        self.assertFalse(mn.msdial_edge_passes(replace(good, score=None), 0))
        self.assertTrue(mn.msdial_edge_passes(replace(good, score=math.nan), 0))


class Preprocessing(unittest.TestCase):
    def test_base_is_before_filter_and_relative_equality(self):
        s = mn.Spectrum(200, ((100, 5), (150, 10), (250, 1000)))
        self.assertEqual(mn.preprocess(s).peaks, ((150, 100),))

    def test_first_apex_and_stable_sort_and_sum(self):
        s = mn.Spectrum(200, ((101.1, 5), (100.9, 5), (100.1, 5)))
        self.assertEqual(mn.preprocess(s).peaks, ((100.1, 100), (101.1, math.sqrt(5) / math.sqrt(10) * 100)))

    def test_refine_boundaries_and_no_precursor_removal(self):
        s = mn.Spectrum(200, ((-1, 2), (100, 4), (200, 4), (200.05, 4), (200.05001, 4), (10001, 4)))
        p = replace(mn.MSDIAL_MN_V1, frame=0, pow=1, scale=4)
        self.assertEqual(mn.preprocess(s, p).peaks, ((100, 4), (200, 4), (200.05, 4)))
        self.assertEqual(mn.preprocess(mn.Spectrum(200, ())).peaks, ())

    def test_python_switches_topk_ties_and_window_equality(self):
        p = replace(mn.MSDIAL_MN_V1, frame=0, topk=1, pow=1, scale=10)
        s = mn.Spectrum(200, ((50, 10), (75, 10), (100, 5), (100.00001, 4), (175, 3)))
        self.assertEqual(mn.preprocess(s, p).peaks, ((50, 10), (75, 10), (175, 3)))
        p = replace(p, precwin=17, topk=0)
        s = mn.Spectrum(200, ((182.99999, 10), (183, 10), (200, 10)))
        self.assertEqual(mn.preprocess(s, p).peaks, ((182.99999, 10),))

    def test_zero_normalization_is_ieee_nan(self):
        result = mn.preprocess(mn.Spectrum(200, ((100, 0),)))
        self.assertTrue(math.isnan(result.peaks[0][1]))

    def test_nonfinite_input_intensities_are_rejected_before_filtering(self):
        for intensity in (math.nan, math.inf, -math.inf):
            # Invalid peaks must be rejected even above the precursor cutoff.
            for peaks in (((250, intensity), (100, 1)), ((100, 1), (250, intensity))):
                with self.subTest(peaks=peaks), self.assertRaisesRegex(ValueError, "finite input intensities"):
                    mn.preprocess(mn.Spectrum(200, peaks))

    def test_frozen_spectrum_copies_mutable_inputs(self):
        peaks = [[100, 1]]
        spectrum = mn.Spectrum(200, peaks)
        peaks[0][1] = 2
        self.assertEqual(spectrum.peaks, ((100, 1),))


class KernelProperties(unittest.TestCase):
    @staticmethod
    def pairs():
        rng = random.Random(20261113)
        for _ in range(500):
            # Distinct integral intensities avoid tied assignment products and keep
            # C# argument-ordered Bonanza sums exact in the bitwise parity check.
            # The m/z grid also stays away from tolerance-window boundaries;
            # the boundary asymmetry has a separate anchor above.
            aa = sorted(rng.sample(range(10, 30), rng.randrange(2, 9)))
            bb = sorted(rng.sample(range(10, 30), rng.randrange(2, 9)))
            a = mn.Spectrum(400, [(x * 10, rng.randrange(1, 100000)) for x in aa])
            b = mn.Spectrum(414, [(x * 10, rng.randrange(1, 100000)) for x in bb])
            yield a, b

    @staticmethod
    def kernels():
        return [mn.msdial_moddot_legacy, mn.msdial_bonanza_v1,
                lambda a, b: mn.modified_cosine(a, b, tol=.025, assign="greedy"),
                lambda a, b: mn.modified_cosine(a, b, tol=.025, assign="optimal"),
                mn.entropy_li2021]

    def test_500_pairs_bounds_identity_symmetry_count_and_no_mutation(self):
        for a, b in self.pairs():
            snapshot = (a.precursor_mz, a.peaks, b.precursor_mz, b.peaks)
            for index, kernel in enumerate(self.kernels()):
                result = kernel(a, b)
                reverse = kernel(b, a)
                self.assertGreaterEqual(result.score, 0)
                self.assertLessEqual(result.score, 1 + 1e-12)
                self.assertLessEqual(result.matched_peak_count, min(len(a.peaks), len(b.peaks)))
                if index < 2:
                    self.assertEqual(result.score, reverse.score)
                else:
                    self.assertAlmostEqual(result.score, reverse.score, delta=1e-12)
                self.assertAlmostEqual(kernel(a, a).score, 1, delta=1e-12)
                self.assertEqual(kernel(a, b), result)  # call order has no effect
            self.assertEqual(snapshot, (a.precursor_mz, a.peaks, b.precursor_mz, b.peaks))
            greedy, optimal = [k(a, b) for k in self.kernels()[2:4]]
            self.assertGreaterEqual(optimal.score, greedy.score - 1e-12)

    def test_moddot_single_match_is_one_documented_defect(self):
        n = 0
        for a, b in self.pairs():
            result = mn.msdial_moddot_legacy(a, b)
            if result.matched_peak_count == 1:
                self.assertEqual(result.score, 1.0)
                n += 1
        self.assertGreater(n, 0)

    def test_unmatched_invariance_and_full_norm_penalties(self):
        for a, b in self.pairs():
            extra = mn.Spectrum(b.precursor_mz, b.peaks + ((9000, 100000),))
            for index, kernel in enumerate(self.kernels()[:4]):
                old, new = kernel(a, b), kernel(a, extra)
                self.assertEqual(old.matched_peak_count, new.matched_peak_count)
                if index == 0:
                    self.assertEqual(old.score, new.score)
                elif old.score > 0:
                    self.assertLess(new.score, old.score)

    def test_hungarian_against_exhaustive_assignments(self):
        rng = random.Random(20261113)
        for n in range(1, 6):
            for _ in range(20):
                weights = [[rng.randrange(30) for _ in range(n)] for _ in range(n)]
                actual = sum(weights[i][j] for i, j in mn._hungarian(weights))
                expected = max(sum(weights[i][j] for i, j in enumerate(p))
                               for p in itertools.permutations(range(n)))
                self.assertEqual(actual, expected)

    def test_component_assignment_rectangular_and_duplicate_edges(self):
        edges = [(0, 0, 10), (0, 0, 9), (0, 1, 9), (1, 0, 9), (3, 4, 5), (4, 4, 4)]
        accepted = mn._optimal(edges)
        self.assertEqual(sum(w for _, _, w in accepted), 23)
        self.assertEqual(len(accepted), 3)


class IdentityGrammar(unittest.TestCase):
    def test_all_canonical_ids_round_trip_and_dispatch(self):
        a = mn.Spectrum(250, ((100, 3), (150, 4)))
        for identity in mt.kernel_ids() + mt.kernel_ids(mn.MSDIAL_MN_V1) + mt.smoke_ids():
            base, params = mn.parse_kernel_id(identity)
            self.assertEqual(mn.parse_kernel_id(mn.format_kernel_id(base, params)), (base, params))
            self.assertEqual(mn.canonical_kernel_id(identity), identity)
            self.assertEqual(mn.score(identity, a, a).kernel_id, identity)
            keys = [part.split("=")[0] for part in identity.split("?", 1)[1].split("&")]
            self.assertEqual(keys, sorted(keys))

    def test_aliases_and_decimal_canonicalization(self):
        for alias in ("msdial.moddot.legacy", "msdial.bonanza.v1", "entropy.li2021"):
            self.assertEqual(mn.canonical_kernel_id(alias), mn.format_kernel_id(*mn.parse_kernel_id(alias)))
        value = "gnps.modified_cosine?tol=0.0250&prep=none&mzpow=-0&ipow=1.0&assign=greedy"
        self.assertEqual(mn.canonical_kernel_id(value), mt.kernel_ids()[3])
        params = {"prep": "none", "tol": 0.00000001}
        self.assertEqual(mn.format_kernel_id("msdial.bonanza.v1", params), "msdial.bonanza.v1?prep=none&tol=0.00000001")
        self.assertEqual(mn._literal(1e20), "100000000000000000000")

    def test_malformed_ids_raise(self):
        bad = ["unknown", "gnps.modified_cosine", "msdial.cosine_matched.legacy",
               "msdial.bonanza.v1?prep=none", "msdial.bonanza.v1?prep=none&tol=0&x=1",
               "msdial.bonanza.v1?prep=none&tol=1&tol=1", "msdial.bonanza.v1?prep=none&tol=NaN",
               "msdial.bonanza.v1?prep=none&tol=inf", "msdial.bonanza.v1?prep=none&tol=none",
               "msdial.bonanza.v1?prep=none&tol=-1", "msdial.bonanza.v1?prep=none&tol=0",
               "msdial.bonanza.v1?prep=none&tol=1e-2", "entropy.li2021?prep=none&bin=0.05&weighting=true"]
        for value in bad:
            with self.subTest(value=value), self.assertRaises(ValueError):
                mn.parse_kernel_id(value)
        for key, value in [("prep.topk", 1.5), ("prep.rel", True), ("prep.abs", math.inf),
                           ("prep.minmz", 20000)]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                mn.format_kernel_id("msdial.bonanza.v1", {**mt.prep_params(), "tol": .025, key: value})

    def test_assign_is_required(self):
        with self.assertRaises(TypeError):
            mn.modified_cosine(mn.Spectrum(200, ()), mn.Spectrum(200, ()), tol=.025)

    def test_python_310_syntax_and_imports(self):
        allowed = {"__future__", "dataclasses", "decimal", "math", "re", "similarity", "parsers",
                   "argparse", "hashlib", "itertools", "json", "pathlib", "sys", "time", ""}
        for module in (mn, mt):
            source = Path(module.__file__).read_text(encoding="utf-8")
            tree = ast.parse(source, feature_version=(3, 10))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    self.assertTrue(all(n.name in allowed for n in node.names))
                elif isinstance(node, ast.ImportFrom):
                    self.assertIn(node.module or "", allowed)


if __name__ == "__main__":
    unittest.main()
