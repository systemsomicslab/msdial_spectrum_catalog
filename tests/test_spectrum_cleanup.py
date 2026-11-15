"""Tests for the spectral-entropy peak cleanup."""

from __future__ import annotations

import itertools
import math
import random
import unittest

from msdial_spectrum_catalog.spectrum_cleanup import clean_by_entropy, spectral_entropy


def brute_force(peaks, ceiling):
    """Drop the weakest intensity level, one level at a time, until the entropy is at most the ceiling."""
    peaks = [p for p in peaks if p[1] > 0]
    while spectral_entropy(peaks) > ceiling and len({i for _, i in peaks}) > 1:  # one level is kept whole
        weakest = min(i for _, i in peaks)
        peaks = [p for p in peaks if p[1] != weakest]
    return peaks


class SpectralEntropyTest(unittest.TestCase):
    def test_values(self):
        self.assertEqual(spectral_entropy([]), 0.0)
        self.assertEqual(spectral_entropy([(100.0, 5)]), 0.0)
        self.assertAlmostEqual(spectral_entropy([(m, 7) for m in range(4)]), math.log(4))
        self.assertAlmostEqual(spectral_entropy([(1, 1), (2, 3)]), -(0.25 * math.log(0.25) + 0.75 * math.log(0.75)))
        self.assertAlmostEqual(spectral_entropy([(1, 1), (2, 0)]), 0.0)  # zero intensity carries nothing

    def test_a_low_entropy_spectrum_is_unchanged(self):
        peaks = [(50.0, 1), (60.0, 1000), (70.0, 3)]
        result = clean_by_entropy(peaks, 4.5)
        self.assertEqual(result.peaks, peaks)
        self.assertEqual(result.removed, 0)
        self.assertIsNone(result.intensity_threshold)

    def test_noise_is_dropped_and_order_is_kept(self):
        signal = [(100.0, 1000), (150.0, 600), (200.0, 300)]
        noise = [(100.5 + k, 20 + k % 41) for k in range(300)]
        peaks = sorted(signal + noise)  # m/z order, as written
        result = clean_by_entropy(peaks, 4.5)
        self.assertGreater(result.entropy_before, 4.5)
        self.assertLessEqual(result.entropy_after, 4.5)
        self.assertTrue(set(signal) <= set(result.peaks))
        self.assertEqual(result.peaks, sorted(result.peaks))
        self.assertEqual(result.removed, len(peaks) - len(result.peaks))
        self.assertTrue(all(i >= result.intensity_threshold for _, i in result.peaks))

    def test_equal_intensities_go_together(self):
        peaks = [(float(m), 10) for m in range(200)] + [(500.0, 100)]
        result = clean_by_entropy(peaks, 4.5)
        self.assertGreater(result.entropy_before, 4.5)
        self.assertEqual(result.peaks, [(500.0, 100)])  # the 200 equal ions go together
        alone = clean_by_entropy([(float(m), 10) for m in range(200)], 4.5)
        self.assertEqual(len(alone.peaks), 200)  # one level above the ceiling cannot be split; it is kept

    def test_matches_the_one_level_at_a_time_rule(self):
        rng = random.Random(7)
        for _ in range(200):
            n = rng.randint(1, 400)
            peaks = [(100.0 + k, rng.choice([rng.randint(1, 20), rng.randint(1, 100000)])) for k in range(n)]
            ceiling = rng.choice([1.0, 2.5, 4.5])
            result = clean_by_entropy(peaks, ceiling)
            self.assertEqual(result.peaks, brute_force(peaks, ceiling) if result.removed else peaks)
            if len({i for _, i in result.peaks}) > 1:
                self.assertLessEqual(result.entropy_after, ceiling + 1e-9)

    def test_dropping_the_weakest_ion_never_raises_the_entropy(self):
        rng = random.Random(11)
        for _ in range(300):
            values = sorted((rng.uniform(0.01, 1000) for _ in range(rng.randint(2, 30))), reverse=True)
            entropies = [spectral_entropy([(0, v) for v in values[:k]]) for k in range(1, len(values) + 1)]
            for smaller, larger in itertools.pairwise(entropies):
                self.assertLessEqual(smaller, larger + 1e-12)

    def test_the_ceiling_must_be_positive(self):
        with self.assertRaises(ValueError):
            clean_by_entropy([(1, 1)], 0)


if __name__ == "__main__":
    unittest.main()
