"""Optional reference cross-check; matchms is never a runtime dependency."""
import importlib.util
import inspect
import json
import math
from pathlib import Path
import random
import unittest

from msdial_spectrum_catalog import similarity_mn as mn


@unittest.skipUnless(importlib.util.find_spec("matchms"), "matchms is not installed; reference check deferred")
class MatchmsParity(unittest.TestCase):
    def setUp(self):
        import numpy as np
        from matchms import Spectrum
        from matchms.similarity import ModifiedCosineGreedy, ModifiedCosineHungarian
        self.np, self.Spectrum = np, Spectrum
        self.classes = {"greedy": ModifiedCosineGreedy, "optimal": ModifiedCosineHungarian}

    def reference_spectrum(self, spectrum):
        peaks = sorted(spectrum.peaks, key=lambda p: p[0])
        kwargs = {"metadata": {"precursor_mz": spectrum.precursor_mz},
                  "mz": self.np.array([p[0] for p in peaks], dtype=float),
                  "intensities": self.np.array([p[1] for p in peaks], dtype=float)}
        if "metadata_harmonization" in inspect.signature(self.Spectrum).parameters:
            kwargs["metadata_harmonization"] = False
        return self.Spectrum(**kwargs)

    def compare(self, a, b, params):
        # MN zero/zero weights are outside the finite matchms/SciPy domain.
        # Raw all-zero spectra are still checked (both give zero norm).
        if not all(math.isfinite(i) for _, i in a.peaks + b.peaks):
            return
        cls = self.classes[params["assign"]]
        kwargs = {"tolerance": params["tol"], "mz_power": params["mzpow"],
                  "intensity_power": params["ipow"]}
        signature = inspect.signature(cls).parameters
        if "noise_cutoff" in signature:
            kwargs["noise_cutoff"] = 0
        if "remove_precursor" in signature:
            kwargs["remove_precursor"] = False
        expected = cls(**kwargs).pair(self.reference_spectrum(a), self.reference_spectrum(b))
        result = mn.modified_cosine(a, b, tol=params["tol"], assign=params["assign"],
                                    mz_power=params["mzpow"], intensity_power=params["ipow"])
        self.assertAlmostEqual(result.score, float(expected["score"]), delta=1e-9)
        if params["assign"] == "greedy":
            self.assertEqual(result.matched_peak_count, int(expected["matches"]))

    def test_all_fixture_pairs(self):
        root = Path(__file__).resolve().parents[1] / "similarity-fixtures"
        for path in sorted(root.glob("kernels*.json")):
            for case in json.loads(path.read_bytes())["cases"]:
                for expected in case["expected"]:
                    if not expected["kernel_id"].startswith("gnps.modified_cosine?"):
                        continue
                    _, params = mn.parse_kernel_id(expected["kernel_id"])
                    a, b = [mn.Spectrum(case[k]["precursor_mz"], case[k]["peaks"]) for k in ("a", "b")]
                    if params["prep"] == "mn":
                        prep = mn.MnPrep(**{k[5:]: v for k, v in params.items() if k.startswith("prep.")})
                        a, b = mn.preprocess(a, prep), mn.preprocess(b, prep)
                    with self.subTest(case=case["id"], identity=expected["kernel_id"]):
                        self.compare(a, b, params)

    def test_200_random_pairs(self):
        rng = random.Random(20261113)
        for index in range(200):
            precursors = [300 + rng.randrange(6) * 14 for _ in range(2)]
            spectra = [mn.Spectrum(p, [(mz + rng.uniform(-.01, .01), rng.uniform(.01, 1))
                       for mz in sorted(rng.sample(range(80, 250), rng.randrange(1, 15)))])
                       for p in precursors]
            for assign in ("greedy", "optimal"):
                for mzpow, ipow in ((0, 1), (1, .5)):
                    with self.subTest(index=index, assign=assign, mzpow=mzpow):
                        self.compare(*spectra, {"tol": .025, "assign": assign, "mzpow": mzpow, "ipow": ipow})


if __name__ == "__main__":
    unittest.main()
