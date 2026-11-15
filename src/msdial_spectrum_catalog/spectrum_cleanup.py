"""MS/MS peak cleanup by spectral entropy.

Spectral entropy (Li et al., Nat. Methods 18, 1524-1531 (2021)) is the Shannon entropy of a spectrum's
intensities, normalised to sum 1, with the natural logarithm: S = -sum(p_i ln p_i). Library spectra mostly
have S of 1 to 3, and S above about 5 points to many low-intensity (noise) ions. The cleanup drops the
weakest ions until S is at most a ceiling (4.5 by convention here), so noisy DDA spectra reach the
in-silico tools with their informative ions only.

Dropping the weakest ion always lowers S (adding an ion no stronger than every other one raises it), so
"remove the weakest ions until S <= ceiling" keeps the largest set of strongest ions whose entropy is at
most the ceiling. Ions of equal intensity are kept or dropped together, so the result does not depend on
peak order. A spectrum already at or below the ceiling is returned unchanged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

DEFAULT_MAX_ENTROPY = 4.5


def spectral_entropy(peaks) -> float:
    """S = -sum(p ln p) over positive intensities normalised to sum 1; 0 for an empty spectrum."""
    intensities = [float(intensity) for _, intensity in peaks if float(intensity) > 0]
    total = sum(intensities)
    if total <= 0:
        return 0.0
    return -sum((value / total) * math.log(value / total) for value in intensities)


@dataclass(frozen=True)
class CleanupResult:
    peaks: list
    entropy_before: float
    entropy_after: float
    removed: int
    intensity_threshold: float | None  # the weakest intensity kept, when ions were removed


def clean_by_entropy(peaks, max_entropy: float = DEFAULT_MAX_ENTROPY) -> CleanupResult:
    """Drop the weakest ions (ties together) until the spectral entropy is at most max_entropy.

    The kept peaks stay in their input order. Peaks with zero or negative intensity carry no information
    and are dropped only when a cleanup happens.
    """
    if max_entropy <= 0:
        raise ValueError("max_entropy must be positive")
    peaks = list(peaks)
    before = spectral_entropy(peaks)
    if before <= max_entropy:
        return CleanupResult(peaks, before, before, 0, None)
    # Entropy of the k strongest ions: S_k = ln T_k - (sum_{i<=k} I_i ln I_i) / T_k, with T_k their sum.
    intensities = sorted((float(i) for _, i in peaks if float(i) > 0), reverse=True)
    total = weighted = 0.0
    threshold = intensities[0]
    index = 0
    while index < len(intensities):
        level = intensities[index]
        group_total = group_weighted = 0.0
        while index < len(intensities) and intensities[index] == level:
            group_total += level
            group_weighted += level * math.log(level)
            index += 1
        candidate_total, candidate_weighted = total + group_total, weighted + group_weighted
        if math.log(candidate_total) - candidate_weighted / candidate_total > max_entropy:
            break
        total, weighted, threshold = candidate_total, candidate_weighted, level
    if total == 0:
        # Even the strongest tie group alone is above the ceiling (many ions of one intensity): keep it.
        threshold = intensities[0]
    kept = [peak for peak in peaks if float(peak[1]) >= threshold]
    return CleanupResult(kept, before, spectral_entropy(kept), len(peaks) - len(kept), threshold)
