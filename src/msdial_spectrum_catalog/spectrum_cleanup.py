"""MS/MS peak cleanup by spectral entropy.

Spectral entropy (Li et al., Nat. Methods 18, 1524-1531 (2021)) is the Shannon entropy of a spectrum's
intensities, normalised to sum 1, with the natural logarithm: S = -sum(p_i ln p_i). The cleanup drops the
weakest ions until S is at most a ceiling, so noisy DDA spectra reach the in-silico tools with fewer
low-intensity ions.

Scale. Li et al. compute S after their own preprocessing (peaks within 0.05 Da merged, ions below 1% of
the base peak and the precursor region removed); their ranges (library spectra mostly 1 to 3, above about
5 suggesting noise) describe that scale. Here S is computed on the spectrum as exported, without that
preprocessing, so split peaks and weak ions raise it and those ranges are not a calibration of the ceiling.
The ceiling (4.5 by convention here) is an empirical choice to be judged by its effect on annotation.

The cut is relative. The threshold depends on the spectrum's own intensity distribution, not on an
absolute noise floor: in a flat spectrum with many similar ions it can remove ions at a sizeable fraction
of the base peak. The result reports the kept threshold and the largest removed ion relative to the base
peak, so an evaluation can be stratified by them.

Rule. Dropping the weakest ion never raises S (an ion no stronger than every other one adds entropy), so
the cleanup keeps the strongest ions, adding one intensity level at a time while S stays at or below the
ceiling. Ions of equal intensity are kept or dropped together, so the result does not depend on peak
order: it is the largest set of strongest ions, closed under equal intensities, whose entropy is at most
the ceiling. When the strongest intensity level alone is above the ceiling (many ions of one intensity),
no such set exists: that level is kept, the rest is dropped, and the result says so (status
"above_ceiling"). Entries with zero or negative intensity carry no entropy and are left as they are, so a
spectrum the cleanup touches differs from the original only by the positive ions it dropped. A spectrum
already at or below the ceiling is returned unchanged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

DEFAULT_MAX_ENTROPY = 4.5
UNCHANGED, CLEANED, ABOVE_CEILING = "unchanged", "cleaned", "above_ceiling"


def check_max_entropy(max_entropy) -> float:
    """Return max_entropy as a float, or raise ValueError unless it is a finite positive number."""
    try:
        value = float(max_entropy)
    except (TypeError, ValueError):
        raise ValueError(f"max_entropy must be a finite positive number, not {max_entropy!r}") from None
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"max_entropy must be a finite positive number, not {max_entropy!r}")
    return value


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
    removed: int  # positive ions dropped
    intensity_threshold: float | None  # the weakest intensity kept, when ions were dropped
    base_peak: float | None  # the strongest intensity, when there is one
    largest_removed: float | None  # the strongest intensity dropped, when ions were dropped
    status: str  # UNCHANGED, CLEANED or ABOVE_CEILING


def clean_by_entropy(peaks, max_entropy: float = DEFAULT_MAX_ENTROPY) -> CleanupResult:
    """Drop the weakest positive ions (equal intensities together) until the spectral entropy is at most
    max_entropy. The kept peaks stay in their input order; non-positive entries are kept as they are."""
    max_entropy = check_max_entropy(max_entropy)
    peaks = list(peaks)
    before = spectral_entropy(peaks)
    intensities = sorted((float(i) for _, i in peaks if float(i) > 0), reverse=True)
    base = intensities[0] if intensities else None
    if before <= max_entropy:
        return CleanupResult(peaks, before, before, 0, None, base, None, UNCHANGED)
    # Entropy of the k strongest ions: S_k = ln T_k - (sum_{i<=k} I_i ln I_i) / T_k, with T_k their sum.
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
    status = CLEANED if total > 0 else ABOVE_CEILING  # ABOVE_CEILING: the strongest level alone is too much
    kept = [peak for peak in peaks if float(peak[1]) >= threshold or float(peak[1]) <= 0]
    dropped = [float(i) for _, i in peaks if 0 < float(i) < threshold]
    return CleanupResult(kept, before, spectral_entropy(kept), len(dropped), threshold, base,
                         max(dropped) if dropped else None, status)
