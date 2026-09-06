"""Spectral similarity measures for reference-versus-reference comparison.

These are ports of the MS-DIAL kernel in
`MsdialWorkbench/src/Common/CommonStandard/Algorithm/Scoring/MsScanMatching.cs`, with one deliberate
departure documented below.

TWO CONVENTIONS, NEVER MIXED. MS-DIAL's GetSimpleDotProduct and GetWeightedDotProduct return
`covariance^2 / scalarM / scalarR`, which is the SQUARE of a cosine, and store it in
`MsScanMatchResult.Squared*`. Those fields are read only for threshold comparison; both text exporters
write the non-squared computed getters, so the values in a `.mdpeak` or `.mdalign` column are plain
cosines. Everything here reports cosine scale and names the squared value explicitly where it is kept.

SYMMETRY. MS-DIAL's weighted dot product is asymmetric in two places: `peakCountPenalty` is derived from
`lSpectrumCounter`, which counts peaks on the reference side only, and the 0.01 low-intensity cutoff in
its final loop is applied to the measured side only. Asymmetry is acceptable for query-versus-library
scoring, where one side really is the reference, but an ambiguity class asserts a mutual relation, so
`sim(A, B)` must equal `sim(B, A)`. `weighted_cosine` is therefore a symmetric variant: the intensity
cutoff applies only where both sides are below it, and no peak-count penalty is applied. The penalty is
redundant here because the caller's admissibility gate already excludes low-information spectra
explicitly, and it would otherwise make the relation depend on argument order.
`msdial_weighted_dot_product` is the faithful asymmetric port, kept so a class edge can be compared
against an MS-DIAL annotation score.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

Peaks = list[tuple[float, float]]

# MsScanMatching.cs derives peakCountPenalty from the count of reference peaks above 0.1 relative
# intensity. Retained for the faithful port only.
_PEAK_COUNT_PENALTY = {1: 0.75, 2: 0.88, 3: 0.94, 4: 0.97}

# The intensity floor applied inside GetWeightedDotProduct's accumulation loop.
_INTENSITY_CUTOFF = 0.01

# GetSpectralEntropySimilarity bins on a fixed grid, int(mass / bin), rather than by sliding window.
# Two peaks 0.002 apart can therefore fall in different frames. That is MS-DIAL's behaviour and is
# reproduced rather than improved on, so entropy values stay comparable with MS-DIAL's own.
_ENTROPY_BIN = 0.05


@dataclass(frozen=True)
class SpectrumComparison:
    """One reference-versus-reference comparison, with every convention named."""

    comparable: bool
    reason: str | None = None
    weighted_cosine: float | None = None
    entropy_similarity: float | None = None
    matched_peak_count: int = 0
    informative_peak_count_a: int = 0
    informative_peak_count_b: int = 0
    msdial_weighted_squared: float | None = None
    msdial_peak_count_penalty: float | None = None


def normalize_to_base_peak(peaks: Peaks) -> Peaks:
    """Scale intensities so the base peak is 1.0, dropping non-positive peaks."""
    positive = [(mz, intensity) for mz, intensity in peaks if intensity > 0]
    if not positive:
        return []
    base = max(intensity for _, intensity in positive)
    return [(mz, intensity / base) for mz, intensity in positive]


def informative_peak_count(peaks: Peaks, relative_floor: float) -> int:
    """Count peaks at or above a relative-intensity floor, after base-peak normalization."""
    return sum(1 for _, intensity in normalize_to_base_peak(peaks) if intensity >= relative_floor)


def _align(peaks_a: Peaks, peaks_b: Peaks, tolerance: float) -> list[tuple[float, float, float]]:
    """Pair peaks within a +/- tolerance window, returning (mz, intensity_a, intensity_b) frames.

    Both inputs must be sorted by m/z. A frame carries the summed intensity of every peak of each
    spectrum inside the window, so a frame with one side zero is an unmatched peak. The walk is
    symmetric in the two arguments, unlike MS-DIAL's, whose index advancement differs between the
    measured and reference loops.
    """
    frames: list[tuple[float, float, float]] = []
    i = j = 0
    while i < len(peaks_a) or j < len(peaks_b):
        if i >= len(peaks_a):
            centre = peaks_b[j][0]
        elif j >= len(peaks_b):
            centre = peaks_a[i][0]
        else:
            centre = min(peaks_a[i][0], peaks_b[j][0])
        lower, upper = centre - tolerance, centre + tolerance
        sum_a = 0.0
        while i < len(peaks_a) and peaks_a[i][0] < upper:
            if peaks_a[i][0] >= lower:
                sum_a += peaks_a[i][1]
            i += 1
        sum_b = 0.0
        while j < len(peaks_b) and peaks_b[j][0] < upper:
            if peaks_b[j][0] >= lower:
                sum_b += peaks_b[j][1]
            j += 1
        frames.append((centre, sum_a, sum_b))
    return frames


def _restrict(peaks: Peaks, mass_begin: float, mass_end: float) -> Peaks:
    return sorted((mz, intensity) for mz, intensity in peaks if mass_begin <= mz <= mass_end)


def weighted_cosine(
    peaks_a: Peaks,
    peaks_b: Peaks,
    tolerance: float,
    *,
    mass_begin: float = 0.0,
    mass_end: float = float("inf"),
    intensity_cutoff: float = _INTENSITY_CUTOFF,
) -> tuple[float | None, int]:
    """Symmetric m/z-weighted cosine on cosine scale, plus the matched frame count.

    Weighting follows MS-DIAL: the covariance term is sqrt(intensity_a * intensity_b) * m/z and each
    scalar term is intensity * m/z, which is a cosine between vectors weighted by sqrt(intensity * m/z).
    Returns None when either spectrum has no usable peak in range.
    """
    a = normalize_to_base_peak(_restrict(peaks_a, mass_begin, mass_end))
    b = normalize_to_base_peak(_restrict(peaks_b, mass_begin, mass_end))
    if not a or not b:
        return None, 0
    covariance = scalar_a = scalar_b = 0.0
    matched = 0
    for mz, intensity_a, intensity_b in _align(a, b, tolerance):
        # Symmetric noise floor: drop a frame only when NEITHER spectrum has signal in it. A frame where
        # one side is strong and the other is noise is exactly a discriminating fragment and is kept.
        if intensity_a < intensity_cutoff and intensity_b < intensity_cutoff:
            continue
        covariance += math.sqrt(intensity_a * intensity_b) * mz
        scalar_a += intensity_a * mz
        scalar_b += intensity_b * mz
        if intensity_a > 0 and intensity_b > 0:
            matched += 1
    if scalar_a <= 0 or scalar_b <= 0:
        return None, matched
    # Divide by the product, not successively: (x / a) / b and (x / b) / a are not always the same
    # double, which would make the relation asymmetric by one ulp and could flip a comparison sitting
    # exactly on the threshold. IEEE-754 multiplication is commutative, so the product is exact here.
    squared = covariance * covariance / (scalar_a * scalar_b)
    # Floating-point accumulation can still push a perfect self-comparison a hair above 1.
    return math.sqrt(min(squared, 1.0)), matched


def simple_cosine(
    peaks_a: Peaks,
    peaks_b: Peaks,
    tolerance: float,
    *,
    mass_begin: float = 0.0,
    mass_end: float = float("inf"),
) -> float | None:
    """Symmetric unweighted cosine of the sqrt-intensity vectors, on cosine scale."""
    a = normalize_to_base_peak(_restrict(peaks_a, mass_begin, mass_end))
    b = normalize_to_base_peak(_restrict(peaks_b, mass_begin, mass_end))
    if not a or not b:
        return None
    covariance = scalar_a = scalar_b = 0.0
    for _, intensity_a, intensity_b in _align(a, b, tolerance):
        covariance += math.sqrt(intensity_a * intensity_b)
        scalar_a += intensity_a
        scalar_b += intensity_b
    if scalar_a <= 0 or scalar_b <= 0:
        return None
    return math.sqrt(min(covariance * covariance / (scalar_a * scalar_b), 1.0))


def _binned(peaks: Peaks, bin_width: float) -> Peaks:
    """Fixed-grid binning as in SpectrumHandler.GetBinnedSpectrum: sum per frame, keep the apex m/z."""
    frames: dict[int, list[tuple[float, float]]] = {}
    for mz, intensity in peaks:
        frames.setdefault(int(mz / bin_width), []).append((mz, intensity))
    result = []
    for members in frames.values():
        # Sort before reducing so the summation order, and the apex chosen among tied intensities, do
        # not depend on the order the peaks arrived in.
        ordered = sorted(members)
        apex = max(ordered, key=lambda peak: (peak[1], peak[0]))[0]
        result.append((apex, math.fsum(intensity for _, intensity in ordered)))
    return result


def _combined(peaks_a: Peaks, peaks_b: Peaks, bin_width: float) -> Peaks:
    """SpectrumHandler.GetCombinedSpectrum: merge both spectra per frame and halve the summed intensity."""
    frames: dict[int, list[tuple[float, float]]] = {}
    for mz, intensity in list(peaks_a) + list(peaks_b):
        frames.setdefault(int(mz / bin_width), []).append((mz, intensity))
    result = []
    for members in frames.values():
        # Same reason as _binned, and it matters more here: this helper concatenates both spectra, so
        # without an order-independent reduction the entropy similarity would depend on which spectrum
        # was passed first, and the relation must be symmetric.
        ordered = sorted(members)
        apex = max(ordered, key=lambda peak: (peak[1], peak[0]))[0]
        result.append((apex, math.fsum(intensity for _, intensity in ordered) * 0.5))
    return result


def spectral_entropy(peaks: Peaks) -> float:
    """Shannon entropy in bits of the intensity distribution, ignoring non-positive peaks."""
    # Drop non-positive peaks BEFORE totalling. Dividing by a signed total lets a share exceed 1 and
    # returns a negative "entropy", which is not a degenerate value but a wrong one.
    positive = sorted(intensity for _, intensity in peaks if intensity > 0)
    total = math.fsum(positive)
    if total <= 0:
        return 0.0
    entropy = 0.0
    for intensity in positive:
        share = intensity / total
        entropy -= share * math.log2(share)
    return entropy


def entropy_similarity(peaks_a: Peaks, peaks_b: Peaks, bin_width: float = _ENTROPY_BIN) -> float | None:
    """Li-Fiehn spectral entropy similarity, matching MsScanMatching.GetSpectralEntropySimilarity.

    The published formula divides by ln(4); MS-DIAL computes entropy in bits and multiplies by 0.5,
    which is the same normalization. Symmetric by construction.
    """
    if not peaks_a or not peaks_b:
        return None
    total_a = sum(intensity for _, intensity in peaks_a)
    total_b = sum(intensity for _, intensity in peaks_b)
    if total_a <= 0 or total_b <= 0:
        return None
    unit_a = [(mz, intensity / total_a) for mz, intensity in peaks_a]
    unit_b = [(mz, intensity / total_b) for mz, intensity in peaks_b]
    merged = spectral_entropy(_combined(unit_a, unit_b, bin_width))
    entropy_a = spectral_entropy(_binned(peaks_a, bin_width))
    entropy_b = spectral_entropy(_binned(peaks_b, bin_width))
    # Sum the two single-spectrum entropies before subtracting, so argument order cannot change the
    # result by one ulp. The bracket is twice the Jensen-Shannon divergence in bits, which lies in
    # [0, 1], so the similarity is provably in [0, 1] and any excursion outside it is accumulation
    # noise -- clamp rather than report a negative similarity.
    similarity = 1.0 - (2.0 * merged - (entropy_a + entropy_b)) * 0.5
    return min(max(similarity, 0.0), 1.0)


def msdial_weighted_dot_product(
    measured: Peaks,
    reference: Peaks,
    tolerance: float,
    *,
    mass_begin: float = 0.0,
    mass_end: float = float("inf"),
) -> tuple[float | None, float]:
    """Faithful port of GetWeightedDotProduct: returns (squared value, peak-count penalty).

    Asymmetric on purpose, matching MS-DIAL: the penalty counts reference peaks above 0.1 relative
    intensity, and the 0.01 intensity cutoff is applied to the measured side only. Use it to compare a
    class edge against an MS-DIAL annotation score, never as the class relation itself.
    """
    a = normalize_to_base_peak(_restrict(measured, mass_begin, mass_end))
    b = normalize_to_base_peak(_restrict(reference, mass_begin, mass_end))
    if not a or not b:
        return None, 1.0
    frames = _align(a, b, tolerance)
    reference_peaks_above_floor = sum(1 for _, _, intensity_b in frames if intensity_b > 0.1)
    penalty = _PEAK_COUNT_PENALTY.get(reference_peaks_above_floor, 1.0)
    covariance = scalar_measured = scalar_reference = 0.0
    for mz, intensity_a, intensity_b in frames:
        if intensity_a < _INTENSITY_CUTOFF:
            continue
        covariance += math.sqrt(intensity_a * intensity_b) * mz
        scalar_measured += intensity_a * mz
        scalar_reference += intensity_b * mz
    if scalar_measured <= 0 or scalar_reference <= 0:
        return 0.0, penalty
    return covariance * covariance / scalar_measured / scalar_reference * penalty, penalty


def compare(
    peaks_a: Peaks,
    peaks_b: Peaks,
    *,
    tolerance: float,
    minimum_informative_peaks: int,
    minimum_matched_peaks: int,
    relative_floor: float,
    mass_begin: float = 0.0,
    mass_end: float = float("inf"),
) -> SpectrumComparison:
    """Compare two reference spectra, applying the admissibility gate before reporting any similarity.

    The gate is not a refinement. Half of the records in the lab's public positive library carry fewer
    than ten peaks, and two spectra with two peaks each will look identical without that meaning
    anything. A pair that fails the gate is reported as not comparable with a reason, which is a distinct
    outcome from being indistinguishable and must never be merged into a class.
    """
    informative_a = informative_peak_count(_restrict(peaks_a, mass_begin, mass_end), relative_floor)
    informative_b = informative_peak_count(_restrict(peaks_b, mass_begin, mass_end), relative_floor)
    if informative_a < minimum_informative_peaks or informative_b < minimum_informative_peaks:
        return SpectrumComparison(
            comparable=False,
            reason=(
                f"insufficient_informative_peaks: {informative_a} and {informative_b} peaks at or above "
                f"{relative_floor:g} relative intensity, minimum {minimum_informative_peaks}"
            ),
            informative_peak_count_a=informative_a,
            informative_peak_count_b=informative_b,
        )
    # The cutoff shares the gate's floor so the two cannot drift apart: a peak the gate counts as
    # informative is a peak the similarity must see.
    cosine, matched = weighted_cosine(
        peaks_a, peaks_b, tolerance, mass_begin=mass_begin, mass_end=mass_end,
        intensity_cutoff=relative_floor,
    )
    if cosine is None:
        return SpectrumComparison(
            comparable=False,
            reason="no_usable_peaks_in_range",
            matched_peak_count=matched,
            informative_peak_count_a=informative_a,
            informative_peak_count_b=informative_b,
        )
    if matched < minimum_matched_peaks:
        return SpectrumComparison(
            comparable=False,
            reason=f"insufficient_matched_peaks: {matched}, minimum {minimum_matched_peaks}",
            matched_peak_count=matched,
            informative_peak_count_a=informative_a,
            informative_peak_count_b=informative_b,
        )
    squared, penalty = msdial_weighted_dot_product(
        peaks_a, peaks_b, tolerance, mass_begin=mass_begin, mass_end=mass_end
    )
    return SpectrumComparison(
        comparable=True,
        weighted_cosine=cosine,
        entropy_similarity=entropy_similarity(peaks_a, peaks_b),
        matched_peak_count=matched,
        informative_peak_count_a=informative_a,
        informative_peak_count_b=informative_b,
        msdial_weighted_squared=squared,
        msdial_peak_count_penalty=penalty,
    )
