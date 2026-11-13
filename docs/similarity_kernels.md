# Molecular-networking similarity kernels (version 1)

`similarity_mn` uses only the standard library. Inputs are immutable `Spectrum`
objects; low-level scoring functions consume already processed peaks and return
`KernelResult`. `score(identity, a, b)` applies the named preprocessing first.
All scores are unrounded. Stored method_version is `similarity_mn/1`.

## Definitions and source parity

The source references below are in MsdialWorkbench (reviewed at the supplied
source revision, 7bff7bf9c). Paths are relative to that repository.

* `msdial.moddot.legacy`: sum of matched intensity products divided by the
  product of matched-only Euclidean norms. `Algorithm/Scoring/MsScanMatching.cs`
  882-904. No matches gives 0; a zero denominator gives NaN. A single positive
  match scores 1 regardless of unmatched signal. This defect is why it is legacy.
* `msdial.bonanza.v1`: P / (P + Ua + Ub), where P is the matched product sum
  and U is each argument's unmatched squared-intensity sum; lines 906-927.
  It is a ratio (`bonanza_ratio`), not a cosine.
* `msdial.cosine_matched.legacy`: product-ion matches only, matched-only norms,
  with zero denominator mapped to 0. `msdial_all_v1` also exposes Bonanza and
  the zero-guarded ModDot of lines 929-964. Cosine conventions are
  `cosine_matched_only`; no squared scores are stored.
* `gnps.modified_cosine`: sum of assigned weighted products divided by the
  full-spectrum norms. Weight per peak is mz**mzpow * intensity**ipow.
  `assign=greedy` is primary; `assign=optimal` is a sensitivity check. The
  reference is matchms 0.34.0 `ModifiedCosineGreedy`, `ModifiedCosineHungarian`
  and `spectrum_similarity_functions.find_matches/score_best_matches`.
  The convention is `cosine`. No implicit filtering is applied by this kernel.
* `entropy.li2021`: the existing unweighted entropy_similarity, equal to
  1 - (2 H((a+b)/2) - H(a) - H(b))/2 in bits, using fixed bins. No precursor
  shift; empty inputs return None rather than C#'s -1. Convention:
  `entropy_similarity`. This wrapper does not change similarity.py.

The non-mutating MS-DIAL matcher ports lines 991-1046. The smaller precursor is
`e`; equal precursors make the second argument `e`. The outer loop visits `r`
in list order and chooses the unused `e` with strictly smallest intensity
difference; ties keep the first index. Product matching is tested first using
`abs(e_mz-r_mz) < tol`, then neutral loss using
`abs((r_prec-e_prec + e_mz)-r_mz) < tol`, preserving expression order.
Local matched flags replace C# input mutation. Only absolute Da tolerance is
exposed. Bonanza sums unmatched squares in argument order as C# does.

The modified cosine sorts peaks stably by m/z. It uses inclusive windows,
shift = Pa - Pb added to b's m/z, direct candidates followed by shifted
candidates, and no shifted candidates if abs(shift) <= tol. Greedy sorts
stably ascending by weight then reverses (including ties). Optimal de-duplicates
edges and uses a pure Python Hungarian solver per connected component. Tied
optimal assignments may have different matched counts; only their score is a
parity contract. Zero full norm returns 0. Scores are not clamped.

## Preprocessing

MN defaults (`MSDIAL_MN_V1`) follow MsScanMatching.cs 1048-1093 and
`Algorithm/Function/SpectrumHandler.cs` 129-162, 63-70 and 26-30:

1. Take the base intensity before filtering. Keep mz in [0,10000], mz <= P+0.05,
   I >= 0 and I >= base*1*0.01. The precursor stays; isotope removal is off.
2. Optional precursor window removes abs(mz-P) <= precwin.
3. Optional top-k retains p when fewer than k other peaks within +/-topkwin
   have strictly greater intensity. Equal intensities all survive. This is our
   reading of the GNPS text, not a port of GNPS source.
4. Stable m/z order, frame = int(mz/frame), sequential intensity sums. Apex is
   the first strictly greater intensity (lower m/z on ties).
   `Algorithm/Function/FrameBinnedSpectrum.cs` 27-54.
5. Scale by (I**pow)/(max_I**pow)*scale; defaults pow=0.5, scale=100.
   IEEE zero/zero is preserved as NaN rather than raising ZeroDivisionError.

The precursor-window and sliding top-k switches are Python extensions and run
before framing. With frame=0, framing is disabled. `prep=none` sorts only.
The default callers are `MoleculerNetworking.cs` 94-100; MN parameters are
`ParameterBase.cs` 1591-1611. `msdial_edge_passes` follows the gate order in
`MoleculerNetworking.cs` 177-198: precursor difference >400, matched <2,
score <50*0.01. It keeps IEEE NaN comparison behavior and omits GUI rounding
at line 166. None (entropy not comparable) fails the gate.

## Canonical identities

Format: `<base>?<ASCII-sorted key=value pairs joined by &>`. Every applicable
key is required; duplicate, unknown, missing and invalid keys are rejected.
Numbers are finite plain decimal shortest round-trip values without trailing
zeros, decimal point or exponent; negative zero becomes 0. None is `none` and
booleans have lowercase spelling (no current key accepts booleans).

Prep is `none` or `mn`; mn requires all of `prep.abs`, `prep.above`,
`prep.frame`, `prep.maxmz`, `prep.minmz`, `prep.pow`, `prep.precwin`, `prep.rel`,
`prep.scale`, `prep.topk`, `prep.topkwin`. Only prep.above accepts none.
The kernel parameters are:

| Base | Parameters | Bare alias |
| --- | --- | --- |
| msdial.moddot.legacy | tol | MN defaults, tol=0.025 |
| msdial.bonanza.v1 | tol | MN defaults, tol=0.025 |
| msdial.cosine_matched.legacy | tol | none |
| gnps.modified_cosine | assign, ipow, mzpow, tol | none |
| entropy.li2021 | bin, weighting=none | bin=0.05, prep=none |

`assign` is required, with values greedy or optimal. Tolerance and bin width
must be positive; all other numerical controls are nonnegative, minmz <= maxmz,
and topk is an integer. `format_kernel_id`, `parse_kernel_id` and
`canonical_kernel_id` share this validation. Comparison inclusivity belongs to
the base, not a parameter. The full canonical string is the storage method;
`kernel-ids.json` pins examples shared with C#.

## Fixtures and operational checks

`similarity-fixtures/README.md` defines the versioned parity files, public
provenance, bounded input precision and deterministic regeneration. Optimal,
precwin and topk expectations are tagged Python only. Other expectations are
for Python and C#. Fixtures include raw and MN-processed comparisons.

`python -m msdial_spectrum_catalog.similarity_mn_tools smoke --msp PATH`
retains only precursor, peaks and COMMENT PEAKID during parsing, emitting only
aggregates. K1/K2 are default ModDot/Bonanza; K3 is greedy modified cosine with
MN prep; K4 uses frame=0, pow=0.5, precwin=17, topk=6, topkwin=50 and otherwise
MN defaults. K4's square root is an explicit assumption, not verified GNPS
behavior. K5 is entropy with prep none. Quantiles use linear interpolation;
null/nonfinite scores are excluded from quantiles and threshold fractions and
counted separately. Matched-one fractions use all eligible pairs; entropy
reports 0 matched because that wrapper does not define peak assignment.

The optional matchms test makes filtering explicit where supported. matchms
is not a runtime dependency and is never installed by this module.
