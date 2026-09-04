# Ambiguity classes

## The question

> Given that a query matched reference library entry A, which other entries would be indistinguishable
> from A under the study's measurement conditions?

That set is A's ambiguity class. If entry B's spectrum is itself highly similar to A's, then matching A
did not distinguish A from B, and reporting "A" alone overstates the result; the honest report is
"A or B". This is a property of the library, computable without seeing any query.

It matters because a spectrum is redistributed with its metadata, and an annotation that names one
compound where the spectrum cannot separate it from three others is metadata that misleads whoever
reuses it.

## What it is not

- Not library de-duplication. Nothing is deleted; the classes are an annotation-reporting aid.
- Not a global partition of the library. Classes are anchored per entry (see below).
- Not a replacement for MS-DIAL's per-query candidate list. MS-DIAL applies `FilterByThreshold` *before*
  capping at 3 hits per annotator, so a near-identical entry that narrowly fails the threshold vanishes
  from the candidate list. The library-side analysis is the only thing that catches that.
- Not the same relation as a learned-embedding neighbourhood, which contains structural analogs. That
  belongs to the separate proposed `SM` evidence tag.

## The computation

| Step | What | Why |
| --- | --- | --- |
| 1 | Normalize instrument strings to FT / TOF / IT / QQQ / UNKNOWN and collision energy to a value plus a bin | Free-text reality. In the lab's public positive release 24.9% of records carry an unusable `COLLISIONENERGY`, and values mix `20`, `20.0` and `45HCD`. A missing CE becomes `CE_UNKNOWN`, never 0.0 |
| 2 | Group by (InChIKey first block, ion mode, precursor type, instrument class, CE bin) and build one consensus spectrum per group | 501,264 public positive records carry only 32,264 unique first blocks. Without this the classes are dominated by replicate measurements of the *same* compound |
| 3 | Compare only within (ion mode, precursor type, precursor m/z within max(0.01 Da, 10 ppm)) | Anything outside that window is separable by MS1 mass alone. Measured 2,661-fold reduction in pair count |
| 4 | Admissibility gate before any similarity: both spectra need at least 6 peaks at or above 1% relative intensity, and the pair needs at least 4 matched peaks | 50% of records carry under 10 peaks. Two spectra with two peaks each look identical without that meaning anything. A pair failing the gate is `insufficient_evidence`, **a distinct outcome from indistinguishable**, and is never merged into a class |
| 5 | Require **both** the symmetric weighted cosine and the Li-Fiehn entropy similarity to pass | Reverse dot product is forbidden here: it is an asymmetric query-versus-reference quantity, and a class relation must be symmetric |
| 6 | Require formula agreement; otherwise tag `isobaric_not_isomeric` and refuse the edge | 6.4 million of 26.8 million cross-skeleton pairs in the window have *different* formulas. Those are artifacts of the window width, resolvable by formula/mass evidence alone |
| 7 | Anchored neighbourhood: `N(A) = {A} + {X : sim(A,X) >= threshold}`, plus a flag for whether `N(A)` is a clique | Defined relative to A, so there is no seed and the result is deterministic. Greedy clustering would be seed-dependent. Similarity is not transitive, so no transitive closure is ever taken |
| 8 | Report the product ions that would separate the members, and which evidence would break the tie | Turns "A or B" into "A or B, separable by m/z X", which is actionable |

### Symmetry is enforced, not assumed

MS-DIAL's weighted dot product is asymmetric in two places: `peakCountPenalty` counts peaks on the
reference side only, and its 0.01 intensity cutoff applies to the measured side only. Measured
consequence: `msdial_weighted_dot_product(three_peaks, one_peak)` returns 0.385 one way and 0.88 the
other. Asymmetry is fine for query-versus-library scoring, where one side really is the reference, but a
class asserts a mutual relation. `similarity.weighted_cosine` is therefore a symmetric variant with no
penalty and a cutoff that fires only where both sides are below it, and it divides by the product of the
two scalars rather than dividing twice, so `sim(A,B)` equals `sim(B,A)` bit-identically — a one-ulp
difference could otherwise flip a comparison sitting exactly on the threshold. Verified over 20,000
random pairs for all three measures.

### Anchoring produces one class row per member

A mutually indistinguishable group of N members yields N `ambiguity_class` rows, one per anchor, with
identical membership. That is forced by the non-transitivity requirement: `N(B)` can be `{A,B,C}` while
`N(A)` is only `{A,B}`, and that is only representable per anchor. Collapse rows by membership set when
counting groups for a human.

### The threshold is a convention, not a measurement

There is no ground-truth set of isomer pairs known to be separable by RT, CCS or an authentic standard,
so no error rate can be attached to the threshold. Every parameter therefore travels with the class in
`condition_scope_json` (`definition_rules`), so a class can be recalibrated and recomputed. Provisional
default: weighted cosine >= 0.90 and entropy similarity >= 0.85, at 0.025 Da.

A class computed on public spectra alone is **systematically incomplete** — the lab's public positive
release holds 501,264 records against 1,908,422 private ones — so `library_scope_json` records exactly
which libraries and versions were compared. A class id without a library scope is not reproducible.

## Measured behaviour

Public negative VS20 release, 44,353 records, 43 seconds end to end on one core (21 s ingest,
22 s classification), pure standard library.

| Quantity | Value |
| --- | --- |
| Records read / skipped | 44,353 / 0 |
| Consensus spectra | 32,403 |
| Comparison blocks | 4,005 |
| Pairs compared | 711,380 |
| Refused: condition mismatch | 621,361 |
| Refused: insufficient evidence | 52,644 |
| Refused: isobaric, not isomeric | 4,573 |
| Edges admitted | 5,521 |
| Anchored class rows | 5,216 |
| **Distinct membership groups** | **1,846** |
| Singletons | 39,137 |

Distinct groups by tier: 1,129 `same_skeleton`, **438 `different_skeleton`**, 272 `mixed`,
7 `unknown_skeleton`. No class survived with a formula disagreement, which is the formula guard working.
4,729 of 5,216 anchored classes are cliques; the remaining 487 are the real non-transitive cases that a
forced partition would have misrepresented. 4,984 of 5,216 carry discriminating ions.

Condition mismatch dominates the refusals at 87%, which is the honest consequence of requiring an
established instrument class and CE bin. 1,022 classes hold only under `CE_UNKNOWN` and the run warns
about them, because a class asserted under an unestablished condition is weaker than the rest.

### Real examples, all from the public release

| Min pairwise cosine | Members | Condition | Resolved by |
| --- | --- | --- | --- |
| 0.9999 | Diethyl phthalate / Monobutyl phthalate (C12H14O4) | FT, CE 30-40 | IM, RS |
| 0.9999 | N-Acetylserine / O-Acetylserine (C5H9NO4) | CE 30-40 | RT, RS |
| 0.9999 | Adenosine 5'-monophosphate / Adenosine 3'-monophosphate (C10H14N5O7P) | FT, CE 20-30 | RS only |
| 0.9998 | Diphenyl isophthalate / Phenolphthalein (C20H14O4) | FT, CE 60-70 | IM, RS |
| 0.9403 | Four distinct C34H44O19 glycosides, four distinct skeletons | unestablished | RS |

3'-AMP against 5'-AMP is the case worth dwelling on: a well-known real confusion in metabolomics that,
without this machinery, is reported under one of the two names with no indication that the spectrum did
not choose between them.

## Reporting a class

`ambiguity_class_for(database, reference_spectrum_id)` returns the class anchored on the entry a query
matched, with members ordered anchor-first and shaped to drop straight into
`annotation.CandidateInput`. Three tiers behave differently and must not be collapsed:

- **`same_skeleton`** — the 2D constitution is determined and only stereochemistry or isotopic labelling
  is open. Report one structure flagged stereo-unresolved; do not spend an ambiguity class on it.
- **`different_skeleton`, `same_formula`** — genuine constitutional isomers. The target case. Report
  `L3-SP` with the ordered candidate list, or step down to `L3-SC` / `L3-SI` when only a shared
  substructure survives.
- **`isobaric_not_isomeric`** — an artifact of the precursor window. Split it; formula/mass evidence
  already resolves it.

## Command line

```powershell
msdial-spectrum-catalog ingest-reference-library catalog.sqlite library.msp `
  --library-name "MSMS-Public_all-neg-VS20" --library-version VS20
msdial-spectrum-catalog compute-ambiguity catalog.sqlite
msdial-spectrum-catalog show-ambiguity catalog.sqlite "urn:msdial:reference-spectrum:..."
```

`--limit` and `--precursor-mz-range` bound the ingest for a pre-test. `--allow-condition-mismatch`
compares across instrument classes and CE bins and weakens every class it produces, which is why it is
off by default.

## Open

- No ground-truth separable-isomer set, so the threshold carries no error rate.
- Classes computed on the public releases only; including the private libraries will change them.
- The GC-EI libraries have no `PRECURSORMZ`, `PRECURSORTYPE` or `IONMODE`, so the whole blocking scheme
  needs replacing for them — block on exact mass and retention index instead.
- `SPLASH` is present on every record of the public negative release but every value is unique, so it
  does not work as a cheap exact-duplicate pre-filter here.
