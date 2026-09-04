# Annotation Evidence Model

The model follows the Shin-MassBank annotation metadata proposal and keeps three dimensions separate:

1. Annotation level: `L1` to `L5`.
2. Claim specificity for Level 3.
3. Evidence inventory: which kinds of evidence support the claim.

A claim is never promoted because it has many evidence tags, and the level is never derived from the tag
set. The published use cases assign the same bracket `[SL,FM]` to `L1`, `L2` or `L3-SP` depending on
curator judgement, so an automatic level assigner would disagree with the proposal's own examples. The
level is always an explicit input.

## Vocabulary versions

The proposal exists in more than one generation and they disagree. The vocabulary is therefore data, not
code: one JSON file per version under `src/msdial_spectrum_catalog/vocabulary/registry/`, resolved by
`load_vocabulary(version)`.

| Version | Status | Claim tokens | Source |
| --- | --- | --- | --- |
| `smb-v2-consensus` | accepted, **default** | `SP` `SC` `SI` `CP` | V2 document sections 2 and 4, agreed internationally |
| `smb-v1-draft` | superseded | `SP` `CP` `MO` `CL` | earlier document sections 3/6/8, Table S9, and the slide deck |
| `smb-v2.1-proposal` | draft | same as the consensus | a strict extension adding one proposed evidence tag, `SM` |

`smb-v2-consensus` is the default because the V2 document's sections 1 to 4 are the agreed consensus.
Sections 3, 5, 6 and 8 of that same document were not updated and still print the superseded tokens, which
is why no notation in its use-case table resolves under a single version: the claims already read `SC`/`SI`
while the brackets still read `CX`. Those strings are recorded verbatim in `registry/use_cases.json` with
`parseable: false` rather than silently reconciled.

### The token `CP` is re-bound, not aliased

    smb-v1-draft      CP  =  "component proposed"  =  smb:claim/substructure_complete
    smb-v2-consensus  CP  =  "class proposed"      =  smb:claim/class

This is the single most important property of the vocabulary. A flat alias table cannot express it, and a
column storing the bare token is permanently ambiguous. Therefore:

- stored rows key on `(vocab_version, concept_id)`; the two-letter token is produced only at emit time;
- `parse_notation` requires an explicit `vocab_version` with no default;
- `parse_notation_any` returns one reading per version that resolves, so a legacy string can be inspected
  rather than guessed at.

Concept identifiers are fixed once and never renamed, because they are the join key across versions.

```
$ msdial-spectrum-catalog notation 'L3-CP[FM,DF]'
  smb-v1-draft      -> Component proposed
  smb-v2-consensus  -> Class proposed
```

## Level 3 claim tags (`smb-v2-consensus`)

| Tag | Term | Claim |
| --- | --- | --- |
| `L3-SP` | Structure(s) proposed | one 2D constitutional structure, not confirmed by a reference standard; structures may be proposed as a list |
| `L3-SC` | Substructure complete | components, building blocks or modifications are defined, but their positions or arrangements are unresolved |
| `L3-SI` | Substructure incomplete | a diagnostic substructure or motif is proposed, without a component-level representation |
| `L3-CP` | Class proposed | a chemical or metabolite class only |

Specificity, not confidence: `SP` is more specific than `SC`, which is more specific than `SI`, which is
more specific than `CP`. Strength of support is carried entirely by the evidence tags and the criteria set.
`CP` may be appended to `SC` or `SI`; `SP+CP` is exceptional and `SC+SI` is avoided. `validate_combination`
warns rather than raising, because the combination rules are guidance and the workshop is still running.

## Evidence inventory (`smb-v2-consensus`, 12 tags)

`RS` reference standard, `FM` formula/mass, `SL` spectral library, `DF` diagnostic fragment,
`RT` retention, `IM` ion mobility, `IS` in silico, `MN` molecular network, `HO` homologue,
`CO` contextual, `UN` unclassified additional, `OS` other spectroscopy.

Emitted in that order, which is the order printed in the source table. Contextual evidence is `CO`, not the
earlier draft's `CX`.

`OS` and `UN` overlap in the printed text: `UN`'s definition also names UV and IR, and section 6 use case 12
tags UV evidence as `UN`. Ruled on 2026-09-03: **UV, IR and ECD evidence is `OS`**, and the use-case map is
outdated rather than the section 4 definitions. `UN` stays the fallback for evidence that fits no other tag,
such as 2D chromatography.

### `SM`, proposed: similarity without peak correspondence

`SM` was dropped when spectral similarity was folded into `SL`, which leaves representation-based
similarity — a DreaMS embedding cosine, for instance — with no home among the twelve consensus tags.
Folding it into `IS` instead would stretch that tag to cover both computation-from-structure and
comparison-of-two-measurements. `smb-v2.1-proposal` therefore adds `SM` with status `proposed` and the
concept `smb:evidence/model_spectral_similarity`, distinguished from its neighbours by what is compared and
whether any fragment can be pointed at:

| Tag | Compares | Claims | Peak correspondence |
| --- | --- | --- | --- |
| `SL` | query against a reference spectrum of a known compound | the reference's identity | yes, peak by peak |
| `MN` | query against another observed spectrum | a relationship: an edge plus an interpretable mass shift | yes, shared peaks |
| `IS` | a structure, computationally | a ranking from prediction | via the prediction |
| `SM` | two spectra in a learned representation | resemblance only | **none** |

That absence is the point: no fragment can be named as the reason for the similarity, which is why the
evidence cannot be reported as `SL` or `DF`. Two constraints come with the tag. A subtype qualifier is
mandatory — `SM:dreams-emb-cos`, say — because a learned-embedding 0.8 is not comparable to a peak-matching
cosine 0.8 and must never be read against a cosine threshold; and the record must carry the model identity,
since embeddings from different checkpoints are not interchangeable. Contrastive spectrum models are
typically optimized for analog retrieval, so a high value supports a class or substructure claim far more
often than a proposed structure: `SM`-only evidence should not carry an `L3-SP` claim.

Because `SM` is `proposed`, strict mode rejects it and only permissive mode resolves it. Migrating a reading
back to `smb-v2-consensus` drops the tag and marks the reading unresolved, rather than silently recoding it
as `SL` or `IS`.

## Notation

    notation := level [ '-' claim ] [ '[' evidence ']' ]
    level    := 'L1' | 'L2' | 'L3' | 'L4' | 'L5'
    claim    := claim_token { '+' claim_token }
    evidence := [ ev_token { ',' ev_token } ]

A claim token is legal only at `L3`. The bracket is always emitted, so `L5` round-trips as `L5[]`.
Alternation is not part of the grammar: `L3-SC/SI[...]`, `L3-SP[...] or L3-SC[...]` and prose
conditionals are rejected rather than guessed at, because an alternation is two claims and belongs in
structured fields. A candidate set is expressed by `annotation_candidate` rows, not by a string.

`|` is never a delimiter: mzTab-M treats it as the multi-value separator inside a cell.

## Record shape

`annotation_assertion` is the claim, anchored on one spectrum. `annotation_claim_component` makes a
combination queryable instead of hiding it in a `+`-joined string. `annotation_candidate` is the ordered
candidate list that makes an honest "A or B" representable, with `UNIQUE(assertion_id, rank)` and
contiguous ranks enforced at write time. `annotation_evidence` carries, per tag, the verbatim token, the
resolved concept, the measured value, the threshold it was compared against, and whether it passed, so
evidence can be re-tagged when the criteria change without re-running annotation. `criteria_set` plus
`criteria_rule` hold the study-level thresholds, versioned and stamped with the vocabulary version they
were authored against. `annotation_tool_run` records which tool, version and parameters produced a
candidate.

Scores always travel with their convention, because MS-DIAL squares internally and un-squares on the way
out. The chain, verified in source:

| Layer | Value | Where |
| --- | --- | --- |
| `GetSimpleDotProduct` return | cos-squared of the sqrt-intensity vectors | `Math.Pow(covariance,2)/scalarM/scalarR` |
| `GetWeightedDotProduct` return | cos-squared of sqrt(intensity x m/z) vectors, times a peak-count penalty of 0.75 / 0.88 / 0.94 / 0.97 for 1 to 4 matched peaks | same function; its own doc comment says "the square of a typical dot product" |
| `MsScanMatchResult.Squared*` | stores those values as-is | `MsScanMatching.cs:657`, `MassAnnotator.cs:103` |
| threshold comparison | `Squared* >= Squared*CutOff`, i.e. 0.36 and 0.64 | `MsScanMatching.cs:442`, `MassAnnotator.cs:192` — no exporter reads these fields |
| parameter file and UI cutoffs | non-squared, 0.6 and 0.8; `WeightedDotProductCutOff` is a `sqrt` getter | `MsRefSearchParameterBase.cs:33-45` |
| **`.mdpeak` and `.mdalign` columns** | **plain cosine** | `IAnalysisMetadataAccessor.cs:123-125` and `IMetadataAccessor.cs:114-116` write the non-squared computed properties |

So everything user-facing is already consistent in cosine terms and the squaring never leaks; ingested
rows are therefore recorded as `score_convention = 'cosine'`. Two consequences matter for anything that
calls the scoring kernel directly rather than reading a text export: `GetWeightedDotProduct` and friends
return cos-squared, so a 0.9 cosine threshold must be compared against 0.81 or the value must be
square-rooted first; and `GetBatchSimpleDotProduct` returns a dense matrix pre-filled with `-1` for pairs
it did not compute (`IsAvailableSpectrum` failed), where `-1` is a not-computed sentinel and not a
similarity of minus one.

## What MS-DIAL already asserts, and what it does not

`msdial_annotation_result` holds MS-DIAL's own per-feature annotation block, which the ingester previously
discarded entirely. It is not a Level-3 claim; it is the raw material for one. MS-DIAL distinguishes four
outcomes, and collapsing them would overstate the evidence:

| `annotation_kind` | MS-DIAL name | Meaning | May support |
| --- | --- | --- | --- |
| `msms_matched` | no prefix | reference match with a product-ion spectrum | `SL`, and `DF` from matched peaks |
| `low_score` | `low score: ` | product-ion spectrum acquired, search criteria failed | `SL` with `passed = 0`, as a broad candidate |
| `precursor_only` | `no MS2: ` | no product-ion spectrum was acquired at all | `FM` only, never `SL` |
| (no row) | `Unknown` | no candidate | nothing |

`MsdialCore/Utility/DataAccess.cs` `SetMoleculeMsPropertyAsSuggested` writes `"no MS2: "` when
`MS2RawSpectrumID < 0` and `"low score: "` otherwise. `"w/o MS2: "` is the MS-DIAL 4 spelling and is
commented out in MS-DIAL 5, but MS-DIAL 4 exports are still ingestable so both are recognized.

Two MS-DIAL export artefacts are normalized on the way in, because leaving them would read as
measurements. Both trace to one mechanism: `MsScanMatching`'s scoring functions return **-1** when there
is nothing to compare, that -1 is stored in `MsScanMatchResult.Squared*`, and the non-squared getters
clamp it with `Math.Max(value, 0f)`.

- `Matched peaks count` and `Matched peaks percentage` have no clamping getter, so the raw `-1` reaches
  the column. A count cannot be negative, so `-1` becomes `NULL`.
- The dot-product getters do clamp, so an older `.mdpeak` carries `0.000` on the same row -- a clamped
  sentinel, neither a measurement nor an unset field. Left as `0.0` it reads as "compared, scored zero",
  a stronger statement than "never compared", so an exact `0.0` on a `precursor_only` row becomes
  `NULL`. Only `precursor_only` rows are normalized: a `low_score` row did have a spectrum compared, so
  a zero there is a real result and survives.

Both were fixed upstream in MsdialWorkbench PR #785, which routes every score column of both exports
through `AnnotationScoreFormat.Score` and writes `null` when
`MsScanMatchResult.IsSpectrumComparisonPerformed` is false. New outputs carry `null` directly, so the
normalization above is a no-op on them; it stays necessary for outputs produced before that change.

One thing cannot be recovered from an older `.mdalign`: its exporter bound a float overload that printed
`null` for any value within 1e-10 of zero, so a *genuine* cosine of zero -- a real comparison that found
no overlapping fragment -- was written as `null` too and is indistinguishable from a comparison that
never happened. The related name convention was fixed in PR #784, which adds
`CommonStandard/Utility/AnnotationName.cs` and uses it in place of the stale `w/o MS2` literal.

`candidate_name` holds the name with the status prefix stripped, and `candidate_is_named` is `0` when the
reference record's name is an in-house identifier rather than a compound name. Only a named
`msms_matched` row is eligible to support `SL` evidence. On the reference LC-MS negative demo that is
1,639 of 7,510 annotation rows.

`Total score` is MS-DIAL's unnormalized weighted composite, not a cosine — it reaches 2.611 on that demo —
so `score_convention` describes the three dot-product columns only. Note also that `TotalScoreCutoff`
exists in the parameter object but is not applied by `MsScanMatchResultEvaluator`, so it must not be
published as an applied threshold.

## Known gaps

- Candidate lists cannot yet be recovered from MS-DIAL's text exports. `StandardAnnotationProcess` keeps
  `NUMBER_OF_ANNOTATION_RESULTS = 3` hits per annotator and persists them into `.pai`/`.arf`, but every
  text exporter reads only `MatchResults.Representative`, and mzTab-M hardcodes `rank = "1"` and
  `SME_ID_REF_ambiguity_code = "null"`. Interoperable "A or B" reporting needs a Console-side
  per-candidate exporter.
- Ambiguity classes are computed and persisted; see `docs/ambiguity_classes.md`. What is still missing
  is a ground-truth set of isomer pairs known to be separable by RT, CCS or an authentic standard, so
  the similarity threshold remains a convention with no error rate attached.
- In-silico predicted spectra have no producer yet, although `reference_library.library_kind` already
  distinguishes them so one can never be reported as spectral-library evidence.
- Concept identifiers are stored as free text with no vocabulary table to constrain them; `vocab_version`
  is the only provenance.
