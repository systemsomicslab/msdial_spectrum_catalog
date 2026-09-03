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
| `smb-v2-consensus` | active, default | `SP` `SC` `SI` `CP` | V2 document sections 2 and 4, agreed internationally |
| `smb-v1-draft` | superseded | `SP` `CP` `MO` `CL` | earlier document sections 3/6/8, Table S9, and the slide deck |

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

Emitted in that order, which is the order printed in the source table. Two changes from the earlier draft
are worth noting: contextual evidence is `CO`, not `CX`; and `SM` (spectral similarity) is absent from the
consensus set, so a learned-embedding similarity such as DreaMS has no settled home among these twelve.
`OS` and `UN` overlap in the agreed text — `UN`'s definition also names UV and IR, and use case 12 tags UV
evidence as `UN` — so curators need an explicit rule for choosing between them.

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

Scores always travel with their convention. MS-DIAL reports dot products as **squared** cosines — its C#
cutoffs are literally `.6F * .6F` and `.8F * .8F`, and `WeightedDotProductCutOff` is a computed property
returning `sqrt` of the stored value. A cosine of 0.9 is 0.81 in MS-DIAL's own numbers, which is precisely
the range where a merge decision lives, so `score_convention` is stored explicitly and never inferred.

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
measurements:

- `Matched peaks count` and `Matched peaks percentage` of `-1` mean "not applicable". A count cannot be
  negative, so `-1` becomes `NULL`.
- For a `precursor_only` row, `.mdalign` writes `null` dot products while `.mdpeak` writes `0.000` — the
  unset default of `MsScanMatchResult`'s float fields. Left as `0.0` it reads as "compared, scored zero",
  which is a stronger statement than "never compared", so an exact `0.0` becomes `NULL`. A non-zero value
  is kept so a future MS-DIAL change surfaces as an anomaly instead of being discarded.

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
- Reference and in-silico predicted spectra are first-class (`reference_library`, `reference_spectrum`)
  but nothing populates them yet, so `spectrum_similarity` and `ambiguity_class` have no producer.
- Concept identifiers are stored as free text with no vocabulary table to constrain them; `vocab_version`
  is the only provenance.
