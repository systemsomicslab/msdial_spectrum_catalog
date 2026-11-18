# WS8 P3: candidate normaliser

`normalise_candidates.py` turns the P3 literature rule-candidate records (NPR/PlantMAT, stage 2, stage 3) into
one normalised table for drafting diagnostic-ion rules. It maps the two `kind` vocabularies and the family
labels onto one vocabulary, normalises ion mode, adduct and instrument, splits each `formula` into term rows,
recomputes every mass from its formula, flags printed values that disagree, and groups term rows across papers.
It is stdlib only (Python 3.10). Its outputs are local and are never committed.

```powershell
python -B tools/ws8_p3/normalise_candidates.py `
  --inputs npr=<literature>/rule_candidates_npr_plantmat.json `
           stage2=<literature>/review/stage2/rule_candidates_oa.json `
           stage3=<literature>/review/stage3/rule_candidates_manual.json `
  --out <annotation-root>/ws8-p3/candidates-<yyyymmdd>
```

Input labels must appear in `kind-map-p3.tsv`. The output folder must sit in a folder named `ws8-p3`, be named
`candidates-<yyyymmdd>`, and be new or empty. Nothing is dropped: a value that cannot be mapped goes to
`unmapped.tsv`, and the run fails if `records.tsv` does not hold every input record exactly once (checked in
memory and again after re-reading the file).

## Outputs

| File | Content |
|---|---|
| `records.tsv` | One row per input record: every normalised column next to the original field |
| `terms.tsv` | One row per formula token (per polarity for `both` records), with mass, flags and group id |
| `groups.json` | Term rows grouped by (family, kind, ion mode, adduct, formula), with record ids, citations, paper counts |
| `counts_family_kind.tsv` | Records per unified family and kind, per source, with a total row |
| `counts_family_group.tsv` | Records, term rows, groups, papers and standards-based papers per family group and family |
| `manual_queue.tsv` | Records to read by hand: stage-3 `other`, and any kind not in the map |
| `mass_flags.tsv` | Terms flagged `printed_mismatch` or `lowres_printed` |
| `unmapped.tsv` | Field, value and record ids of everything that could not be mapped |
| `run.json` | Input paths and sha256, code commit, script and mapping sha256, parameters, counts, reconciliation |

## Mapping files (labels only)

- `kind-map-p3.tsv`: (source, kind_original) to kind: `product_ion`, `neutral_loss`, `ion_ratio`, `sequence`,
  `precondition`, `adduct_behaviour`, `absence`, `veto`, `manual`. Stage-3 `diagnostic_absence` is `absence`
  (an ion stated to be absent). Stage-2 `veto` stays `veto`: the stage-3 README says the two correspond only
  roughly, and stage-2 vetoes also hold confounders, isobars and stated limits of MS/MS, which stage 3 put under
  `other`. Both have `evidence_sense = negative`. Stage-3 `other` is `manual` and is never guessed.
- `family-map-p3.tsv`: family_original to a lower-case family, a `family_group` for drafting and a scope
  (`in_scope`, `out_of_scope` for lipids and acyl-lipid classes, `generic` for tool notes, heuristics and
  substructure priors). Only spelling, case and synonyms are merged; distinct classes keep their label.
- `instrument-map-p3.tsv`: ordered regexes to `instrument_class` (`qtof`, `orbitrap_ft`, `ion_trap`,
  `qqq_or_qtrap`, `mixed`, `compilation`, `unknown`) and `instrument_strat`. Ion traps are their own class.
  **A mixed record counts toward the instrument that produced the fragmentation** (`instrument_strat`): an ion
  trap used for MSn with a QTOF or Orbitrap only for accurate mass counts as `ion_trap`. When both instruments
  fragmented, or the record does not say which one produced its ions, `instrument_strat` stays `mixed` and the
  record belongs to no single-instrument stratum.
- `paper-map-p3.tsv`: DOIs of records without a `paper_id` (the NPR review and PlantMAT).

## Rules

- Ion mode: `positive`, `negative`, `both` or `unknown`; a leading "not stated" stays `unknown` even when the
  extractor added an inference, and `ion_mode_inferred` marks a stated polarity that was inferred.
- Adduct: a list of canonical strings (`[M+H]+`, `[M-H]-`, `[M]+`, `[M+HCOO]-`, `[M-H]-` and `[M-2H]2-` as two
  entries, ...); unicode minus and dashes become `-`, `M+` becomes `[M]+`, a missing or "not stated" adduct is
  `unstated`. Parenthetical remarks are ignored, text after an arrow is a product, and in "X from Y" the
  precursor is Y. `both` records give one term row per polarity, each with the adducts of that sign; a charged
  token keeps only the polarity of its own charge.
- Terms: tokens of the `formula` field (separated by `;` or `,`, with an optional parenthetical label, a charge
  `+`, `-`, ` 2-` or superscript, and a radical dot). Tokens inside ion notations such as `[M-H-CO2]-` are not
  terms. When `formula` is empty, tokens are taken from `value` only if they start with carbon and contain
  hydrogen (`formula_from = value`); otherwise the record has no term (`no_formula`). A charged token is an ion,
  an uncharged one a neutral; `kind_conflict` marks a charged neutral loss or an uncharged product ion. The
  canonical formula is Hill order, then the charge (`[C21H18O12]2-` for multiply charged ions).
- Masses: element masses equal `tools/ws8_p2/p2_diagnostics.py:76` (RDKit most-common isotope) and the electron
  mass 0.000549 Da its line 77. m/z = (sum of masses - z * electron) / |z|. C, H, N, O, P, S, Na, K and Cl only;
  any other element gives `unsupported_element` and no mass.
- Printed values: decimals within 0.5 Da of the computed mass are compared with the product tolerance
  max(0.005 Da, 10 ppm), plus half the last printed digit. Numbers in `value` are used first; `note`, which
  often holds recomputed masses, is used only when `value` has none, so a misprint is not hidden by its
  correction. With no match, a value printed with fewer than 3 decimals whose half last digit is not below the
  tolerance gives `lowres_printed`; any other gives `printed_mismatch`.
- Groups keep every record id and citation (DOIs as `doi:10.xxx` in lower case; NPR `primary_refs` as secondary
  citations) and count `n_papers` (distinct `paper_id`, or DOI when there is none) and `n_standards_papers`
  (papers with `standards_based` true). These counts are evidence, not votes.

Tests use synthetic records only: `tests/test_ws8_p3_normalise.py`.
