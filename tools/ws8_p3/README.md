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
| `terms.tsv` | One row per formula token, polarity (`both` records) and adduct, with mass, role, flags, group and support ids |
| `groups.json` | Grouped term rows by (family, kind, ion mode, adduct, formula), with record ids, citations, paper counts and the roll-up counts of their term |
| `term_support.tsv` | The roll-up over adducts: one row per (family, kind, ion mode, formula) with the union of papers and standards-based papers |
| `counts_family_kind.tsv` | Records per unified family and kind, per source, with a total row |
| `counts_family_group.tsv` | Records, term rows (all and grouped), groups, terms, papers and standards-based papers per family group and family |
| `manual_queue.tsv` | Records to read by hand: stage-3 `other`, and any kind not in the map |
| `mass_flags.tsv` | One row per flagged formula token: `printed_mismatch`, `lowres_printed`, and the advisory `note_mismatch`, `note_lowres` |
| `unmapped.tsv` | Field, value and record ids of everything that could not be mapped or was not made a term |
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
- `family-override-p3.tsv`: record ids whose family label is broader than their class scope. The stage-3 records
  labelled `flavonoid_C_glycoside` whose scope is di-C-glycosides only go to `flavonoid_di_c_glycoside`, the
  family of the `flavone_CC_glycoside` label. `family_original` is kept and `family_override` gives the reason;
  a row whose `family_original` no longer matches the record stops the run.
- `instrument-map-p3.tsv`: ordered regexes to `instrument_class` (`qtof`, `orbitrap_ft`, `ion_trap`,
  `qqq_or_qtrap`, `mixed`, `compilation`, `unknown`) and `instrument_strat`. Ion traps are their own class.
  **A mixed record counts toward the instrument that produced the fragmentation** (`instrument_strat`): an ion
  trap used for MSn with a QTOF or Orbitrap only for accurate mass counts as `ion_trap`. When both instruments
  fragmented, or the record does not say which one produced its ions, `instrument_strat` stays `mixed` and the
  record belongs to no single-instrument stratum.
  Orbitrap and FT-ICR hybrids with an ion trap (LTQ Orbitrap, LTQ FT, LIT-FTICR, Orbitrap Elite and Velos, and
  the tribrids Fusion, Lumos, Eclipse, Ascend) are `mixed`, and their stratum comes from the activation the
  record names (`instrument_activation`), read from `collision`, from explicit phrases in the instrument text
  ('CID in the linear ion trap'), and from `condition` only when those say nothing: trap CID, LIT or CAD gives
  `ion_trap`; HCD or beam-type CID gives `orbitrap_ft`; IRMPD or SORI-CID in the FT cell gives `orbitrap_ft`;
  both kinds give `mixed`; nothing gives `mixed`. In-source and pre-trap activation and negated mentions
  ('not CID') are ignored. A bare normalised or percentage energy ('NCE 35%') means trap CID on the LTQ-based
  hybrids, where it is the default, and stays `mixed` on the tribrids.
- `paper-map-p3.tsv`: DOIs of records without a `paper_id` (the NPR review and PlantMAT).

## Rules

- Ion mode: `positive`, `negative`, `both` or `unknown`; a leading "not stated" stays `unknown` even when the
  extractor added an inference, and `ion_mode_inferred` marks a stated polarity that was inferred.
- Adduct: a list of canonical strings (`[M+H]+`, `[M-H]-`, `[M]+`, `[M+HCOO]-`, `[M-H]-` and `[M-2H]2-` as two
  entries, ...); unicode minus and dashes become `-`, `M+` becomes `[M]+`, a missing or "not stated" adduct is
  `unstated`. Parenthetical remarks are ignored, text after an arrow is a product, and in "X from Y" the
  precursor is Y.
- Term rows: one per formula token, polarity and adduct. A `both` record gives one row per polarity, each with
  the adducts of that sign (a charged token keeps only the polarity of its own charge); a record with several
  adducts of one polarity gives one row per adduct (`expanded_from_adducts`).
- Terms: tokens of the `formula` field (separated by `;` or `,`, with an optional parenthetical label, a charge
  `+`, `-`, ` 2-` or superscript, and a radical mark: `•`, or `*` or `.` before the charge, or `*` after it).
  Tokens inside ion notations such as `[M-H-CO2]-` are not terms. When `formula` is empty, tokens are taken from
  `value` only if they start with carbon and contain hydrogen (`formula_from = value`); otherwise the record has
  no term (`no_formula`). A charged token is an ion, an uncharged one a neutral; `kind_conflict` marks a charged
  neutral loss or an uncharged product ion, and `odd_electron_neutral` an odd-electron neutral loss (a radical
  loss such as CH3•, or, above about 100 Da, usually an ion written without its charge). The canonical formula
  is Hill order, then the charge (`[C21H18O12]2-` for multiply charged ions).
- Formula fragments that are not terms go to `unmapped.tsv`: `formula_fragment` for multiplied formulas (`2H2O`),
  R-group templates (`RCOOH`) and symbols that are not elements; `formula_composite` for expressions that join
  formulas into one loss (`CH2O+2H2O`, `RCOOH + NH3`). The record gets `formula_flags`.
- Roles and grouping: a term row is kept out of groups (`group_exclusion`, `group_id` empty) when its token is
  the precursor or intact molecule (`role = precursor`: its label or the words just before it say precursor,
  parent, molecular ion, `[M-H]-`, `M+`, `MH+` ...), when it belongs to an `absence` or `veto` record
  (`role_unverified`: the token may be the absent ion, the observed one, an isobar or a competitor, and a reader
  must say which), or when it occurs only inside a composite expression (`composite_part`). These rows stay in
  `terms.tsv` with the same flags.
- Masses: element masses equal `tools/ws8_p2/p2_diagnostics.py:76` (RDKit most-common isotope) and the electron
  mass 0.000549 Da its line 77. m/z = (sum of masses - z * electron) / |z|. C, H, N, O, P, S, Na, K and Cl only;
  any other element gives `unsupported_element` and no mass.
- Printed values: decimals within 0.5 Da of the computed mass are compared with the product tolerance
  max(0.005 Da, 10 ppm), plus half the last printed digit. A number written next to a different formula
  ('C30H46NaO4+ 493.3288') is that formula's mass and is not compared. Numbers in `value` are used first; `note`,
  which often holds recomputed masses, is used only when `value` has none, so a misprint is not hidden by its
  correction, and a disagreement found there is advisory (`note_mismatch`, `note_lowres`). With no match, a value
  printed with fewer than 3 decimals whose half last digit is not below the tolerance gives `lowres_printed`; any
  other gives `printed_mismatch`.
- Groups keep every record id and citation (DOIs as `doi:10.xxx` in lower case; NPR `primary_refs` as secondary
  citations) and count `n_papers` (distinct `paper_id`, or DOI when there is none) and `n_standards_papers`
  (papers with `standards_based` true). The same term is split over adduct groups (`unstated` included), so
  `term_support.tsv` and the `support_*` fields of each group count the papers of the term over every adduct.
  These counts are evidence, not votes.

## Known limitations

- Family labels map per label. Records under `alkaloid_general` that name one subclass stay `generic` until a
  reader relabels them. Records labelled as C-glycosides whose scope is di-C only are moved by the id list in
  `family-override-p3.tsv`; records whose scope covers mono- and di-C-glycosides together keep the mono label.
- `acylcarnitine` (and its confounder `betaine`) is `in_scope` because the directive names it among the P3 v1
  families, although it is a fatty-acyl class (CAR in the MS-DIAL lipid model); whether "lipids out of scope"
  covers it is the author's decision, and a one-line change in `family-map-p3.tsv`.
- The precursor role reads only the words next to a token; a precursor named elsewhere in the record, or an
  uncharged product ion in a neutral-loss record, is not detected (the second is flagged `odd_electron_neutral`
  when its formula is odd-electron).

Tests use synthetic records only: `tests/test_ws8_p3_normalise.py`.
