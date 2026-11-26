# WS8 P3: candidate normaliser

`normalise_candidates.py` turns the P3 literature rule-candidate records (NPR/PlantMAT, stage 2, stage 3) into
one normalised table for drafting diagnostic-ion rules. It maps the two `kind` vocabularies and the family
labels onto one vocabulary, normalises ion mode, adduct, instrument and collision regime, splits each `formula`
into term rows, recomputes every mass from its formula, flags printed values that disagree, and groups term rows
across papers. It is stdlib only (Python 3.10). Its outputs are local and are never committed.

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
| `records.tsv` | One row per input record: every normalised column next to the original field, with its record flags and the citation label (`citation`, `citation_original`, `citation_fix`) |
| `terms.tsv` | One row per formula token, polarity (`both` records) and adduct, with mass, role, flags, group and support ids |
| `groups.json` | Grouped term rows by (family, kind, ion mode, adduct, formula), with record ids, citations, paper counts, instrument strata, collision regimes, `kind_consistent` and the roll-up counts of their term |
| `term_support.tsv` | The roll-up over adducts: one row per (family, kind, ion mode, formula) with the union of papers and standards-based papers |
| `counts_family_kind.tsv` | Records per unified family and kind, per source, with a total row |
| `counts_family_group.tsv` | Records, term rows (all and grouped), groups, terms, papers and standards-based papers per family group and family |
| `manual_queue.tsv` | Records to read by hand: stage-3 `other`, and any kind not in the map |
| `mass_flags.tsv` | One row per flagged formula token: `printed_mismatch`, `lowres_printed`, and the advisory `note_mismatch`, `note_lowres` |
| `unmapped.tsv` | Field, value and record ids of everything that could not be mapped or was not made a grouped term |
| `run.json` | Input paths and sha256, code commit, script and mapping sha256, parameters, counts, override, paper-map and record-flag use, reconciliation |

## Mapping files (labels, record ids and formulas only)

- `kind-map-p3.tsv`: (source, kind_original) to kind: `product_ion`, `neutral_loss`, `ion_ratio`, `sequence`,
  `precondition`, `adduct_behaviour`, `absence`, `veto`, `manual`. Stage-3 `diagnostic_absence` is `absence`
  (an ion stated to be absent). Stage-2 `veto` stays `veto`: the stage-3 README says the two correspond only
  roughly, and stage-2 vetoes also hold confounders, isobars and stated limits of MS/MS, which stage 3 put under
  `other`. Both have `evidence_sense = negative`. Stage-3 `other` is `manual` and is never guessed.
- `family-map-p3.tsv`: family_original to a lower-case family, a `family_group` for drafting and a scope
  (`in_scope`, `out_of_scope` for lipids and acyl-lipid classes, `generic` for tool notes, heuristics and
  substructure priors). Only spelling, case and synonyms are merged; distinct classes keep their label.
- `family-override-p3.tsv`: record ids whose family label is broader than their class scope. Stage-3 records
  labelled `flavonoid_C_glycoside` whose scope is di-C-glycosides only go to `flavonoid_di_c_glycoside`; a
  `steryl_glycoside` record whose scope is steryl esters only goes to `steryl_ester` (out of scope, an acyl-lipid
  class), while a record covering steryl glycosides and esters together keeps `steryl_glycoside`.
  `family_original` is kept and `family_override` gives the reason; a row whose `family_original` no longer
  matches the record stops the run.
- `instrument-map-p3.tsv`: ordered regexes to `instrument_class` (`qtof`, `orbitrap_ft`, `ion_trap`,
  `qqq_or_qtrap`, `mixed`, `compilation`, `unknown`) and `instrument_strat`, and the activation patterns.
  Ion traps are their own class.
  **A mixed record counts toward the instrument that produced the fragmentation** (`instrument_strat`): an ion
  trap used for MSn with a QTOF or Orbitrap only for accurate mass counts as `ion_trap`. When both instruments
  fragmented, or the record does not say which one produced its ions, `instrument_strat` stays `mixed` and the
  record belongs to no single-instrument stratum.
  Orbitrap and FT-ICR hybrids with an ion trap (LTQ Orbitrap, LTQ FT, LIT-FTICR, Orbitrap Elite and Velos, and
  the tribrids Fusion, Lumos, Eclipse, Ascend, also when named without "Orbitrap") are `mixed`, and their stratum
  comes from the activation the record names (`instrument_activation`), read from `collision`, from explicit
  phrases in the instrument text ('CID in the linear ion trap'), and from `condition` only when those say
  nothing: trap CID, LIT or a bare CID or CAD gives `ion_trap`; HCD, beam-type CID, or CID with N2 as collision
  gas (the linear trap uses helium) gives `orbitrap_ft`; IRMPD or SORI-CID in the FT cell gives `orbitrap_ft`;
  both kinds give `mixed`; nothing gives `mixed`. In-source and pre-trap activation and negated mentions
  ('not CID') are ignored. A bare normalised or percentage energy ('NCE 35%') means trap CID on the LTQ-based
  hybrids, where it is the default, and stays `mixed` on the tribrids.
- `stratum-override-p3.tsv`: (source, id regex) to `instrument_strat` and `collision_regime` for the records of
  one paper, where one record names the activation and the others share its method without restating it. The
  computed stratum stays in `instrument_strat_computed` and the reason in `instrument_override`. Rows are
  inferences from the records until a reader confirms them against the paper.
- `term-override-p3.tsv`: (source, id, canonical formula) to a token role (`precursor`, `not_observed`,
  `printed_as`) where the record says it but the tokeniser cannot read it. A row naming a formula that the
  record does not hold stops the run.
- `paper-map-p3.tsv`: (source, id regex) or (source, DOI) to a DOI and a journal reference; DOIs and journal
  references only, never author names. It gives the DOI of records without a `paper_id` (the NPR review and
  PlantMAT), and, with a `reference`, replaces a citation label that names a wrong author, no author or an inverted
  name by `<reference>, <doi>` (`citation`; the input stays in `citation_original`, the row's note in
  `citation_fix`). A row without an id pattern applies to the records of its source whose own citation carries its
  DOI; the first matching row wins, and a matching record that cites another DOI stops the run. Papers, groups and
  counts do not change. The six label fixes are provisional WS0 fixes of 2026-11-19 (Stop A handback, section 6.1).
- `record-flags-p3.tsv`: (source, id regex, optional canonical formula) to a flag, `msn_level_misplaced` (an ion or
  loss of a later MSn stage given as MS2) or `source_misread` (the record misreads its source), with a short note.
  A flagged record keeps its term rows in `terms.tsv`, but they stay out of groups (`record_flag`) until a reader
  corrects the record; with a formula, only that token's rows. `records.tsv` lists the flags and notes; a row whose
  formula is not a token of a record it matches stops the run, and rows that match nothing are listed in `run.json`.
  The three rows are provisional WS0 fixes of 2026-11-19.

## Rules

- Ion mode: `positive`, `negative`, `both` or `unknown`; a leading "not stated" stays `unknown` even when the
  extractor added an inference, and `ion_mode_inferred` marks a stated polarity that was inferred.
- Adduct: a list of canonical strings (`[M+H]+`, `[M-H]-`, `[M]+`, `[M+HCOO]-`, `[M-H]-` and `[M-2H]2-` as two
  entries, ...); unicode minus and dashes become `-`, `M+` becomes `[M]+`, a radical mark after the charge
  (`•`, or `.` ending the field or a phrase: `[M+H]+.`) becomes `•`, a missing or "not stated" adduct is
  `unstated`. Parenthetical remarks are ignored, text after an arrow is a product, and in "X from Y" the
  precursor is Y.
- Collision regime (`collision_regime`): `in_source` when the collision field names in-source activation and,
  once those phrases are removed, no collision-cell or trap activation, or names only a source voltage (cone,
  fragmentor, orifice); otherwise the regime of the fragmenting instrument: `beam` (QTOF, triple quadrupole,
  Orbitrap HCD), `trap` (ion trap, trap CID on hybrids), `ft_cell` (IRMPD, SORI), `mixed` or `unknown`.
  `in_source_mentioned` marks a record that also names in-source activation (pseudo-MS3). Strata for drafting are
  instrument stratum x collision regime.
- Term rows: one per formula token, polarity and adduct. A `both` record gives one row per polarity, each with
  the adducts of that sign. A charged token keeps only the polarity of its charge; an uncharged token takes the
  polarity the text gives it ('C14H9O2 and C13H9O2 (negative)', a 'positive:' heading of its clause,
  `polarity_from_text`), else both polarities, flagged `polarity_unverified` when the text names polarities
  elsewhere. A record with several adducts of one polarity gives one row per adduct (`expanded_from_adducts`).
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
  formulas into one loss (`CH2O+2H2O`, `RCOOH + NH3`); `formula_isotope_label` for isotope-labelled notations
  (`C10H10(15N)+`), which give no term. The record gets `formula_flags`.
- Roles and grouping: a term row stays in `terms.tsv` but is kept out of groups (`group_exclusion`, `group_id`
  empty) when
  - its token is the precursor or intact molecule (`precursor`): its label or the words just before it say
    precursor, parent, molecular ion, `[M-H]-`, `M+`, `MH+` ... ('[M+H]+. at m/z 369.1686 (C18H27NO7)'), a heading
    of its clause introduces precursors ('Precursor formulas from Table 1: ...'), or it is the ion before an arrow
    in a neutral-loss record ('C24H39O5- 407.2803 -> C23H35O2-'); 'precursor m/z 353: C7H11O6-' names the
    precursor by its number, so the token after the colon is not one;
  - the text names it as not observed (`not_observed`: 'the analogues would be C7H7O4- ...');
  - it is a notation quoted as printed (`printed_as`: 'C7H11O6- (printed as C7H12O6)'), which the extractor does
    not vouch for;
  - its record is an `absence` or `veto` (`role_unverified`: the token may be the absent ion, the observed one,
    an isobar or a competitor, and a reader must say which);
  - its record has kind `manual` (stage-3 `other`), until a reader assigns a kind;
  - it occurs only inside a composite expression (`composite_part`);
  - `record-flags-p3.tsv` flags its record, or this token of it (`record_flag`; the flag name is in `flags`).
  `not_observed` and `printed_as` tokens are also listed in `unmapped.tsv`. A group is `kind_consistent` when no
  member carries `kind_conflict`.
- Masses: element masses are RDKit's most-common-isotope masses, which `tools/ws8_p2/p2_diagnostics.py:73-77`
  calls at run time (RDKit 2026.09.2; the copied values were checked equal, difference 0.0, on 2026-11-18), and
  the electron mass 0.000549 Da is the literal of its line 77. m/z = (sum of masses - z * electron) / |z|.
  C, H, N, O, P, S, Na, K and Cl only; any other element gives `unsupported_element` and no mass.
- Printed values: decimals within 0.5 Da of the computed mass are compared with the product tolerance
  max(0.005 Da, 10 ppm), plus half the last printed digit. A number written right next to a different formula
  whose own mass lies within 0.5 Da of it ('C30H46NaO4+ 493.3288') is that formula's mass and is not compared;
  one-element tokens, composite parts and a ':' between formula and number never pair. Numbers in `value` are
  used first; `note`, which often holds recomputed masses, is used only when `value` has none, so a misprint is
  not hidden by its correction, and a disagreement found there is advisory (`note_mismatch`, `note_lowres`).
  With no match, a value printed with fewer than 3 decimals whose half last digit is not below the tolerance
  gives `lowres_printed`; any other gives `printed_mismatch`.
- Groups keep every record id and citation (DOIs as `doi:10.xxx` in lower case; NPR `primary_refs` as secondary
  citations) and count `n_papers` (distinct `paper_id`, or DOI when there is none) and `n_standards_papers`
  (papers with `standards_based` true). The same term is split over adduct groups (`unstated` included), so
  `term_support.tsv` and the `support_*` fields of each group count the papers of the term over every adduct.
  These counts are evidence, not votes.

## Known limitations

- Family labels map per label. Records under `alkaloid_general` that name one subclass stay `generic` until a
  reader relabels them. Records whose class scope is narrower than their label are moved only by the id list in
  `family-override-p3.tsv`.
- `acylcarnitine` (and its confounder `betaine`) is `in_scope` because the directive names it among the P3 v1
  families, although it is a fatty-acyl class (CAR in the MS-DIAL lipid model); whether "lipids out of scope"
  covers it is the author's decision, and a one-line change in `family-map-p3.tsv`.
- The stratum override in `stratum-override-p3.tsv` is inferred from one record's collision text and awaits a
  reader's check of the paper's methods.
- Token roles read the words around a token; a role stated elsewhere in the record needs a row in
  `term-override-p3.tsv`. An uncharged product ion in a neutral-loss record is not detected (it is flagged
  `odd_electron_neutral` when its formula is odd-electron).

Tests use synthetic records only: `tests/test_ws8_p3_normalise.py`.

# WS8 P3: rule-file format reference (`rule_format.py`)

`rule_format.py` is the stdlib Python reference of the MS-FINDER diagnostic-rule format in the WS8-P3 C# design
(section 3 and the arithmetic of sections 4.1-4.7). It reads and checks rule sets; it never evaluates a spectrum
(that is the C# engine, followed later by a Python port).

```powershell
python -B tools/ws8_p3/rule_format.py validate <RuleSet> [--draft] [--rules <glob>]
python -B tools/ws8_p3/rule_format.py fmt <file>
python -B tools/ws8_p3/rule_format.py sha <RuleSet> [--draft]
python -B tools/ws8_p3/rule_format.py masses <RuleSet> --rule <rule_id> [--draft]
```

- `validate` prints every `RuleLoadError` as `code<TAB>file<TAB>RFC 6901 pointer`, sorted by (file, pointer,
  code), then `OK n_rules=<n>` (exit 0) or `ERRORS n=<n>` (exit 1). `--rules` checks only the matching rule files;
  the tables are always checked.
- `--draft` treats the directory as an unpinned draft: it skips the version pin, accepts a missing
  `tables/motifs.json` and `tables/motif_nesting.json`, and merges in memory, after the table rows and in path
  order, `_motifs/*.json` and `_nesting/*.json` (top-level arrays of rows) and `_requests/*.json`
  (`{"units": [rows], "adducts": [rows]}`). An id defined twice across files is `duplicate_id`. Staging files are
  never hashed, and `_audit/` is ignored. Without `--draft`, a `*.json` under a `_` folder is `unexpected_file`.
- `fmt` rewrites a file in storage form (`json.dumps(obj, sort_keys=True, indent=1, ensure_ascii=True) + "\n"`,
  LF), after checking that the strict reader accepts the result.
- `sha` prints the rule-set sha256, each hashed file's canonical sha256 and each loaded rule's
  `rule_semantic_sha256` (the value a review block must carry).
- `masses` prints, for one rule, every hypothesis (composition x adduct in `adduct_ambiguous` order) with its
  neutral formula and precursor m/z, and every term member target; members within 1.0 Da of the precursor are
  marked `precursor_region`.

Format decisions this reference makes where the design is silent or inconsistent (the C# loader mirrors them):

- **Tables are objects.** Design 3.4 calls `adducts.json`, `units.json`, `motifs.json` and `motif_nesting.json`
  arrays, but its reader refuses a top-level array (`top_level_not_object`) and every `*.json` is hashed through
  that reader. Each is therefore `{"format": "msfinder-diagnostic-<adducts|units|motifs|motif-nesting>",
  "format_version": 1, "<adducts|units|motifs|nesting>": [rows]}`, like `elements.json`.
- **Error codes the design leaves open:** `bad_type`, `bad_value`, `kind_field_choice`, `missing_file`,
  `unexpected_file`, `duplicate_id`, `adduct_polarity`, `adduct_name_mismatch`, `missing_required_unit`,
  `nesting_cycle`, `mode_mismatch`, `unknown_variable`, `series_empty`, `review_not_allowed`, `history_format`.
  `ERROR_CODES` in the module lists every code with its meaning. Reader errors carry the pointer `""`.
- A float field needs a decimal point (`0.0`, not `0`). A canonical adduct name must agree with its multiplier,
  `add` and `remove` (`[M+HCOO]-` has `add` `CHO2`); aliases are free spellings, unique with spaces stripped.
- A `scaffold_ions` key must be present, with a non-empty list for the rule's mode, in every entry of the list.
- A direct term child of `k_of` must be `counted`; a `not` inside a `not`, and a `k_of` inside a `not`, are
  `role_placement`. `than` is required for `ordering` and `veto_stronger_than` and null otherwise; `ratio` needs
  `than`, and `ratio_source` is set exactly when `ratio` is.
- `claim_terms_too_few` does not stop the two-ion passes (design 3.6 item 6 expects it with `claim_single_ion`).
  The precursor region of the load-time check is |target - precursor| <= 1.0 Da + 1e-10 Da (100 pDa).
- `mass_range` also refuses an ion or loss whose mass is not positive (an empty composition).
- `rule_semantic_sha256`'s `elements` covers the symbols written in the formulas the rule resolves; the hydrogen
  that the Y0 forms add or remove counts only when a formula holds H.

Tests use synthetic rule sets only: `tests/test_ws8_p3_rule_format.py` (one single-defect case per error code,
canonical-JSON vectors, masses against hand sums, the two-ion and precursor-target examples of design 3.6,
`parameters_origin` derivation, review staleness, draft merging and the version pin).

# WS8 P3: shared tables (`build_tables.py`)

`build_tables.py` holds the rows it wrote for the draft rule set (`ruleset-draft-20261118`): `ruleset.json`,
`tables/elements.json` (the RDKit masses and electron of `rule_format.py`), `tables/adducts.json` (10 adducts, from
the MS-DIAL adduct resource and, for `[M]+`, two anthocyanin papers), `tables/units.json` (12 units) and
`scaffolds/flavonoid_aglycones.json` (51 aglycones with their A-ring RDA ions). Rows and citations were transcribed
from the normalised candidates (`candidates-20261118`: DOI and locator of each record, the source field of the NPR
records). Before writing, it checks each aglycone against the ions the records give: the NPR `[A+H]+` formulas, the
nominal aglycone ions printed elsewhere, and the accurate m/z of the RDA ions.

```powershell
python -B tools/ws8_p3/build_tables.py write <RuleSet> [--rule-set-version 0.1.0-draft]
python -B tools/ws8_p3/build_tables.py compare <RuleSet>
```

- `write` creates the five files in storage form and prints `write`, `unchanged` or `differs` per file. A file that
  holds other content is never overwritten: the run stops with exit 1 before writing anything.
- `compare` prints, per table, how many built rows are equal, differ or are missing, the rows the builder does not
  produce (`extra`), and whether the built rows still lead the table (the unit order is the canonical order). Exit 0
  when every built row is present and equal; extra rows are allowed.
- **The builder reproduces only its own rows.** The tables were later extended by the merged unit and adduct
  requests (26 units, and `[M+H-2H2O]+`; the request files now sit under the draft's `_staging/_requests`). Those
  rows exist only in the rule set's own `tables/`, and `compare` lists them as `extra`. The scaffold lists of the
  individual families are written by their drafters, not by this script.

Tests use temporary directories and synthetic rows only: `tests/test_ws8_p3_build_tables.py`.

# WS8 P3: C# parity vectors (`cjson_vectors.py`)

`cjson_vectors.py` writes `cjson-vectors.json`, the CPython-oracle vectors the C# port of the rule-file format must
reproduce (C# design sections 2.5, 3.1-3.2, 4.9, 4.11 and 7.2). Everything is synthetic; the file is generated,
never committed, and is itself strict JSON in storage form, so the C# strict reader can read it.

```powershell
python -B tools/ws8_p3/cjson_vectors.py write <dir>/cjson-vectors.json [--source-md]
python -B tools/ws8_p3/cjson_vectors.py check <dir>/cjson-vectors.json
```

| Section | Vectors |
|---|---|
| `strict_json` | input bytes (base64) to the canonical text and its sha256, or the reader's error code: key order by code point, nesting to depth 128, whitespace, escapes, every accepted integer and float form and boundary, and every refused form (BOM, invalid UTF-8, grammar, duplicate keys, NaN, depth 129, integer and float ranges, `-0`, exponents, more than 15 digits, control and non-ASCII characters) |
| `float_lexeme` | accepted lexeme to IEEE bits, with the mantissa and scale of the one division of design 3.1 |
| `fixed_decimal`, `result_double` | double and decimals to exact half-even text (no sign on all-zero text); `null` when not finite |
| `string_escape`, `ieee_bits`, `instrument_class` | UTF-16 text to `ensure_ascii` JSON; NaN canonicalisation; the ASCII-only instrument class |
| `file_hashes`, `ruleset_load`, `rule_semantic` | per-file canonical sha256 and the rule-set sha256 of small file maps; two clean synthetic rule sets with every hash and each rule's `resolved` set; `rule_semantic_sha256` cases |
| `masses` | element lexemes to integer pDa (`mass_precision` above 12 fraction digits), formulas to pDa, ions of both charges with the electron, adduct precursors (multiplier, add, remove), refused formulas |
| `loader_errors` | a clean base file map and, for every code of `rule_format.ERROR_CODES`, a patch (null removes a path) with the exact sorted error list of `load_files(files, draft=False)`, version pin included |

`write` prints the sha256; `--source-md` writes `SOURCE.md` beside the file (generator, catalog commit, interpreter,
sha256; repository-relative paths only). The output is the same bytes for any hash seed and interpreter version
(checked with 3.12 and 3.14). The C# step copies both files into
`tests/Common/CommonStandardTests/Resources/DiagnosticRuleFixtures/` of MsdialWorkbench and does not edit them.
Design 7.2 also plans `engine/loader-errors.json` owned by MsdialWorkbench; the loader vectors here are the CPython
side of the same codes, and the C# step decides whether they replace or complement it.

Tests: `tests/test_ws8_p3_cjson_vectors.py` (stable output under two hash seeds, every error code with a loader
vector that a fresh load reproduces, and the other sections recomputed).

# WS8 P3: rule-level measurement (`p3_measure.py`)

`p3_measure.py` measures every diagnostic-ion rule of the pinned rule set 0.2.0
(`ecf7988dea90992eb7233fc856ce2cb19a8afe4f1f5683628a4b35a459c93569`) on library spectra with known structures,
by running the C# console's `diagrules` subcommand and joining its evidence to structure-derived truth. The frozen
VS20+VS21 pass of 2026-11-24 ran rule set 0.1.0 (`af4a7ebf...b2eb`). Sources: VS20 (development), and the held-out
VS21 and C26 (the compound-disjoint part of CASMI26 `train.parquet`, handback section 22). It reuses P2's pure helpers by import (`p2_diagnostics.py` is not modified) and the
Python rule-format reference (`rule_format.py`) for rule metadata and hypothesis masses.

Stages (positional), all writing into one run directory:

| Stage | Needs | Does |
|---|---|---|
| `prepare` | `--ruleset-dir` | P2's record filters and counters, P1 standardisation, truth motif sets from the shipped 106 motif SMARTS (`tables/motifs.json`), record units, within-source clusters, `prepared.parquet`, and the allowlisted MSP shards `msp/vs20/representatives-NN.msp` with `representatives-index.tsv` |
| `inventory` | `--ruleset-dir` | truth only, after prepare, never the console: `inventory<suffix>.json` and `.tsv` with the prepare counters, exclusions, stats and P2 parity per source, the distinct C26 `instrument_type` strings with their record counts and P2's instrument class, and the power count (for every rule with `measure: true` and `source_checked: true`, per source and pooled in the rule's ion mode, the representatives and clusters of the universe and those whose truth motif set holds the rule's motif); no evidence, no precision |
| `fire` | `--console`, `--ruleset-dir` | the synthetic chimeric set (`msp/vs20/chimeric-NN.msp`, `chimeric-index.tsv`) and the in-source recast set (`insource-NN.msp`, `insource-index.tsv`), then one console run per shard (all 13 passes) at BelowNormal, at most `--workers` at once, into `fire/vs20/<set>/<shard>/`; every manifest is checked (pinned rule set, input sha256 and record count, the 13 passes); all shards must report the same `build_identity` and `engine_semantics` (and, for VS21, the frozen ones) and no record errors in any pass, or `fire` stops; an existing shard is reused only when its manifest validates |
| `measure` | `--ruleset-dir` | joins evidence to truth, counts, intervals, the gate, `fires.parquet` and the aggregates; marks the run complete |
| `all` | both | prepare, fire, measure |
| `freeze` | both, `--freeze` | writes the freeze file the VS21 pass must match |

Options: `--out` (run directory), `--sources vs20` (or `vs20,vs21` or `vs20,vs21,c26`, see VS21 and C26),
`--c26-train` (CASMI26 `train.parquet`; default `<annotation-root>\external\casmi26\...\train.parquet`), `--workers` 1..8 (default 8:
chemistry processes, MSP shards and concurrent console processes), `--limit` (deterministic per-file fraction,
the 1 % smoke), `--seed` (default 20261118; chimeric partners and bootstrap), `--bootstrap` (default 1000),
`--detach` (P2's WMI launch, log in `run.log`), `--vs21`, `--freeze`.

Run directories are `<annotation-root>\ws8-p3\p3cs-<yyyymmdd>-<sha>` (beside P2's `ws8-p2`) and
`p3cssmoke-<yyyymmdd>-<sha>`. A smoke name needs `--limit` and a full name refuses it; `<sha>` must be a prefix of
the checked-out commit; the tree must be clean and the commit pushed; a completed run is immutable, and later
stages must use the parameters recorded by the first.

```powershell
$py = '<annotation-py312>\Scripts\python.exe'   # the RDKit environment P2 uses
# 1 % smoke, from the pin's root
& $py -B tools/ws8_p3/p3_measure.py all --sources vs20 --workers 8 --limit 0.01 `
  --ruleset-dir <MsdialWorkbench>\src\Common\CommonStandard\StructureFinder\DiagnosticRules\RuleSet `
  --console <MsdialWorkbench>\src\MSFINDER\MsfinderConsoleAppCore\bin\Release\net472\MsfinderConsoleApp.exe `
  --out <annotation-root>\ws8-p3\p3cssmoke-<yyyymmdd>-<shortsha> --detach
# full VS20 pass: the same without --limit, into p3cs-<yyyymmdd>-<shortsha>
# held-out pass with C26: freeze first, then prepare, inventory (truth only), fire and measure
& $py -B tools/ws8_p3/p3_measure.py prepare --sources vs20,vs21,c26 --vs21 `
  --freeze <annotation-root>\ws8-p3\freeze-<yyyymmdd>-<sha>.json --ruleset-dir <RuleSet> `
  --console <MsfinderConsoleApp.exe> --c26-train <CASMI26>\train.parquet `
  --out <annotation-root>\ws8-p3\p3cs-<yyyymmdd>-<shortsha>
& $py -B tools/ws8_p3/p3_measure.py inventory --sources vs20,vs21,c26 --vs21 --freeze <the same> `
  --ruleset-dir <RuleSet> --console <MsfinderConsoleApp.exe> --out <the same run directory>
```

`inventory` needs `--console` only to check the freeze (it hashes the bin folder and never runs the console).

## Outputs

| File | Content |
|---|---|
| `prepared.parquet` | one row per representative: record id, first block, Hill formula, motif ids, cluster id, mode, curated adduct, strata, provenance bucket, peaks (local row file) |
| `msp/`, `fire/` | the console inputs, index tables (file name, record index, record id) and evidence (local) |
| `fires.parquet` | one row per (set, level, record, pass, rule) where the rule had a hypothesis: `n_hypotheses`, `fired`, `veto_blocked`, the chosen hypothesis' adduct, `motif_match` |
| `rules_precision.csv` | per (rule, version, ion mode, pass, source, instrument class, energy bin): universe and conditional counts, records and clusters, precision, Wilson, recall within the precondition, coverage, LR+, prior stress, the conditional stress grid, cluster-bootstrap intervals, and for the base pass the decoy fires and FDR in records (`decoy_*`, report only) and clusters (`c_decoy_*`, gated) |
| `rules_competitors.csv` | per rule and competitor motif: conditional competitor records and clusters, fires, FPR |
| `rules_synthetic.csv` | per rule and synthetic set (`set`): each chimeric level, and each in-source transformation row (`transform`, `change`, `product_adduct`) plus their pool (`transform` = `all`, gated): records not evaluated by the console, fires, TP/FP, precision and Wilson in records (report only) and clusters (gated) |
| `rules_occupancy.csv` | per rule and precondition mode: both occupancy definitions and `non_conditional` |
| `rules_precision.json` | label `ws8-p3-class-rules`, criteria, per rule: parameters, status, every criterion with its value, the pooled row of every pass, synthetic, occupancy, competitors, provenance breakdown |
| `review_proposal-<run-id>.json` | `passed_strict` (tier-1 numbers) and `passed_provisional_only` (tier-2 numbers) as separate lists, both empty in a VS20-only run, and the `incomplete_strict` and `incomplete_provisional` rules; the author chooses which tier to accept |
| `report.md`, `run.json` | summary; code commit, pin path, inputs and sha256, rule-set files (raw and canonical sha256), console path and sha256, parameters, exclusions, P2 parity, file hashes |

Aggregates carry ids only (`rule_id`, `motif_id`); `measure` and `inventory` scan them as P2 does
(`audit_aggregates`).

**Naming.** One rule (`suffix()`): every file derived from VS21 carries `.local-vs21`, every file derived from C26
carries `.local-vs21.local-c26` (C26 records are filtered against the VS21 keys of the same run), and the pooled
outputs of a run carry the tags of all its sources, so a `vs20,vs21,c26` run writes
`rules_precision.local-vs21.local-c26.csv`, `fires.local-vs21.local-c26.parquet`,
`inventory.local-vs21.local-c26.json` and so on; VS20 files have no tag. The console tags its own files only with
`.local-vs21` (it does so when "vs21" is in the input path, as it is in the C26 shard names); C26 shard folders
under `fire/c26/` carry `.local-c26` in their names.

## Criteria

Each rule is gated twice, on the same primary-mode pooled cell, counting clusters. Tier 1 is the table registered
in MsdialWorkbench PR #908 (`CRITERIA`). Tier 2 (`CRITERIA_PROVISIONAL`) was registered by the author on
2026-11-24, before any VS21 pass, after the VS20 development pass had been seen. Counting unit, primary pass,
eligibility and the pending logic are the same in both tiers.

| Criterion | Tier 1 (strict) | Tier 2 (provisional) |
|---|---|---|
| Fired clusters (pooled) | >= 100 | >= 30 |
| Cluster Wilson lower bound | >= 0.95 | >= 0.90 |
| Per-source point precision, sources with >= 30 clusters | >= 0.95 | >= 0.90 |
| Prior stress at half the library-wide prevalence | >= 0.90 | report-only, not gated |
| Cluster decoy FDR (+5, +9, -7 Da), decided exactly | <= 0.05 | <= 0.05 |
| Occupancy (axis and records) | <= 5 % | <= 5 % |
| Synthetic chimeric and in-source cluster precision | >= 0.95 | >= 0.90 |
| `tuned_vs20` rules | the tier's Wilson bound and prior stress on VS21 alone | the tier's Wilson bound on VS21 alone |
| `tuned_pooled` rules | never pass (external set) | never pass (external set) |
| `as_published` rules | gated on the run's sources pooled (VS20 + VS21, or VS20 + VS21 + C26) | the same |

Both tiers are written to `run.json` and the freeze file; the VS21 pass refuses a change in either.

## Status

Per tier: `not_measured` (outside `measure: true` and `source_checked: true`, or refused by the console);
`not_estimable` (fewer fired clusters than the tier's minimum in the primary pass); `fail` (a gated criterion
failed); `incomplete` (every gated criterion this run can measure is met, others pending); `pass` (every gated
criterion met). The combined label is `pass_strict` when tier 1 passes, `pass_provisional_only` when only tier 2
passes, `incomplete_provisional` when tier 2 is `incomplete` while tier 1 failed or was not estimable, and
otherwise the tier-1 status. `rules_precision.json`, the gate row of `rules_precision.csv` and `report.md` give
both statuses with their failed and pending lists. The held-out criterion is met only with VS21 in the run, so a
VS20-only run has no `pass` in either tier.

## Decisions where the specification was silent

- **Record unit.** The curated precursor type is the rule set's canonical adduct when its alias table
  (`tables/adducts.json`, spaces stripped, as the console maps it) knows the spelling; otherwise the stripped
  spelling. Unmapped records stay in the universe and are counted (`adduct_unmapped_records`); the console skips
  them in `mz_only_curated`, where their precondition therefore does not hold. Charge >= 2 (the adduct table, or a
  `]n+`/`]n-` suffix) and the declared in-silico tag are excluded after P2's counters and counted in `p3_exclusions`.
- **P2 parity.** The P2-named counters are compared with P2's `run.json` `exclusions` for every source of a full
  run (a mismatch stops the run); a smoke records `not_checked`.
- **Reader.** P3 carries a copy of P2's `iter_msp` that also reads the allowlisted `COMMENT` and reduces it to a
  provenance bucket at once; a test checks it yields exactly what P2's reader yields.
- **Provenance.** `provenance-p3.tsv` declares ordered regexes on `COMMENT`. A present but unmatched field is
  `other`, an absent one `unknown`; if no record is bucketed the report says `provenance_unavailable`. The
  in-silico detector is the declared tag only (no exact-mass detector). The buckets were checked against the
  VS20 origin counts on 2026-11-23 (aggregate counts only): every VS20 record has a COMMENT carrying
  `origin=<source DB>`; the patterns give massbank 31,765, gnps 14,830 and riken 9,078 records in positive mode and
  massbank 13,714, gnps 3,669 and riken 7,602 in negative mode; the large origins MassSpecGym, BMDMS-NP and the
  Fiehn libraries fall into `other`; `unknown` is 0; no origin looked in-silico.
- **Per-source criterion.** It uses the library sources (VS20, VS21 and, in a three-source run, C26), as P2 does. Provenance buckets are reported
  per rule (`by_provenance`, `provenance_flag` when a bucket with >= 30 fired clusters is below 0.95), not gated.
- **Primary pass.** A rule's declared `primary_precondition_mode` picks the pass that decides its status:
  `adduct_ambiguous` is `base`, `mz_only` is `mz_only_curated`. Decoy FDR always comes from the base-pass
  population, as the decoy passes run without adduct evidence.
- **Counting unit of the decoy and chimeric criteria.** PR 908 makes clusters the counting unit, so the gate
  reads the cluster-level decoy FDR (mean over the three shifts of decoy-fired clusters / fired clusters) and the
  cluster-level chimeric precision. The record-level values stay beside them as report-only columns, and the
  JSON criteria carry them as `record_level`.
- **Populations.** The universe is every representative of the ion mode; records a pass did not evaluate count as
  precondition-false (`n_not_evaluated`). The conditional population is the records where the rule has at least
  one hypothesis (the console's compact result lists exactly those).
- **Clusters.** A cluster fires when any member fires; truth is uniform within it. The gate's n is fired
  clusters, and prior stress is P2's `statistics()` on the cluster universe table. With zero false positives the
  zero-cell correction means prior stress >= 0.90 needs roughly 145 true-positive clusters, more than n >= 100.
- **Bootstrap.** Conditional-population clusters are resampled (multinomial by conditional binomials) with a
  seed string per cell; intervals for record- and cluster-level precision, for every stratum of `base` and
  `mz_only_curated` and the pooled cells of the other passes.
- **Occupancy.** Axis: the fraction of the nominal-mass span of the rule's hypothesis precursors,
  `[floor(min), floor(max) + 1)` Da, within the rule's precursor tolerance of a hypothesis (common adducts
  included for `adduct_ambiguous`, the rule's adducts only for `mz_only`). Records: the fraction of motif-negative
  representatives of the mode whose precondition holds. Either above 5 % makes a rule `non_conditional`.
- **Chimeric set.** Every representative gets one seeded partner (same mode, another first block, precursor
  within +-12.5 Da), used at all three levels; the partner's peaks are scaled so its highest peak is 10, 30 or
  50 % of the primary's base peak (precursor +-1 Da excluded) and merged; the primary's precursor and adduct are
  kept. Every record is mixed, so the "motif-negative records first" order does not change the set. A level
  without fires does not fail the criterion; an empty chimeric set leaves it pending. The console runs all 13
  passes on the chimeric shards (it has no pass selection); only `base` and `mz_only_curated` are read.
- **In-source recast set.** `insource-p3.tsv` (parent_motif, parent_adduct, change, product_adduct, source,
  note; rows awaiting the author's confirmation) defines the transformations; the row number in file order is the
  transformation id. A representative is recast by a row when its curated adduct is the row's parent adduct (alias-
  mapped), its ion mode is that adduct's polarity, its truth motif set holds `parent_motif` (any record for `*`),
  and its precursor ion (truth Hill formula x multiplier, plus the adduct table's added atoms, minus its removed
  ones) holds every atom of the loss; the truth formula is used for this check only and never reaches the console.
  The new precursor m/z is the observed one minus the exact loss mass from the rule set's `tables/elements.json`.
  Peaks above the new precursor + 0.01 Da are dropped; peaks within 0.01 Da of it are replaced by one precursor
  peak at the base-peak intensity of the kept peaks; a record with no kept peak is counted (`no_peaks_left`) and
  dropped. PRECURSORTYPE is the product adduct (canonical when the alias table maps it, else the stripped spelling,
  counted as `adduct_unmapped`; the base pass ignores it, `mz_only_curated` skips it). Truth stays the parent's.
  Losses are read in any element order (`NH3`, `SO3`, `H3PO4`), each symbol once, every element in the element
  table; a row whose parent motif or parent adduct is unknown, or whose product adduct has the other polarity,
  stops the run. The set is built deterministically; above 500,000 records it is cut to a seeded sample of
  500,000 and `run.json` says so. The gate reads the pooled cluster precision (>= 0.95, or no fires); a mode
  without recast records leaves the criterion pending. The recast approximates the product's own MS/MS, which the
  library does not hold.
- **Not evaluable.** A criterion whose value cannot be computed (no fires) is pending, never passed or failed;
  such a rule is `not_estimable` anyway.
- **Strata.** Instrument class and energy bin are P2's `instrument_bin` and `energy_bin`, with `all` margins; the
  console's own instrument class is cross-tabulated in `run.json` and the report.
- **Console.** `diagrules.ini` sets MS-FINDER's product tolerance to 0.01 Da (every measured rule uses its own
  `p2` tolerance). Representatives are split into `--workers` contiguous MSP shards; the join key is (shard file,
  record index) through the index table, and `NAME` carries only the record id. `fires.parquet` has the chosen
  hypothesis only (compact console results), not one row per hypothesis. `fire` stops when any pass of any shard
  reports record errors; skips are allowed and `report.md` lists them per source, set and pass with their reasons.
  `run.json` records the sha256 of every `.exe`, `.dll` and `.config` in the console's bin folder and the build
  identity the shards share.
- **VS21.** `--sources vs20,vs21` runs only with `--vs21` and a freeze file
  (`ws8-p3\freeze-<yyyymmdd>-<sha>.json`, written by `freeze`; the name and folder are checked) whose content
  matches now: the sha256 of every RuleSet file, of the console and of every `.exe`, `.dll` and `.config` in its
  bin folder, of the measurement code (`p3_measure.py`, `rule_format.py`, `p2_diagnostics.py`, `mcs_core.py`,
  `provenance-p3.tsv`), of `insource-p3.tsv`, and the parameters `seed`, `bootstrap` and the console INI text. `freeze` also runs the
  console once on a synthetic record and stores its `build_identity` and `engine_semantics`, which every VS21
  shard must match. VS21 needs VS20 in the same run (deduplication by first block). Its row files and every
  aggregate are named `*.local-vs21.*`. The freeze also pins the C26 parameters (below).
- **C26.** `--sources vs20,vs21,c26` runs behind the same `--vs21` flag and freeze. CASMI26 `train.parquet` is
  identified by its sha256 (`9423f90e...bcd98`); it is CC BY-NC 4.0 with a no-redistribution rule, so no row, key,
  SMILES or spectrum leaves the `*.local-c26.*` files, and tests use synthetic parquet files. `iter_c26` reads
  batches of a column subset and yields the tuples `iter_msp_p3` yields, once per ion mode (rows whose
  `ionization_mode` is neither positive nor negative come with the positive input and fail P2's mode check);
  `record_id` is the parquet row ordinal x 2 (+1 in negative mode). Before P2's checks, and counted in
  `p3_exclusions`: `c26:excluded_ingest_lib` (enveda-180), `c26:missing_structure` (no `normalized_smiles` or
  `inchikey14`), `c26:precursor_error_ppm_missing` and `c26:precursor_error_ppm_above_limit` (|ppm| > 10). Spectra
  merged from several collision energies are kept, counted (`c26:merged_collision_energy_records`) and form the
  energy stratum `merged`. Records whose structure key (the same standardisation as VS20/VS21) is among the VS20 or
  VS21 keys of the run are excluded (`c26:vs20_overlap_records`, `c26:vs21_overlap_records`). `ingest_lib` is the
  provenance bucket (report only, `c26:<lib>` in `by_provenance`). Instrument strings are classified by P2's
  `instrument_bin`; prepare writes them with counts to `instrument_types.local-vs21.local-c26.tsv` (never to
  `run.json`) for `inventory`. P2 parity is `not_checked` (there is no P2 run). The freeze pins the train sha256,
  the ppm limit, the excluded libraries and the classifier.
- **Held-out small cells.** Only cells of a held-out source (VS21, C26) are suppressed, and only where their own
  denominator is below 5: a denominator of 1-4 (fires, fired clusters, the conditional population, the universe,
  positives, the occupancy negatives) is blanked with the counts it bounds and every rate over it, and the
  likelihood-ratio and stress values of that row go with it. VS20 cells, pooled cells, and the competitor,
  synthetic and occupancy rows (pooled over sources) are never blanked; pooled criterion values stay in the JSON,
  and per-source criterion values of a held-out source with 1-4 fired clusters are blanked. The count tables of
  `report.md` blank held-out counts of 1-4 only (instrument classes are pooled and kept). A pooled value minus the
  VS20 value can give a held-out value back; the rule accepts that.

## Not done in this phase

Cross-library VS20/VS21 cluster matching by cosine, the all-spectra weighted sensitivity variant, `as_applied`,
the MS-FINDER-level measurement, the runtime measurement, conditional mining, the precondition decoy and the
figures. `report.md` and `rules_precision.json` list them too.

Tests use synthetic spectra, public molecules and hand-built counts only: `tests/test_ws8_p3_measure.py`
(the console integration test runs when `WS8_P3_CONSOLE` and `WS8_P3_RULESET` point to the console and its
RuleSet folder).
