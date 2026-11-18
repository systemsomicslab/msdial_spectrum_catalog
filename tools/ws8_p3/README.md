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
