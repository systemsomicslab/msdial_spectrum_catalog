# WS8 P1 SI substructure-size experiment

Standalone local research CLI; the catalog package and its dependency list stay unchanged.
Run with the existing D5 Python 3.12 environment. Install nothing. Inputs are read-only.
Only source, this document, motif definitions and synthetic/public molecule tests belong in Git.

```
<D5-python> -B tools/ws8_p1/p1_si_bar.py all --pilot --sources vs20,vs21 --workers 16 --seed 20261113 --out <assigned-root>/p1pilot-<date>-<commit> --launch
<D5-python> -B tools/ws8_p1/p1_si_bar.py all --sources vs20,vs21 --workers 16 --pilot-run <pilot-output> --out <assigned-root>/p1-<date>-<commit> --launch
```

Stages: `prepare`, `neighbours`, `sets`, `mcs`, `curves`, `choose`, `all`.
Each stage updates run.json. Resume a failed stage explicitly, then run subsequent stages.
`--reuse-prepare <run>` copies the prepared immutable inputs and selected truths for a pilot rerun;
neighbours, sets and all measurements are recomputed. Use a different valid pilot run directory.
`--launch` requests Windows job breakaway, detached hidden execution, and BelowNormal priority;
failure to break away stops the launch. All workers verify BelowNormal independently and use one
numba/BLAS thread. No ports, services, ledgers or campaign processes are accessed.

Standardisation and the MCS identity live in mcs_core.py for reuse by WS8. The DB checksum is verified
before use. DB columns are positional. All structures use recomputed constitution keys and formulas.
Library representatives maximize usable positive-intensity peak count; m/z values are sorted and
raw positive finite intensities are retained, with no additional peak filtering or transformation.
The installed matchms 0.34.0 names its greedy modified cosine `ModifiedCosineGreedy`; its `pair`
method supplies the score and match count (tolerance 0.02). Both relaxed variants retain the six-match
minimum. Candidate ranking ties break on constitution key. Each source searches its own library.
VS21 truth exclusion uses all usable standardized VS20 truth constitutions, before truth sampling.

The candidate sets are nested in k. Random sets incrementally draw through SeedSequence children for
each k; missing exact-formula entries are supplemented by nominal mass 0 then +/-1. Each row records
relaxation. Scaffolds preserve the original nearest-neighbour order; an empty or truth-identical
scaffold makes the corresponding set unavailable. Scaffold member keys identify the actual scaffold.
Motif oversamples are additional library truths with usable spectra, deduplicated across motifs,
excluded from size/chance curves, and included in motif precision. Unconditional chance formulas
are drawn uniformly within heavy-atom bins before their members are drawn uniformly.

Timeouts discard RDKit's partial result. Empty and timed-out MCS are excluded from size curves.
Motif matching is independent of MCS completion and includes all constructed sets. Curve bins and
cumulative thresholds are separate: atomic point bins are 1..30, 31..40, 41..60 and 61..80;
atomic cumulative thresholds are every integer 1..80. Fractions use 0.05 steps. Quantiles use linear
interpolation. Curves include pooled, marginal and fully crossed source/lipid/heavy-atom cells.

A bar requires both spectral_nn and same_formula_db_random to pass the target Wilson lower bound,
n minimum, and strict matched-chance95 criterion. The maximum of their first passing thresholds is
rechecked in both families; if confidence becomes nonmonotone, the next common passing threshold is
used. If an axis has no passing value it is null; partial estimability is explicit. Global axes are
the maxima over estimable k, with missing k listed; no value is invented for an unestimable axis.
The evidence objects identify family-specific support on the written cumulative curves.
Truth-fraction diagnostics show thresholds through 1.0; the final point bin also counts values
above 1.0, for which containment is impossible. Library structure rejection counts refer to
unique raw-SMILES/adduct representatives; unusable-spectrum counts refer to input records.

All real-data row files stay beneath the assigned run directory. Combined scratch files and every
VS21 row artifact have `.local-vs21.` in the filename. Deterministic sets artifacts exclude timing;
MCS timing summaries live in run.json. Aggregate files contain no structures or identifiers.
run.json fingerprints every artifact except itself and the still-open run.log; its own checksum
is reported externally. Run logs contain stage counts and timings only. The full run chooses cap
600 when the pilot projection is at most six hours, otherwise 300, recording the pilot decision.

The motif SMARTS are explicit public functional-group/ring definitions, not universal class labels.
Sugar motifs cover pyranosides, ignore stereochemistry, and require the displayed free hydroxyls.
Head-group motifs include the glycerol-phosphate linkage. Their measured precision, not their names,
determines whether they may bypass the size rule. The file is suitable for byte-for-byte reuse.

Tests: the repository unittest suite uses its Python 3.14 interpreter (chemistry skips without RDKit);
`tests.test_ws8_p1` uses D5. P1-A4 additionally requires a pilot rerun; compare deterministic sets bytes
and curve bytes, and list timeout differences if any. P1-A6 shortfalls appear in the chance CSV.

After a completed full run and its pilot rerun, use the lightweight audit (no spectral scoring or MCS):

```powershell
& 'D:\13_MSDIAL_Public_Reanalysis\envs\annotation-py312\Scripts\python.exe' -B -m tools.ws8_p1.audit_p1 --out <full-run-directory> --pilot <pilot-directory> --rerun <pilot-rerun-directory>
```

The audit checks frozen artifact hashes, source splits, truth exclusion, timeout rejection, every
written correctness/chance statistic against measured rows, and per-k/stratum/global threshold choices
against the CSVs. It recomputes motif counts and pilot curve differences. Supplemental aggregate reports,
the SI checksum sidecar, and `claim-rules-v1.source.json` stay in the run directory. `finalization.json`
records their hashes and the audit code identity; the measurement `run.json` and `si_bar.json` remain
unchanged. A timeout exception is reported explicitly when measured pilot files differ. A null size
bar remains unavailable to the claim builder; it must never be interpreted as a zero threshold.
