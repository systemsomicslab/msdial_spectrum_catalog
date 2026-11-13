# WS8 P2: diagnostic-fragment SI evidence

This standalone experiment measures whether a sourced product ion or neutral loss
indicates the linked SMARTS motif in library compounds. It imports P1 chemistry
helpers, leaves P1's motif table unchanged, and adds no package runtime dependencies.
Use the existing D5 interpreter; do not install packages into it.

```powershell
& 'D:\13_MSDIAL_Public_Reanalysis\envs\annotation-py312\Scripts\python.exe' -B tools/ws8_p2/p2_diagnostics.py all --sources vs20,vs21 --workers 8 --limit 0.01 --out F:\MSDIAL_Public_Reanalysis\annotation\ws8-p2\p2smoke-20261113-<shortsha> --detach
# After smoke validation, omit --limit and use p2-20261113-<shortsha>.
```

Stages are `prepare`, `observe`, `precision`, `mine`, `choose`, and `all`. Run
individual stages in that order in the same directory with the same arguments.
`--sources` accepts `vs20`, `vs21`, or `vs20,vs21`; VS20 is always prepared to
establish the complete deduplication universe, even in VS21-only mode. Mining
requires both sources for held-out validation; missing held-out support is explicit.
`--limit` is a deterministic per-file fraction, selecting every round(1/fraction)th
record. All raw structures are still standardized for the VS20 overlap exclusion.
The full run uses every record. Completed runs are immutable.

`--detach` uses the built-in WMI Win32_Process service to launch independently of
nested application jobs, writes the run's `run.log`, and returns its child PID.
The child checks it is outside every Windows job. A direct Win32 breakaway launch
was insufficient on this host and exited before reading any inputs. The run
process and each spawned chemistry worker set and verify BelowNormal via
Win32. At most eight workers are permitted, and BLAS/Numba thread variables are one
before optional imports. No ports or campaign state are accessed.

The source table contains 32 rows from the spec's starting diagnostics. Every mass
is the sum of RDKit's most-common-isotope masses, plus 0.000549 Da for a negative
ion or minus it for a positive ion. Neutral losses have no electron correction.
The TSV has the actual formula-derived values, not the rounded examples in the
spec. Sources are public papers, reviews and the public MS-DIAL lipid rules.
For phosphatidylserine, the 185-Da loss is positive mode and the 87-Da serine loss
is negative mode, following the public rule. Only explicit precursor types in the
table are eligible. Multiply charged precursors are excluded from curated neutral
losses because their m/z loss is charge dependent.

The extra motif table includes O-pentopyranoside (C5 is `[CH2;R]`), carnitine,
N-acyl taurine and adenine. The motifs describe the displayed substructures, not
all possible members of the named class. Sugar SMARTS require the displayed free
hydroxyls and ignore stereochemistry. N-acyl taurine excludes free taurine. The
tests cite hand-picked PubChem examples and show both positive and negative
matches after P1 standardization. Both P1 pentoside and P2 pentopyranoside are
matched and mined; only the stricter motif is linked to the curated 132-Da loss.

Preparation excludes records with no usable peaks, invalid precursor m/z, mismatched
mode, missing precursor type, nominal peaks (every m/z within 0.001 Da of an
integer, or every lexical m/z having fewer than two decimal places), or rejected
structures. The representative is the eligible record with most positive finite
peaks per recomputed first block, mode and precursor type; ties keep the earliest
record. VS21 is removed by first block against *all* standardizable VS20 structures,
including nominal or otherwise spectrally excluded records. Exclusion counts
distinguish record units from unique raw-SMILES structure rejections.

The base peak excludes precursor +/-1 Da. Diagnostic tolerance and intensity
thresholds are inclusive. Neutral losses use the recorded PRECURSORMZ and require
the matched peak to be at least 1 Da below it. Mode-level pools sum the specified
compound/mode/precursor-type units across eligible types. Individual precursor
types are also reported. The all-spectra sensitivity variant gives each eligible
spectrum weight 1/(eligible spectrum count for that unit). Its fractional Wilson
and log intervals are descriptive, not confidence intervals for independent
replicate spectra.

Statistics retain TP/FP/FN/TN: uncorrected precision, Wilson 95%, recall, FPR,
prevalence, LR+ and log 95% interval. When any cell is zero, add 0.5 to all four
cells for LR+ and its interval only (Haldane-Anscombe correction). LR+ and prior
stress are unavailable if a positive or negative population is absent. Prior
stress uses half the uncorrected prevalence and the LR+ lower bound. The rule is
n_obs >=100, pooled Wilson lower >=0.95, precision >=0.95 in every source with
n_obs >=30, and prior stress >=0.90; insufficient n is `not_estimable`.

Variants change one parameter at a time: tolerance 0.005/0.02 Da, intensity 5/10%,
all-spectra weighting, formula gating, instrument type, and collision energy.
The report-only formula gate requires at least the diagnostic residue's counts of
each heavy element; hydrogen is ignored because ionization and adducts change its
stoichiometry. Instrument bins are TOF, Orbitrap/FT, other. Collision energies must
be a single parseable value; explicit eV/V, NCE/% and unspecified units are kept
separate, each in <20 / 20-<40 / >=40 bins. Ranges and stepped fields remain
unparseable. The P1 cross-check includes every family from its aggregate table.

Mining first collapses representatives to one compound per mode, OR-ing bins across
single-charge precursor types. It uses half-up 0.01-Da bins and >=1% intensity;
neutral losses must be >=12 Da. A motif needs >=200 VS20 positives. Bins present in
>=5% of positives are ranked by VS20 precision, then positive support, kind, and
mass bin; the top 20 combined product-ion/loss bins are frozen before VS21 is read.
VS21 evaluation uses these exact bins and the same thresholds, with n_min=30.
Formula assignment searches nonnegative CHNOPS compositions within 3 mDa of the
bin centre (electron correction for ions), without a valence heuristic; it stops
after two formulas prove ambiguity. Mined candidates never enter the allowed list.

All output is beneath the validated run directory. `prepared*.parquet` holds
eligible spectra and motif matches; `observations*.parquet` holds representative
diagnostic observations. Every VS21 row artifact uses `.local-vs21.` in its name.
`confusion_counts.csv`, `diagnostics_precision.csv/.json`, `mined_selection-vs20.json`,
`mined_candidates.csv`, `p1_crosscheck.csv`, `report.md`, and the figure contain
aggregates only. `run.json` records input names and checksums (never private paths),
commit, versions, parameters, process/priority checks, exclusions, stage times,
and hashes of all artifacts except itself and the open log. Report its own hash
externally. `audit_outputs` recomputes JSON statistics and decisions from CSV counts
and scans aggregates for identifiers/forbidden columns (P2-A4/A5).

Validation:

```powershell
$env:PYTHONPATH='src'
& 'C:\Users\Hiroshi Tsugawa\AppData\Local\Python\pythoncore-3.14-64\python.exe' -B -m unittest discover -s tests
& 'D:\13_MSDIAL_Public_Reanalysis\envs\annotation-py312\Scripts\python.exe' -B -m unittest tests.test_ws8_p2 -v
```

P2-A1..A3/A6/A7 use synthetic spectra and cited public molecules. P2-A4 runs the
actual CSV-to-JSON choice path with synthetic counts and detects tampering. Full-run
A4/A5 are checked automatically before completion; the test module never reads a
real library or row artifact.
