# Running MS-Emblator over many units

MS-Emblator (systemsomicslab/msemblator) predicts formulas and structures from MS/MS with SIRIUS,
MS-FINDER, msbuddy and MetFrag, and ranks them with a learned model. In mode 2 its SIRIUS structure step
also runs CANOPUS. It runs outside MS-DIAL, and its results rejoin a unit's alignment by unit and
alignment ID. Three modules carry a run:

| Module | Runs with | Does |
|---|---|---|
| `msemblator_export` | the catalog's Python | selects a unit's representative MS/MS and writes identity-free MSP, a manifest and a selection record |
| `msemblator_batch` | the catalog's Python | exports and chunks many units (`plan`), runs the driver per chunk (`run`), summarises (`status`, `throughput`) |
| `msemblator_driver` | **MS-Emblator's** Python (its `.venv`) | runs MS-Emblator on one chunk with SIRIUS's stored login and keeps everything the run produced |

## Before the first run

1. Install MS-Emblator and its runtime assets as its README says.
2. Log in to SIRIUS once, in your own terminal, with the SIRIUS that MS-Emblator bundles:
   `<msemblator>\external\sirius\sirius.exe login -u <account> -p`. `-p` prompts for the password.
   SIRIUS stores one login per version, so a login to another SIRIUS install does not count.
   Never pass `--sirius_pass`: a password on a command line is visible in the process list.
3. Check the login without printing it:

   ```
   python -c "from msdial_spectrum_catalog.msemblator_driver import stored_login; stored_login(r'<msemblator>\external\sirius')"
   ```

   It is silent when a login is stored and raises otherwise. It never prints `login --show`, which names
   the account.

## Plan

```
python -m msdial_spectrum_catalog.msemblator_batch plan --root <annotation root>\msemblator --run-id <run-id> \
    --chunk-size <N> (--devset <devset.json> | --units <unit_key>=<unit output dir> ...)
```

For each unit, `benchmark` (strict reference matches) and then `unresolved` are exported to
`<run>\input\<unit_key>\`, and each MSP is split into chunks of at most N records in alignment-ID order:
`<unit_key>_<selection>_c000.msp`, `c001`, and so on. The benchmark manifest is named
`<unit_key>_benchmark_manifest.local-vs21.tsv`: it holds the private library's identities and stays
local. `plan.json` is written once and lists, in run order (benchmark first, the largest unit last in each
selection):

- `exports[]`: `unit_key`, `selection`, `msp`, `source_msp_sha256` (equal to the selection record's
  `msp_sha256`), `selection_json`, `manifest`, `records`, and the tier counts;
- `chunks[]`: `chunk_id`, `unit_key`, `selection`, `chunk`, `chunk_sha256`, `records`, `manifest`,
  `output_dir`, `order`.

Paths in `plan.json` are relative to the run folder. The plan refuses a run folder inside the reanalysis
root (any folder holding `_campaigns\`): unit outputs are read-only inputs.

**Chunk size.** MS-Emblator gives the whole SIRIUS structure workflow one 1800 s timeout, whatever the
number of compounds. Measure `t_struct`, the SIRIUS structure seconds per spectrum, on a smoke run and use
`N = floor(0.5 x 1800 / t_struct)`, capped to [20, 300]; `throughput` computes it.

## Run

```
python -m msdial_spectrum_catalog.msemblator_batch run --plan <run>\plan.json --msemblator-root <msemblator>
```

- One SIRIUS at a time: the runner takes `<root>\.sirius.lock` (pid, host, run id, and the pid of the
  driver it last started) and refuses while a live process holds it. The lock is stale only when
  neither the runner nor that driver is alive. **The lock covers only runners that use the same lock
  file:** two annotation roots, or two MS-Emblator installs, each with its own default lock, can run
  SIRIUS at the same time against one login. Keep one root per machine, or give every runner the same
  `--lock`. MS-Emblator's `work\` folder is shared state, and the driver also refuses a non-empty one.
- Each chunk runs `<msemblator>\.venv\Scripts\python.exe <this checkout>\src\msdial_spectrum_catalog\msemblator_driver.py`
  into `<run>\units\<unit_key>\<selection>\cNNN\`, at below-normal priority (charter C8). SIRIUS's Java
  and MS-FINDER inherit the priority class from the driver.
- Resume: a chunk whose `msemblator_run.json` says exit 0 is skipped, after its `input_sha256` is checked
  against the plan's `chunk_sha256`. A mismatch stops the batch, and a `plan.json` whose sha256 differs
  from the one `run.json` recorded is refused. A failed chunk is moved to the next free
  `cNNN.failed-<n>` and run again, up to `--max-attempts` (2) in all. An interrupted attempt counts. An
  attempt that failed on login or the web service does not count, and is not retried at once: after a
  re-login, a new invocation runs it again.
- Stops: after `--stop-after` (3) consecutive chunks failing on login or the web service; when the driver
  refuses (exit 3); when `<run>\STOP` exists; after `--limit` chunks.
- Exit codes: 0 when every chunk the runner reached is done; 5 when it went through but some chunks are
  left failed (`run.json` lists them under `chunks_failed` and, for chunks whose attempts an earlier
  invocation used up, `chunks_exhausted`); 1 when it stopped (STOP file, refusing driver, consecutive
  service failures, interrupt, a record for another input); 3 when it refused to start (lock held, plan
  changed). An unexpected exception is recorded as the outcome `crashed: <type>: <message>` and
  re-raised.
- Interrupts: on Ctrl+C or Ctrl+Break the runner kills the running chunk's process tree, moves
  MS-Emblator's `work\` into that chunk's bundle, appends an event with `exit_code: null`, writes the
  outcome `interrupted during <chunk_id>` to `run.json` and frees the lock before anything slow (the usage
  limits are not queried after an interrupt). A restart sets the partial bundle aside as a failed attempt
  and runs the chunk again. Prefer the `STOP` file: it ends the batch between chunks and loses nothing.
- `run.json` records each invocation: the catalog and MS-Emblator commits (and whether either checkout was
  dirty), the interpreters, the parameter-file and library sha256, the affinity mask, the driver priority,
  `max_attempts`, the usage limits SIRIUS reports (the allowlisted keys only), the start, end and outcome,
  and the failed and exhausted chunks. It is written through a temporary file and a rename, so a kill
  mid-write leaves the previous version whole.
- Every chunk attempt is appended to `batch_events.jsonl`, with `counts_as_attempt` and the end of the
  driver's stderr, from which account details are removed.

Run a long batch from a checkout that does not change while it runs (a detached worktree at the pushed
commit), and outside any process that may be closed under it. On Windows, `start /affinity FFFF` caps it
at 16 logical CPUs; the Java and MS-FINDER children inherit the mask. Add `/belownormal` to `start` so the
runner itself also runs at below-normal priority.

## What a chunk keeps (one importable bundle)

```
units\<unit_key>\<selection>\cNNN\
  msemblator_run.json   the record (below)
  msemblator.log        every line printed, time-stamped, account details removed
  name_map.tsv          msemblator_id, name, unit_id, alignment_id, structure_step
  formula_summary.csv  formula_score.csv  structure_score.csv  structure_summary.csv
  metfrag\              MetFrag's outputs, copied out of the checkout, which the next run clears
  work\                 MS-Emblator's work folder, moved out of the checkout: both SIRIUS projects,
                        the SIRIUS formula, structure and CANOPUS summaries, the MS-FINDER outputs,
                        save_folder\formula_fixed_msp\id_change.msp
```

`msemblator_run.json` holds `msemblator_commit`, `sirius_version`, `mode`, `input`, `input_sha256` (the
chunk's sha256), `input_records`, `manifest`, `parameter_file_sha256`, `started_at`, `finished_at`,
`duration_seconds`, `exit_code`, `error`, `login_ok`, `failed_steps` (for example
`["structure:MetFrag"]`; `[]` when clean), `stage_seconds` (`formula.sirius`, `formula.msfinder`,
`formula.msbuddy`, `formula.summary`, `formula.total`, `structure.sirius`, `structure.metfrag`,
`structure.msfinder`, `structure.summary`), `rankers` (per step, the ranking model MS-Emblator loaded:
`file`, `ranker` as `xgboost`, `catboost` or `lgbm`, and `sha256`; a score's identity is
`msemblator.<step>.rank_score.<ranker>`), `service_error_lines`, `name_map_ok`, `name_map`
(`records`, `sent_to_structure`, `skipped_no_formula`, `problems`), `metfrag_files`, `harvest`
(`{relative path: {sha256, bytes}}`), `canopus_summary_files`, `outputs` (`{csv: sha256}`) and
`environment` (MS-Emblator commit and dirty flag, model sha256, library DB size and mtime, Python, affinity
mask, logical CPUs, `model_assets` as given with `--model-assets`, for example the Zenodo record and md5
of the installed model set, and `third_party_data`: the data files msbuddy fetches into its package folder
on first use, with size, sha256, source URL and fetch time).

**Models must match the code.** MS-Emblator's rankers check feature names: a model set from an older
release fails on newer code with a CatBoost "should be feature with name ..." error. Install the asset
set the README of the checked-out commit names, and pass its identity with `--model-assets`.

**First use downloads data.** msbuddy fetches about 0.86 GB from its GitHub releases into
`site-packages/msbuddy/data` the first time it runs. Run once on a small input before a batch, so the
fetch is a known event rather than a side effect of the first chunk.

`name_map.tsv` reads only `name`, `unit_id` and `alignment_id` from the manifest, never a `reference_*`
column. `structure_step` is `sent`, or `skipped_no_formula` when MS-Emblator chose no formula and so left
the spectrum out of the structure step (and of CANOPUS).

| Driver exit | Meaning |
|---|---|
| 0 | clean |
| 1 | MS-Emblator reported failed steps (`failed_steps`); the CSVs exist |
| 2 | exception, no stored login, or a name map that is not 1:1 over the input |
| 3 | refused to start (busy `work\` or non-empty output); nothing written, no record |
| 4 | the CSVs exist, but no CANOPUS summary although the SIRIUS structure step did not fail |

## Status and throughput

```
python -m msdial_spectrum_catalog.msemblator_batch status --plan <run>\plan.json
python -m msdial_spectrum_catalog.msemblator_batch throughput --plan <run>\plan.json --target devset=<spectra>
```

`status` counts chunks done, failed and pending, with failure reasons, the attempts each failed chunk has
used (`attempts`, all bundles; `attempts_counted`, those that count against `--max-attempts`), the
pending chunks that hold a partial bundle from an interrupted attempt (`partial_bundles`), seconds and
bytes per spectrum, and any overlapping chunk intervals (there must be none). Take final counts from
`status`, which reads every chunk record.

`backfill-rankers --plan <run>\plan.json --msemblator-root <msemblator> --expected-sha256 <list>` adds
`rankers` to records written before the driver recorded it. It writes nothing for a record whose log names
no ranking model, and lists those chunks under `no_model_in_log`. `throughput` writes `throughput.json`: seconds
per spectrum per stage, bytes per spectrum (`work`, SIRIUS projects, whole bundle), the proposed chunk size
and the extrapolation to each target.
