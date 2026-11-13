# Repository rules

Read these rules before editing. Sources are the annotation roadmap charter
(2026-11-13), WS1 catalog foundation specification, and repository contracts.

## Core and schema

- Keep `dependencies = []`. The core uses only the standard library and Python
  3.10 syntax. Never import numpy, pandas, scipy, RDKit or matchms in src/ (D5).
- WS1 owns migrations v11-v13 (C1). Later migrations take the next free version
  after the previous PR merges; never reserve versions in parallel. WS2 declares
  none. Bump SCHEMA_VERSION and test upgrades from the previous main version.
- `class_prediction` is the only per-level class model. WS6 adds exactly
  `tool_subject`, `tool_candidate`, `canopus_prediction_context` and `tool_import`;
  the last is import bookkeeping, not a pre-assertion store.
- Pre-assertion results use class tables, WS6's subject/candidate/context stores
  or WS7's overlay tables (C9). assertion_id is NOT NULL on assertion evidence
  and candidate tables. Only WS8 converts these results into assertions.

## Evidence contracts

- Cosine-family kernels store cosines, never cos² (C4). Other kernels store the
  convention their identity names. Always set score_convention. MS-DIAL exported
  dot products are already cosines; explicitly squared columns retain that name.
- `spectrum_similarity.method` holds the full C3 kernel identity in WS2's
  canonical `<base>?<sorted k=v&...>` form, including preprocessing parameters.
  Preserve existing legacy rows and the similarity UNIQUE key.
- Ingest and adapters never derive levels. Levels are explicit inputs until
  WS8's versioned derivation. SP/SC/SI/CP describe specificity. WS8 derives the
  most specific claim true of every evidence-filtered candidate (charter §6).
  Never store an evidence-completeness grade; ML filters derive it from inventory.
- Vocabulary concepts join by `(vocab_version, concept_id)`.
- Results join by unit + alignment ID (PEAKID in the MSP COMMENT), and through
  `.mdpeakid.tsv` by Peak ID (C5). Tools rewriting NAME must keep a join manifest.

## Private data

- VS21 is identified by name and sha256 only: `MSMS-Private-pos-VS21.msp`
  `0fa0a96b…4804f`, `MSMS-Private-neg-VS21.msp` `7e2581ab…dc4d132f`.
  GLACIER-VS15 and campaign `Msp20261013094720_NCDK-TUAT-LC25_no-unknown-db.lbm2`
  follow the same rules; the public Zenodo LBM2 is a different file.
- Never publish private-library content, hit names, InChIKeys, SMILES, library
  records or library paths in commits, fixtures, PRs, issues or shared logs/reports.
  Never send private reference spectra to SIRIUS. Private catalogs are local only
  and end in `.local-vs21.sqlite`.
- Only `linkage.read_annotator_settings` opens `msp_annotator_settings.tsv`.
  It takes annotator_id, evidence_tier and library_kind by column name only.
  Never read, store, log or hash msp_file_path. Never register the settings file
  as an artifact or run_metadata_source. Share this reader with the export when
  both implementations reach main.
- Only whitelisted manifest fields enter the database; refuse path-like values.
  Never parse workflow-settings.json, method.txt or
  MS_DIAL_Materials_and_Methods.txt into the database.
- D3 permits local reads of F: unit outputs, private libraries and annotation
  outputs and local builds. Publishing remains restricted. PRs carry counts,
  annotation output file names and checksums only; fixtures are synthetic/public.
  Under `D:\13_MSDIAL_Public_Reanalysis\context` read only the charter and your
  task's own spec, never context\prompts or other context files.
- Never commit databases, generated JSON, MS-DIAL outputs, libraries or ICD files.

## Operations

- One worktree per agent (C10); never switch or modify another checkout. Add
  detached build pins only from your worktree, under
  `D:\0_SourceCode\msdial_spectrum_catalog-pins\<sha>` (WS0 dispatch override).
  Build from committed, pushed pins and record commit and pin path in run.json.
- C7 outputs live at `annotation\<ws-short-name>\<run-id>\` with run.json
  recording input checksums, code commit, environment and parameters. WS1 writes
  only annotation\devset, annotation\catalog and annotation\_handbacks on F:.
- Never write under analysis/. Open ledgers only with `file:...?mode=ro`, SELECT
  only. Refuse nonempty WAL without --snapshot; snapshot never uses immutable.
  Never touch ports 8765/8766/8790 or campaign processes. Sequential ingests,
  single thread, at most 16 threads for any job.
- Test exactly in your worktree (never install into the interpreter):
  `$env:PYTHONPATH='src'; & 'C:\Users\Hiroshi Tsugawa\AppData\Local\Python\pythoncore-3.14-64\python.exe' -B -m unittest discover -s tests`
- Every catalog command uses `$env:PYTHONPATH='<worktree-or-pin>\src'; & 'C:\Users\Hiroshi Tsugawa\AppData\Local\Python\pythoncore-3.14-64\python.exe' -B -m msdial_spectrum_catalog.cli <subcommand> ...`.
  Never install or invoke the console script. Keep the CLI __main__ guard.
- Create English commits and PRs. Never merge PRs, force-push, or rebase a pushed
  branch; merge main to resolve integration conflicts.
