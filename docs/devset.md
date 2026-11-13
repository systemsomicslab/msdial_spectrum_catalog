# Development-set selection

Run `-B -m msdial_spectrum_catalog.cli devset build --campaigns-root ROOT
--campaign ID [--campaign ID ...] --exclude UNIT=REASON --out FILE` with the
absolute interpreter and PYTHONPATH documented in AGENTS.md. Builds use a
committed, pushed detached pin. The JSON records that commit.

The generator enumerates the ledger's unit table and cross-checks its state with
the latest transition. Campaign manifests are not consulted, so split parts are
included. Only done units are selected. Acquisition types determine DIA flags;
the ledger determines the Both ion-mode flag. Alignment MSP headers supply total
and nonempty record counts; sample count comes from mdpeak files.

Ledgers open with `mode=ro`, using SELECT statements in a single read transaction.
A nonempty WAL requires `--snapshot`; this mode never uses immutable access and
records ledger/WAL sizes and nanosecond mtimes. Changes during selection fail.
Before/after observations are recorded alongside the maximum transition sequence
and its timestamp, read in the same transaction as the units. A done unit without
a transition fails.
Without snapshot, immutable fallback is allowed only with an empty/absent WAL.
No ledger or campaign file is written.

Use `--snapshot --name stage3b-snapshot` for the live campaign. Snapshot requires
an explicit name and always has
no smoke/small subsets. The standard selection contains the C6 subsets restricted
to selected units. Every consumer uses the generated file rather than selecting
again. Paths remain local; only counts and checksums belong in public reports.

The content digest covers canonical, sorted-key, compact UTF-8 JSON of campaigns,
exclusions, subsets (when present), and units. Creation time and generator fields
are excluded. Existing outputs are refused unless `--force-new-date` is explicitly
provided. Ledger/WAL mtimes are also excluded from the digest because runner
heartbeats can change them without changing the selection. Unknown exclusions fail.
Never write under the analysis root.
