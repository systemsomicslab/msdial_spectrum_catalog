# Class ontology and structure lookup

Migration v13 defines the only class-prediction model. Terms belong to an explicit
ontology/version (`classyfire.chemont` or `npclassifier`; `unknown` is explicit).
Terms use native IDs when available and casefolded names otherwise. Parent edges
permit multiple superclasses for an NPClassifier class.

`ClassTermInput` and `upsert_lineage` accept levels in increasing order, never
mixing ontologies or versions in a lineage. `ClassPredictionInput` carries a
subject, lineage, per-level posteriors, rank and optional formula context/rank and
source-row reference. `record_class_predictions` requires an existing tool run
and spectrum/alignment subject, validates probabilities, stores reported values
without renormalising, and marks the deepest level as most specific. These are
pre-assertion records: they do not derive levels or claims. WS6 maps CANOPUS here.

Use `-m msdial_spectrum_catalog.cli show-classes DB SUBJECT_ID` to inspect them.
All commands use the interpreter/PYTHONPATH invocation in AGENTS.md.

`load-structure-classes --database DB --icd PATH` reads MS-FINDER's tab-separated
`InChIKey`, `Direct parent name`, `Direct parent ID` layout. The library identity
is the file name and SHA-256; its path is never stored or reported. Each valid
lookup assigns an InChIKey to its native ChemOnt term under one
`msfinder.icd_lookup` tool run. Repeating the same file adds no rows. If run.json
exists beside the database, its counts, database checksum and lookup identity
are updated as part of the build.

**Observed layout limitation:** this three-column file has neither ontology
version, parent depth nor hierarchy edges. Lookup terms use version `unknown`
and level `direct_parent`; no kingdom/class depth or parent edges are invented.
`direct_parent` is a lookup-only level, outside the ordered prediction API.
An emitted, versioned lineage is required to provide actual hierarchy edges for
lowest-common-ancestor inference. This is an explicit deviation from the WS1
specification's expectation that the ICD would supply those edges.
An emitted lineage with the same ontology version and native term ID can replace
the lookup-only level with its known depth, preserving all existing assignments;
reloading the lookup never downgrades that depth.

No ontology dictionary or ICD data is bundled. Tests create an invented lookup
with the same layout at runtime.
