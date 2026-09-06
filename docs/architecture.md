# Architecture

## Boundary with Repository Catalog

`msdial_repository_catalog` answers which studies and analysis units exist and
why they are relevant. This project answers exactly which analytical artifact,
sample feature, and spectrum supports a downstream claim. The stable bridge is
`repository + accession + external analysis_unit_id`.

The databases are separate because metadata crawling and spectrum reanalysis have
different volume, update, retention, and release policies.

## Provenance graph

An MS-DIAL run is immutable and identified by a hash of its registered text
artifacts. Re-running with different parameters creates another run. Sample
features link to `.mdpeak` rows; deconvoluted spectra link to individual `.mdmsp`
records; alignment features link to `.mdalign` rows; alignment members primarily
link through the compact `.mdpeakid.tsv` matrix; consensus spectra link to
alignment `.mdmsp` records. Optional `.mdprovenance.tsv` rows add detailed audit
properties without changing the compact source-feature link.

For gap-filled cells, the matrix stores `-1`, `has_source_peak=false`, and
`feature_id` is null. When the detailed audit sidecar is present, corrected RT,
m/z, and intensity remain recorded without inventing a source peak, and
`peak_origin` names why there is no source peak: `detected`, `gap_filled` or
`absent`. The sidecar used to distinguish those by sentinel instead, writing `-2`
for a gap-filled peak id where the matrix wrote `-1`, so the two artifacts
disagreed about the same cell; naming the state removed the disagreement.

Ingestion refuses to read a sentinel as a measurement, and validation refuses to
pass one. A member with a source peak must carry a usable m/z, and a member
without one must carry no raw-spectrum index. Both checks exist because a real
sidecar once stored the inverse of both — every member that had a source peak
reported `mz = -1`, and every gap-filled member reported MS1 and MS2 scan 0,
which is a real scan index pointing at an unrelated spectrum — through a run that
reported valid with no errors and no warnings. The exporter side is fixed in
MsdialWorkbench; the guards here also cover sidecars written before that fix, and
`peak_origin` is read by name so an older file without it still ingests.

## Annotation candidates against the representative

`msdial_annotation_result` records what MS-DIAL's text export said about a feature: one row, the
representative. `msdial_annotation_candidate` records what the search kept: every threshold-passing
result, ranked. The two tables describe one decision from different sides and are checked against each
other rather than merged, because they can disagree only if one of them describes a different run.

Four invariants hold on every candidate set, and each is an error rather than a warning:

- the ranks of one subject are `1..candidate_count`, once each, and every row of the set agrees on the
  count
- exactly one row is the representative, and it is rank 1
- a spectral score exists only where `is_spectrum_comparison_performed` is true
- the rank 1 name equals the `candidate_name` of the `.mdalign` row for the same feature, after
  MS-DIAL's `no MS2: ` and `low score: ` prefixes have been stripped

The last one is the cross-artifact check. It held on all 1919 annotated features of the reference demo
before it was made an error.

No `annotation_kind` is derived for a candidate. The three-way kind comes from a name prefix that only
the representative carries, and the sidecar states `is_reference_matched`, `is_annotation_suggested` and
`is_spectrum_comparison_performed` directly, which is more information than the prefix. Deriving a kind
from those would also mislabel a text-database result, which has no reference spectrum to compare rather
than no product-ion spectrum to compare with.

## Scale path

SQLite stores searchable metadata and compressed spectrum payloads in v0.1.
When corpus size warrants it, `spectrum_blob` can move to repository/run shards
or Parquet without changing any logical identifier or provenance table.
