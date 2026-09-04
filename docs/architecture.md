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

## Scale path

SQLite stores searchable metadata and compressed spectrum payloads in v0.1.
When corpus size warrants it, `spectrum_blob` can move to repository/run shards
or Parquet without changing any logical identifier or provenance table.
