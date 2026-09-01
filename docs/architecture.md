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
m/z, and intensity remain recorded without inventing a source peak.

## Scale path

SQLite stores searchable metadata and compressed spectrum payloads in v0.1.
When corpus size warrants it, `spectrum_blob` can move to repository/run shards
or Parquet without changing any logical identifier or provenance table.
