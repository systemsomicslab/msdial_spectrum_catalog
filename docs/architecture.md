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
records; alignment features link to `.mdalign` rows; alignment members link to
`.mdprovenance.tsv` rows; consensus spectra link to alignment `.mdmsp` records.

For gap-filled cells, `has_source_peak=false` and `feature_id` is null. The
corrected RT, m/z, and intensity remain recorded as an alignment member without
inventing a source peak.

## Scale path

SQLite stores searchable metadata and compressed spectrum payloads in v0.1.
When corpus size warrants it, `spectrum_blob` can move to repository/run shards
or Parquet without changing any logical identifier or provenance table.

