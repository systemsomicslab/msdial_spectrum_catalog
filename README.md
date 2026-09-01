# MS-DIAL Spectrum Catalog

MS-DIAL Spectrum Catalog records how a public repository study becomes a set of
sample features, deconvoluted spectra, alignment features, and consensus spectra.
It is a companion to `msdial_repository_catalog`: repository metadata remains in
that catalog, while high-volume analytical results and spectrum provenance live here.

## Traceability model

```text
repository -> accession -> analysis unit -> MS-DIAL run
                                      |-> sample -> mdpeak feature -> sample mdmsp spectrum
                                      |-> alignment spot -> member feature(s)
                                                        -> alignment mdmsp consensus spectrum
```

Every imported artifact is recorded with SHA-256, byte size, path, and record
location. Spectrum peak arrays are zlib-compressed and deduplicated by content hash.

MS-DIAL Console additionally exports `AlignResult-*.mdprovenance.tsv`. This sidecar
is the authoritative mapping from every alignment spot and file to the original
`MasterPeakID`, raw spectrum indices, and representative feature. Existing
`.mdalign`, `.mdmsp`, `.mdpeak`, and mzTab-M formats remain unchanged.

## Quick start

```powershell
python -m pip install -e .
msdial-spectrum-catalog ingest-run catalog.sqlite D:\analysis\output `
  --repository mb_post `
  --accession MPST000007 `
  --analysis-unit lc_ms_negative_dda `
  --analysis-files D:\analysis\analysis_files.csv `
  --parameter-file D:\analysis\method.txt
```

The command exits with code 2 when an alignment exists without its provenance
sidecar or when an MSP record references a missing feature.

Re-check a stored run at any time:

```powershell
msdial-spectrum-catalog validate-run catalog.sqlite "urn:msdial:run:..."
msdial-spectrum-catalog show-spectrum catalog.sqlite "urn:msdial:spectrum:..."
```

## Identifier policy

The catalog uses deterministic `urn:msdial:*` identifiers as primary keys. A PSI
Universal Spectrum Identifier (USI) is stored separately when the repository,
run name, and an unambiguous native spectrum identifier are available. MS-DIAL
deconvoluted and alignment consensus spectra are derived records, so a guessed
USI is never generated from an array index alone.

## Annotation-ready schema

The schema already separates:

- annotation level (`L1` to `L5`)
- Level 3 claim (`SP`, `CP`, `MO`, `CL`)
- evidence tags (`RS`, `FM`, `SL`, `DF`, `RT`, `IM`, `IS`, `MN`, `CX`)
- versioned criteria sets and measured evidence values

No automatic Level assignment is implemented yet. That policy will be developed
after the provenance foundation has been validated.
