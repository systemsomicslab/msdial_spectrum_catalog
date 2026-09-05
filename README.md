# MS-DIAL Spectrum Catalog

MS-DIAL Spectrum Catalog records how a public repository study becomes a set of
sample features, deconvoluted spectra, alignment features, and consensus spectra.
It is a companion to `msdial_repository_catalog`: repository metadata remains in
that catalog, while high-volume analytical results and spectrum provenance live here.

## Purpose

The programme goal this catalog serves is to **redistribute the world's experimental MS/MS spectra with
their metadata**. A spectrum on its own is of limited use; what makes it reusable is everything attached to
it. Three kinds of metadata are in scope:

- **chemical** — the molecule or candidate set, the annotation claim, and the evidence inventory behind it
- **biological** — organism, tissue, disease, intervention, diet, age, sex and the rest of the study context
- **response** — how the molecule changed, for example an increase under inflammation, tied to an explicit
  contrast

The reason to insist on all three is empirical: mzmine has a large natural-products user base largely
because it reaches GNPS, and what users want from GNPS is metadata on their peaks. Spectra without context
do not get reused.

Licensing follows from the same distinction. A commercial reference library such as NIST cannot be
redistributed, but an *experimental* spectrum annotated with the help of one can be, provided the record
cites the library by identifier, version and checksum rather than carrying its content. The evidence model
is built to make that citation explicit.

Volume is the other constraint. Retaining every raw MS/MS from every study would grow without bound and
without adding information, so spectra are reduced before publication: within a study, a consensus spectrum
per alignment feature across samples; and across a reference library, structures whose spectra cannot be
told apart are reported as a candidate set rather than as separate confident identifications. A learned
representation such as a DreaMS embedding is a complementary reduction, not a substitute — the two answer
different questions and can be carried side by side.

## Traceability model

```text
repository -> accession -> analysis unit -> MS-DIAL run
                                      |-> sample -> mdpeak feature -> sample mdmsp spectrum
                                      |-> alignment spot -> member feature(s)
                                                        -> alignment mdmsp consensus spectrum
```

Every imported artifact is recorded with SHA-256, byte size, path, and record
location. Spectrum peak arrays are zlib-compressed and deduplicated by content hash.

MS-DIAL Console additionally exports `AlignResult-*.mdpeakid.tsv`. This compact
matrix is the authoritative mapping from each alignment `MasterAlignmentID` and
sample to the original mdpeak `MasterPeakID`; `-1` represents a gap-filled or
missing source peak. `AlignResult-*.mdprovenance.tsv` is an optional audit export
that supplements the matrix with raw spectrum indices, RT, m/z, and intensity,
plus a `peak_origin` column naming why a member has no source peak (`detected`,
`gap_filled`, `absent`). Existing `.mdalign`, `.mdmsp`, `.mdpeak`, and mzTab-M
formats remain unchanged. The detailed audit file is disabled by default and can
be requested in an MS-DIAL Console parameter file with
`Detailed alignment provenance: True`.

`AlignResult-*.mdcandidate.tsv` is a second optional export, requested with
`Annotation candidates: True`. MS-DIAL keeps up to three threshold-passing
results per annotator and alignment carries them into the spot, but `.mdalign`
publishes only the winner. The sidecar publishes all of them, one row per spot
and candidate, ranked by the same precedence that picks the representative.

A sentinel is never stored as a measurement. On ingest, a member without a source
peak carries no raw-spectrum index and a negative m/z becomes null; on
validation, a member with a source peak must carry a usable m/z and a member
without one must carry no scan index. Both are errors rather than warnings,
because a false provenance pointer is worse than a missing one. See
`docs/architecture.md` for what made these necessary.

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

The command exits with code 2 when an alignment has neither a compact Peak ID
matrix nor a detailed provenance sidecar, or when an MSP record references a
missing feature.

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

## Level 3 annotation records

An annotation record separates the claim, the evidence inventory that supports it, and the study-level
criteria each evidence tag was judged against. See `docs/annotation_evidence_model.md` for the full model.

- annotation level `L1` to `L5`, always an explicit input and never derived from the evidence tags
- Level 3 claim `SP` structure, `SC` substructure complete, `SI` substructure incomplete, `CP` class
- evidence tags `RS FM SL DF RT IM IS MN HO CO UN OS`
- an ordered candidate list per claim, so an honest "A or B" is representable
- versioned criteria sets with per-tag thresholds, and measured values stored beside the threshold they
  were compared against
- tool-run provenance for every candidate, and first-class reference and in-silico predicted spectra

The controlled vocabulary is versioned data, not code, because the proposal is still under discussion and
one token changed meaning between generations: `CP` means "component proposed" in the earlier draft and
"class proposed" in the agreed consensus. Rows therefore store a concept identifier plus a vocabulary
version, and the two-letter token is produced only at emit time.

```powershell
msdial-spectrum-catalog vocabulary
msdial-spectrum-catalog notation "L3-SC[FM,DF,MN,CO]" --vocabulary smb-v2-consensus
msdial-spectrum-catalog show-annotations catalog.sqlite "urn:msdial:run:..."
msdial-spectrum-catalog validate-annotations catalog.sqlite "urn:msdial:run:..."
```

MS-DIAL's own per-feature annotation block is ingested into `msdial_annotation_result`, with its four
distinct outcomes kept apart: a real MS/MS match, a low-score match, a precursor-only suggestion with no
product-ion spectrum at all, and no candidate. Only a named MS/MS match may support spectral-library
evidence. No automatic Level assignment is implemented; that policy is a curation decision.

`msdial_annotation_candidate` holds what the search actually kept, not only what it published: every
candidate of every alignment feature, ranked. On the reference FastLC demo, of **1919 annotated
alignment features only 335 have a single candidate** — 115 have two and 1469 have three. Four fifths of
the annotations the older catalog stored as one identification were searches that had not chosen. The
representative is always rank 1, `.mdalign` and the sidecar must name the same winner, and a spectral
score is stored only where a spectrum comparison actually happened; all three are validation errors when
they fail, because a reader cannot recover the truth from a stored row that asserts otherwise.

```powershell
msdial-spectrum-catalog show-candidates catalog.sqlite "urn:msdial:run:..."
msdial-spectrum-catalog show-candidates catalog.sqlite "urn:msdial:run:..." --alignment-feature "urn:msdial:alignment:..."
```

## Ambiguity classes

Different structures often produce MS/MS spectra that cannot be told apart, because a product-ion
spectrum carries structural information only indirectly. When a query matches library entry A and entry
B's spectrum is itself indistinguishable from A's, the match did not choose between them, and the honest
report is "A or B". `docs/ambiguity_classes.md` describes how those classes are computed and reported.

Measured on the public negative VS20 release, 44,353 records in 43 seconds on one core: 1,846 distinct
groups of mutually indistinguishable entries, of which 438 are genuine constitutional isomers with
different skeletons and the same formula. Adenosine 3'-monophosphate against adenosine
5'-monophosphate is one of them.

```powershell
msdial-spectrum-catalog ingest-reference-library catalog.sqlite library.msp --library-name NAME
msdial-spectrum-catalog compute-ambiguity catalog.sqlite
msdial-spectrum-catalog show-ambiguity catalog.sqlite "urn:msdial:reference-spectrum:..."
```
