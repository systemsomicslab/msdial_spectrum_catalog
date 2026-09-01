# Identifier Policy

## Internal identifiers

Deterministic `urn:msdial:*` identifiers are the catalog primary keys. They cover
studies, analysis units, immutable processing runs, samples, sample features,
deconvoluted spectra, alignment features, alignment members, and consensus spectra.
A run identifier includes a fingerprint of output and input text artifacts, so a
different method or result set cannot silently reuse the same derived spectrum ID.

## Universal Spectrum Identifiers

The HUPO-PSI USI syntax identifies a spectrum through a collection, MS run, index
type, and index. See the [official USI specification](https://github.com/HUPO-PSI/usi)
and the [USI publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC8405201/).

USI is an optional external identifier here. It is populated only when all of the
following are known:

- a repository collection identifier accepted by the target resolver;
- the deposited MS run name;
- an unambiguous native spectrum identifier, scan, or stable index.

An MS-DIAL `MS2SCAN` value is retained as a raw spectrum index but is not assumed
to be a vendor-native identifier. WIFF and Waters data can require multidimensional
native identifiers, and an MS-DIAL deconvoluted or alignment consensus spectrum is
a derived object rather than a single deposited raw scan. Those records therefore
keep internal IDs plus explicit parent links instead of fabricated USIs.

