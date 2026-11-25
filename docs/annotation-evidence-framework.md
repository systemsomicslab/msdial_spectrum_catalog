# Annotation evidence framework

**Status: working proposal, version 0.1, 25 November 2026.** This page introduces an idea for
discussion. It is not a standard and no standards body has endorsed it. A full description with an
evaluation is in preparation for publication; until then terms may change, and each change is listed
under [Version history](#version-history).

## Purpose

A confidence level alone does not say what an annotation proposes or what supports it. Two
annotations reported as "Level 3" can propose one structure, a lipid species or only a compound class,
and can rest on a library match or on an in-silico prediction. The framework keeps three things apart,
so that each can be read, searched and re-evaluated on its own:

1. **Claim**: what is proposed and how specifically, given as a level, with Level 3 qualified by a
   claim tag.
2. **Evidence inventory**: which kinds of evidence support the claim.
3. **Study-level criteria**: the thresholds each kind of evidence was judged against, stated once per
   study.

The vocabulary grew out of a project-level proposal for Shin-MassBank annotation metadata. Its
machine-readable form is the vocabulary version `smb-v2-consensus` in this repository
([registry file](../src/msdial_spectrum_catalog/vocabulary/registry/smb-v2-consensus.json)), and
[annotation_evidence_model.md](annotation_evidence_model.md) describes how the catalog stores it.

## Claim

Levels are those of Schymanski et al. (2014). Only Level 3 takes a claim tag.

| Notation | Claim | Example |
| --- | --- | --- |
| `L1` | Confirmed structure, by an authentic reference standard measured under comparable conditions | |
| `L2` | Probable structure: one structure is supported, for example by an unambiguous library match or diagnostic evidence, without a reference standard | |
| `L3-SP` | Structure(s) proposed: one two-dimensional constitutional structure, or a short closed list of them, not confirmed by a reference standard. Stereochemistry may remain unresolved unless it is claimed | adenosine 3′- or 5′-monophosphate, whose spectra cannot be told apart |
| `L3-SC` | Substructure complete: the components, building blocks or modifications are proposed, but their positions, connectivity or arrangement are not resolved | PC 16:0_18:1; a known compound plus one oxygen at an unknown position |
| `L3-SI` | Substructure incomplete: a diagnostic substructure or motif is proposed, without a component-level description | a glucuronide recognised by its 176.032 Da neutral loss, aglycone unresolved; PC 34:1 |
| `L3-CP` | Class proposed: a chemical or metabolite class, with no more specific claim | a flavonoid, named in a stated ontology such as ChemOnt (ClassyFire) |
| `L4` | Unequivocal molecular formula | |
| `L5` | Exact mass of interest | |

How the claim is chosen:

- **Specificity, not confidence.** `SP` is more specific than `SC`, which is more specific than `SI`,
  which is more specific than `CP`. How well a claim is supported is said by the evidence inventory and
  the criteria, never by the claim tag.
- **The claim is the most specific statement true of every remaining candidate.** The candidates are
  the structures left after the evidence has excluded the rest. If the evidence cannot tell structure A
  from structure B, the claim is `L3-SP` with the list "A or B"; if A and B share only their components,
  it is `L3-SC`.
- **A list must be closed.** Candidates outside it are excluded by a stated, measured criterion, such as
  a score margin or a decoy-based false discovery rate. The unexamined top hits of a search are not a
  closed list.
- **Context may narrow the list only when it is recorded.** Biological, biosynthetic or dietary context
  may remove candidates, for example compounds the studied organism cannot make, but every removal is
  recorded and `CO` is added to the inventory, so the claim without context can always be recomputed.
- **Combinations.** `CP` may be appended to `SC` or `SI` when the class adds information
  (`L3-SC+CP`). `SP+CP` is reserved for exceptions, and `SC+SI` is avoided.
- **The level is not computed from the evidence tags.** The same evidence can support `L2` or `L3-SP`,
  depending on whether the alternatives are excluded.

## Evidence inventory

| Tag | Evidence | Includes |
| --- | --- | --- |
| `RS` | Reference standard | Authentic standard measured under comparable conditions |
| `FM` | Formula/mass | Accurate mass, isotope pattern, adduct consistency, formula plausibility |
| `SL` | Spectral library | MS/MS library or literature spectrum match |
| `DF` | Diagnostic fragment | Diagnostic product ions, neutral losses or fragmentation rules |
| `RT` | Retention information | Retention time or index match, or retention prediction |
| `IM` | Ion mobility | CCS, drift time, mobility separation; experimental or predicted values |
| `IS` | In silico | SIRIUS, MS-FINDER, MetFrag, ClassyFire, ZBRA, predicted spectra or classes, cheminformatics |
| `MN` | Molecular network | Network edge, known-neighbour relation or interpretable mass shift |
| `HO` | Homologue evidence | Kendrick mass defect, homologous series |
| `CO` | Contextual | Biological, biosynthetic, taxonomic, dietary, microbiome, tissue or perturbation context |
| `UN` | Unclassified additional evidence | Evidence that fits no other tag, such as two-dimensional chromatography |
| `OS` | Other spectroscopy | UV, IR, ECD |

- A tag is listed when evidence of that kind met the study's criteria and supports the claim. Evidence
  that was examined and failed is not tagged, although it may be recorded. A claim is never promoted
  because it carries many tags.
- A match against a spectrum predicted from a structure is `IS`, not `SL`.
- Tags are not independent. In-silico tools are trained on public libraries, and a network neighbour
  can be the library reference itself, so a record should name the tool, library and version behind
  each tag. A library is cited by name, version and checksum.
- A separate tag for similarity in a learned spectrum representation (`SM`) is under discussion and is
  not part of this version.

### Notation

```text
L3-SP[FM,SL,RT]
L3-SC+CP[FM,DF,MN]
L4[FM]
```

`level ['-' claim {'+' claim}] '[' [tag {',' tag}] ']'`. A claim tag is allowed only at Level 3. Tags
are separated by commas and written in the order of the table above. `|` is never used, because
mzTab-M reads it as a multi-value separator. Alternatives such as `L3-SC/SI` are not written: report
the claim that is true of all candidates, and list the candidates of an "A or B" in a separate ranked
field.

## Study-level criteria (recommended)

We recommend that each study state, once, the criteria each evidence tag was judged against, and
publish them in a machine-readable form with an identifier and a version that the annotation records
refer to. Doing so makes the annotations reproducible: a reader can see why a tag was given, and can
re-tag the evidence under other thresholds without repeating the annotation. This is a recommendation,
not a condition for using the notation.

Examples of the form a criterion takes. The values are illustrations, not recommended thresholds;
suitable values depend on the instrument, the method and the library.

| Tag | Example criterion |
| --- | --- |
| `FM` | MS1 mass error < 1 ppm; isotope pattern and adduct consistent; formula passes plausibility rules such as the Seven Golden Rules |
| `SL` | Cosine (dot product) > 0.70 and reverse cosine > 0.80, with matching precursor and adduct type; the score convention is stated |
| `DF` | A predefined diagnostic ion, neutral loss or fragmentation pattern, for example the PC head-group ion m/z 184.0733 in positive mode |
| `RT` | Against a standard: ΔRT ≤ 0.01 min; against a predicted RT: ΔRT ≤ 0.1 min or a platform-specific tolerance |
| `IM` | CCS error ≤ 5 Å², or relative error within a platform-specific tolerance |
| `IS` | A tool-specific score, rank, posterior probability or FDR threshold, for example a ZBRA q-like value < 0.01 |
| `MN` | Edge score > 0.80 and an interpretable mass shift to a known or better-supported node |
| `CO` | An explicit curation rule, for example taxon-specific biosynthesis or a dietary or microbial origin |

A machine-readable criteria set can be as small as this:

```text
criteria_set_id  version  tag  measure             comparator  threshold  unit  score_convention
example-study    1        SL   cosine              >           0.70             cosine
example-study    1        SL   reverse_cosine      >           0.80             cosine
example-study    1        RT   abs_delta_rt_pred   <=          0.1        min
```

Each annotation then cites `criteria_set_id`. Measured values such as Δm/z, ΔRT or a score are
optional per annotation, and are most useful when stored beside the threshold they were compared
against.

## MS-DIAL annotation code ("Annotation tag (VS1.0)")

MS-DIAL 4 and 5 write a three-digit code with each annotation. The code has not been described in a
peer-reviewed publication, so it is documented here as implemented in MS-DIAL 5
([`GetAnnotationCode`](https://github.com/systemsomicslab/MsdialWorkbench/blob/92cbf0d3867948d0668862cfd9ba0756ef6f871d/src/MSDIAL5/MsdialCore/Utility/DataAccess.cs#L1319-L1464),
commit `92cbf0d`). The code records which match verdicts MS-DIAL reached; it is not a claim in the
sense above. From this version on, the claim and the evidence inventory supersede it as the way to
report annotation confidence, and the code stays in the exports for compatibility.

The first digit says what supports the match besides m/z; the second says what the spectrum
established.

| Code | Verdicts reached | Can contribute to the inventory as |
| --- | --- | --- |
| `4x0` | MS/MS | `SL` (x = 3) or `DF` (x = 0–2) |
| `3x0` | RT or RI, and MS/MS | as `4x0`, plus `RT` |
| `2x0` | CCS and MS/MS (ion mobility data) | as `4x0`, plus `IM` |
| `1x0` | RT, CCS and MS/MS (ion mobility data) | as `4x0`, plus `RT` and `IM` |
| `450` | MS/MS matched to an in-silico predicted spectrum | `IS`, never `SL` |
| `440` | EI-MS (GC-MS) | `SL` |
| `340` | RT or RI, and EI-MS (GC-MS) | `SL`, `RT` |
| `500` | RT, CCS and m/z, no MS/MS match (ion mobility data) | `FM`, `RT`, `IM` |
| `510` | CCS and m/z, no MS/MS match (ion mobility data) | `FM`, `IM` |
| `520` | RT and m/z, no MS/MS match | `FM`, `RT` |
| `530` | m/z only | `FM` |
| `999` | No annotation | none |

x = 0: lipid acyl positions resolved; 1: acyl chains resolved; 2: lipid class resolved; 3: MS/MS
matched to an acquired reference spectrum. Every LC-MS code other than `999` also rests on a precursor
m/z within the search tolerance, which can contribute `FM`.

- Codes with x = 0–2 come from MS-DIAL's rule-based lipid annotation, which evaluates diagnostic
  product ions and neutral losses, so they are diagnostic-fragment evidence rather than library
  similarity.
- `450` is never promoted by RT or CCS and never takes a lipid level. Peptide matches keep `430`, and
  GC-MS matches to predicted EI spectra keep `440` or `340`.
- The RT or RI verdict requires agreement within the user's tolerance, capped at 2 min (150 Kovats or
  50,000 Fiehn retention-index units), and gives no verdict when the reference has no retention value.
  The CCS verdict uses the user's tolerance as it is.
- `530` covers three different findings: no MS/MS spectrum, an MS/MS comparison below threshold, and an
  MS/MS comparison that explained nothing. The `Evidence source` column (`PrecursorOnly`,
  `WeakSpectrumMatch`, `UnmatchedSpectrum`) tells them apart.
- The code does not encode isotope agreement, lipid double-bond positions, manual acceptance, or
  whether the reference is an authentic standard or a public library spectrum: an in-house standard
  matched with RT is `330`, like any other. A structure accepted from MS-FINDER is recorded only in the
  `Evidence source` column (`InSilico`).
- On the GC-MS branch, `340` is assigned on retention agreement without re-checking the EI-MS verdict,
  so it should be read together with the spectrum score columns.
- The code is written to the `Annotation tag (VS1.0)` column of the peak (`.mdpeak`) and alignment
  (`.mdalign`) text exports; to the mzTab-M `reliability` column, declared as
  `[MS, MS:1003032, compound identification confidence code in MS-DIAL, ]`; and to
  `annotation_tag_vs1` in the optional per-candidate export `.mdcandidate.tsv`. It is not shown in the
  graphical interface and is not written to MassBank or MSP exports.
- Earlier versions differ. `450` was added in October 2026, and in September 2026 the RT verdict was
  changed so that a reference without a retention time no longer yields a `3x0` code.

## Support in MS-DIAL

The aim is for MS-DIAL to write the claim and the evidence inventory in its standard outputs, so that
the inventory is visible wherever MS-DIAL results are used. MS-DIAL can supply only part of it:
context, other spectroscopy and many in-silico results come from elsewhere. MS-DIAL does not yet write
the claim or the inventory itself. The state of the MS-DIAL 5 source in November 2026:

| Tag | Computed and exported today | Not yet available |
| --- | --- | --- |
| `RS` | none | A declaration that a library holds authentic standards measured on the same method |
| `FM` | Precursor m/z agreement, adduct type, MS1 isotope pattern, the reference's formula | Isotope similarity in the alignment table and mzTab-M; formula assignment, which MS-FINDER provides |
| `SL` | MS/MS similarity scores (cosine, not squared), matched peaks, `Evidence source`; up to three ranked candidates per annotator in `.mdcandidate.tsv` and mzTab-M | Library name, version and checksum per row; reference collision energy |
| `DF` | Rule-based lipid annotation (x = 0–2 above); diagnostic-rule results written by MS-FINDER | Which rule, ion or neutral loss supported a lipid annotation |
| `RT` | RT or RI verdict, reference RT, RT similarity; a retention-order check for lipid annotations | Whether the reference RT was measured on the same method or predicted |
| `IM` | CCS, reference CCS, CCS verdict and similarity | Measured versus predicted reference CCS |
| `IS` | Matches to predicted spectra (`450`); MS-FINDER structures accepted in the graphical interface (`Evidence source`) | Tool version, rank and score type per candidate; in-silico tools run from the console |
| `MN` | Molecular networks from the graphical interface; export for GNPS feature-based molecular networking | Network evidence linked to individual annotations |
| `HO` | none | Homologous-series evidence in the exports |
| `CO`, `UN`, `OS` | none | Outside MS-DIAL today |

## How to cite

- For the Level 3 sublevels `L3-SP` and `L3-CP`, cite Matsuda et al. (2026), which reported them as
  Level 3-SP (tentative annotation by planar structure) and Level 3-CP (tentative annotation by compound
  class) and applied them to ZBRA-based annotations with statistical (q-value) criteria. Its "planar
  structure" corresponds to the two-dimensional constitutional structure here.
- For the framework as a whole (the three layers, `L3-SC` and `L3-SI`, the evidence inventory and the
  study-level criteria) and for MS-DIAL's annotation code, cite this page by its URL at a fixed commit
  or tag, with the version number:

> systemsomicslab. Annotation evidence framework: claim, evidence inventory and study-level criteria.
> Working proposal, version 0.1 (25 November 2026).
> https://github.com/systemsomicslab/msdial_spectrum_catalog/blob/COMMIT/docs/annotation-evidence-framework.md

## Relation to other frameworks

- **Schymanski et al. (2014).** The levels are theirs. Their distinction between Level 2a (library) and
  2b (diagnostic) is carried here by the `SL` and `DF` tags rather than by sublevels, and the claim tags
  subdivide their Level 3 by how specific the proposal is.
- **Metabolomics Standards Initiative (Sumner et al. 2007).** Its levels separate identified compounds,
  putatively annotated compounds, putatively characterised compound classes and unknowns. `L3-CP` is
  close to its class level.
- **Evidence codes attached to a level.** Other proposals add letters or codes for the kinds of evidence
  used to a confidence level, and address the same need. The arrangement here keeps the inventory apart
  from both the level and the claim, and states the thresholds once per study, so that evidence can be
  re-tagged when criteria change.
- Lipid examples use the shorthand notation of Liebisch et al. (2020).

## References

1. Matsuda F, Hayasaka R, Torigoe T, Takahashi M, Takahashi Y, Yoshizawa AC, Tsugawa H, Okuda S,
   Neumann S, Schymanski EL, Izumi Y, Hirayama A. MassBank Human: The annotated human metabolome.
   bioRxiv 2026.10.10.758200 (2026). <https://www.biorxiv.org/content/10.64898/2026.10.10.758200v1>
2. Schymanski EL, Jeon J, Gulde R, Fenner K, Ruff M, Singer HP, Hollender J. Identifying small molecules
   via high resolution mass spectrometry: communicating confidence. Environ Sci Technol 48, 2097–2098
   (2014). <https://doi.org/10.1021/es5002105>
3. Sumner LW, Amberg A, Barrett D, et al. Proposed minimum reporting standards for chemical analysis.
   Metabolomics 3, 211–221 (2007). <https://doi.org/10.1007/s11306-007-0082-2>
4. Liebisch G, Fahy E, Aoki J, et al. Update on LIPID MAPS classification, nomenclature, and shorthand
   notation for MS-derived lipid structures. J Lipid Res 61, 1539–1555 (2020).
   <https://doi.org/10.1194/jlr.S120001025>

## Version history

- 0.1 (25 November 2026): first public summary.
