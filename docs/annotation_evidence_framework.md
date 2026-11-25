# Annotation evidence framework

**Status: working proposal, version 0.1, 25 November 2026.** This page introduces an idea for
discussion, at the level of a working proposal. It is not a standard and no standards body has
endorsed it. A paper describing the full framework is still to be written; until then terms may
change, and each change is listed under [Version history](#version-history).

## Purpose

A confidence level alone does not say what an annotation proposes or what supports it. Two
annotations reported as "Level 3" can propose one structure, a lipid species or only a compound class,
and can rest on a library match or on an in-silico prediction. The framework keeps three things apart,
so that each can be read, searched and re-evaluated on its own:

1. **Claim**: what is proposed and how specifically, given as a level, with Level 3 qualified by a
   claim tag.
2. **Evidence inventory**: which kinds of evidence support the claim.
3. **Study-level criteria** (recommended, optional): the thresholds each kind of evidence was judged
   against, which we recommend stating once per study.

<!-- AUTHOR DECISION: credit wording -->
The vocabulary grew out of a project-level proposal for Shin-MassBank annotation metadata. Its
machine-readable form is the vocabulary version `smb-v2-consensus` in this repository
([registry file](../src/msdial_spectrum_catalog/vocabulary/registry/smb-v2-consensus.json)), and
[annotation_evidence_model.md](annotation_evidence_model.md) describes how the catalog stores it.

What this page takes from elsewhere and what it proposes: the levels are those of Schymanski et al.
(2014) [2]; the Level 3 claim tags, the evidence tags and the rules for combining claim tags are those
of the vocabulary version `smb-v2-consensus`, and `L3-SP` and `L3-CP` correspond to sublevels
introduced by Matsuda et al. (2026) [1]. The rules marked *Proposal of this page* under
[Claim](#claim), the recommendations on recording evidence, the criteria format and the mapping from
MS-DIAL's annotation code to the inventory are proposals of this page, left open for the full paper.

## Claim

Levels are those of Schymanski et al. (2014) [2]. In this proposal only Level 3 takes a claim tag.

| Notation | Claim | Example |
| --- | --- | --- |
| `L1` | Confirmed structure, by an authentic reference standard measured under comparable conditions | |
| `L2` | Probable structure: one structure is supported, for example by an unambiguous library match or diagnostic evidence, without a reference standard | |
| `L3` | Tentative candidate (Schymanski et al.); in this proposal it is normally written with one of the claim tags below | |
| `L3-SP` | Structure(s) proposed: one two-dimensional constitutional structure, or a short ranked list of them, not confirmed by a reference standard. Stereochemistry may remain unresolved unless it is claimed. Level 3-SP of Matsuda et al. (2026) [1] | adenosine 3′- or 5′-monophosphate, whose spectra cannot be told apart |
| `L3-SC` | Substructure complete: the components, building blocks or modifications are proposed, but their positions, connectivity or arrangement are not resolved | PC 16:0_18:1; a known compound plus one oxygen at an unknown position |
| `L3-SI` | Substructure incomplete: a diagnostic substructure or motif is proposed, without a component-level description | a glucuronide recognised by its 176.032 Da neutral loss, aglycone unresolved; PC 34:1 |
| `L3-CP` | Class proposed: a chemical or metabolite class, with no more specific claim. Level 3-CP of Matsuda et al. (2026) [1] | a flavonoid, named in a stated ontology such as ChemOnt (ClassyFire) |
| `L4` | Unequivocal molecular formula | |
| `L5` | Exact mass of interest | |

`L3-SP` and `L3-CP` are the Level 3 sublevels 3-SP (tentative annotation by planar structure) and 3-CP
(tentative annotation by compound class) that Matsuda et al. (2026) [1] introduced under the
Schymanski scheme for ZBRA-based annotations with q-value control. This page writes them as
"structure(s) proposed" and "class proposed" and means the same two claims: its two-dimensional
constitutional structure corresponds to the preprint's planar structure, and its chemical or metabolite
class to the preprint's compound class.

How the claim is chosen:

- **Specificity, not confidence.** `SP` is more specific than `SC`, which is more specific than `SI`,
  which is more specific than `CP`. How well a claim is supported is said by the evidence inventory and,
  where a study states them, the criteria, never by the claim tag.
- **The claim should be true of every candidate still in play.** If the evidence cannot tell structure A
  from structure B, report `L3-SP` with the ranked list "A or B"; if A and B share only their components,
  `L3-SC`. *Proposal of this page.*
- **State how a list was bounded.** We recommend that a list say how candidates outside it were set
  aside, for example by a score margin, a rank cut-off or a false discovery rate (q-value). *Proposal of
  this page.*
- **Context may narrow the list when it is recorded.** Biological, biosynthetic or dietary context may
  remove candidates, for example compounds the studied organism cannot make. We propose that each
  removal be recorded and `CO` added to the inventory, so that the claim without context can always be
  recomputed. *Proposal of this page.*
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

- A tag is listed when evidence of that kind met the criteria the study applied and supports the
  claim. Evidence that was examined and failed is not tagged, although it may be recorded. A claim is
  never promoted because it carries many tags.
- A match against a spectrum predicted from a structure is `IS`, not `SL`.
- Tags are not independent. In-silico tools are trained on public libraries, and a network neighbour
  can be the library reference itself, so we recommend that a record name the tool, library and version
  behind each tag, and cite a library by name, version and checksum.
- A separate tag for similarity in a learned spectrum representation (`SM`) is under discussion and is
  not part of this version.

### Notation

```text
L3-SP[FM,SL,RT]
L3-SC+CP[FM,DF,MN]
L4[FM]
```

A label is written as `level ['-' claim {'+' claim}] '[' [tag {',' tag}] ']'`; the repository's
parser also reads a label without the bracket. A claim tag is used only at Level 3. Tags are separated
by commas and written in the order of the table above. `|` is not used, because mzTab-M reads it as a
multi-value separator. Alternatives such as `L3-SC/SI` are not part of the notation: the claim that is
true of all candidates is reported, and the candidates of an "A or B" go in a separate ranked field.

## Study-level criteria (recommended)

We recommend that each study state, once, the criteria each evidence tag was judged against, and
publish them in a machine-readable form with an identifier and a version that the annotation records
can refer to. Doing so makes the annotations easier to reproduce: a reader can see why a tag was given,
and can re-tag the evidence under other thresholds without repeating the annotation. This is a
recommendation, not a condition for using the notation or the evidence inventory.

The table shows the form a criterion can take. The values are illustrations, not recommended
thresholds; suitable values depend on the instrument, the method and the library.

| Tag | Example criterion |
| --- | --- |
| `FM` | MS1 mass error ≤ 5 ppm; isotope pattern and adduct consistent; formula passes plausibility rules such as the Seven Golden Rules |
| `SL` | Cosine (dot product) > 0.70 and reverse cosine > 0.80, with matching precursor and adduct type; the score convention is stated |
| `DF` | A predefined diagnostic ion, neutral loss or fragmentation pattern, for example the PC head-group ion m/z 184.0733 in positive mode |
| `RT` | ΔRT ≤ 0.1 min against a reference measured on the same method; a wider, stated tolerance against a predicted RT |
| `IM` | CCS error ≤ 5 Å², or relative error within a platform-specific tolerance |
| `IS` | A tool-specific score, rank, posterior probability or false discovery rate threshold, for example a q-value threshold for ZBRA-based annotations |
| `MN` | Edge score > 0.80 and an interpretable mass shift to a known or better-supported node |
| `CO` | An explicit curation rule, for example taxon-specific biosynthesis or a dietary or microbial origin |

A machine-readable criteria set can be as small as the table below. The column names follow the
`criteria_set` and `criteria_rule` tables of this repository, and `-` marks an empty cell.

| criteria_set_id | version | evidence_token | metric | comparison | threshold_value | threshold_unit | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| example-study | 1 | FM | mass_error_ppm | <= | 5 | ppm | - |
| example-study | 1 | SL | cosine | > | 0.70 | - | plain cosine |
| example-study | 1 | SL | reverse_cosine | > | 0.80 | - | plain cosine |
| example-study | 1 | RT | abs_delta_rt | <= | 0.1 | min | reference measured on the same method |

Each annotation can then cite the `criteria_set_id`. Measured values such as Δm/z, ΔRT or a score are
optional per annotation, and are most useful when stored beside the threshold they were compared
against.

## MS-DIAL annotation code ("Annotation tag (VS1.0)")

MS-DIAL 4 and 5 write a three-digit code with each annotation. The code has never been published in a
paper, so it is documented here as implemented in MS-DIAL 5
([`GetAnnotationCode`](https://github.com/systemsomicslab/MsdialWorkbench/blob/92cbf0d3867948d0668862cfd9ba0756ef6f871d/src/MSDIAL5/MsdialCore/Utility/DataAccess.cs#L1319-L1464),
commit `92cbf0d`; the function is unchanged since release 5.5.261105). The code records which match
verdicts MS-DIAL reached and is MS-DIAL's way of reporting annotation confidence; it is not a claim in
the sense above. The claim and the evidence inventory are intended to take over this role once MS-DIAL
can write them (see [Support in MS-DIAL](#support-in-ms-dial)); until then the code remains MS-DIAL's
exported confidence value.

The first digit says what supports the match besides m/z (retention, CCS, or neither); the second says
what the spectrum established. Which spectrum tag a match contributes is read from the
`Evidence source` column, not from the code: `ReferenceSpectrum` → `SL`, `RuleBased` → `DF`,
`InSilico` → `IS`.

| Code | Verdicts reached | Can contribute to the inventory as |
| --- | --- | --- |
| `4x0` | MS/MS | spectrum tag from `Evidence source` (x = 0–2 is always `DF`) |
| `3x0` | RT or RI, and MS/MS | as `4x0`, plus `RT` |
| `2x0` | CCS and MS/MS (ion mobility data) | as `4x0`, plus `IM` |
| `1x0` | RT, CCS and MS/MS (ion mobility data) | as `4x0`, plus `RT` and `IM` |
| `450` | MS/MS matched to a spectrum from a library declared as predicted (in silico) | `IS`, never `SL` |
| `440` | EI-MS (GC-MS) | spectrum tag from `Evidence source` (`SL`, or `IS` for a library declared as predicted) |
| `340` | RT or RI, and EI-MS (GC-MS) | as `440`, plus `RT` |
| `500` | RT, CCS and m/z, no MS/MS match (ion mobility data) | `FM`, `RT`, `IM` |
| `510` | CCS and m/z, no MS/MS match (ion mobility data) | `FM`, `IM` |
| `520` | RT and m/z, no MS/MS match | `FM`, `RT` |
| `530` | m/z only | `FM` |
| `999` | No m/z or MS/MS verdict: no annotation, or a structure accepted from MS-FINDER | none, or `IS` from `Evidence source` |

x = 0: lipid acyl positions resolved; 1: acyl chains resolved; 2: lipid class resolved; 3: MS/MS matched
with no lipid level (usually an acquired reference spectrum; peptide matches to predicted b/y ladders,
and lipids whose rules resolved no class, chain or position level, also take 3). Every LC-MS code other
than `999` also rests on a precursor m/z within the search tolerance, which can contribute `FM`.

- Codes with x = 0–2 come from MS-DIAL's rule-based lipid annotation, which evaluates diagnostic
  product ions and neutral losses, so they are diagnostic-fragment evidence rather than library
  similarity.
- `450` is never promoted by RT or CCS and never takes a lipid level. Peptide matches keep `430`
  (`Evidence source` `InSilico`), and GC-MS matches to predicted EI spectra keep `440` or `340`.
- `450` requires the library to be declared as predicted: a predicted-spectrum MSP loaded as such, or
  spectra from MS-DIAL's EAD, OAD and EIEIO lipid generators in a run whose target is not lipidomics. A
  predicted library loaded as an ordinary MSP is coded `430`, with `Evidence source`
  `ReferenceSpectrum`, so neither the code nor the column shows that its spectra were predicted.
- The RT or RI verdict and the CCS verdict are reached only when the run uses retention or CCS,
  respectively, for scoring. The RT or RI verdict requires agreement within the user's tolerance,
  capped at 2 min (150 Kovats or 50,000 Fiehn retention-index units), and gives no verdict when the
  reference has no retention value. The CCS verdict uses the user's tolerance as it is.
- `530` (and likewise `520`, `510`, `500`) covers three different findings: no MS/MS spectrum, an
  MS/MS comparison below threshold, and an MS/MS comparison that explained nothing. The
  `Evidence source` column (`PrecursorOnly`, `WeakSpectrumMatch`, `UnmatchedSpectrum`) tells them
  apart.
- The code does not encode isotope agreement, lipid double-bond positions, manual acceptance, or
  whether the reference is an authentic standard or a public library spectrum: an in-house standard
  matched with RT is `330`, like any other. A structure accepted from MS-FINDER in the graphical
  interface is coded `999`; only `Evidence source` (`InSilico`) and, in mzTab-M, the database
  identifier prefix `MS-FINDER:` show it.
- On the GC-MS branch the `340` test checks retention agreement only. Automatic annotation keeps only
  EI-MS matches, so this matters only for a reference assigned by hand from the compound search.
- The code is written to the `Annotation tag (VS1.0)` column of the alignment (`.mdalign`) text export
  and, except for GC-MS, of the per-file peak (`.mdpeak`) text export; to the mzTab-M `reliability`
  column, declared as `[MS, MS:1003032, compound identification confidence code in MS-DIAL, ]`, where a
  hit from a user-defined text library is written as "annotated by user-defined text library" instead;
  and to `annotation_tag_vs1` in `.mdcandidate.tsv`, a per-candidate export that only the LC-MS console
  writes and that is off by default. It is not shown in the graphical interface and is not written to
  MassBank or MSP exports.
- Earlier releases differ: `450` first appears in MS-DIAL 5.5.261105, and releases before 5.5.261018
  could give `3x0` or `520` for a reference with no retention time. MS-DIAL 4 capped the retention
  verdict at 0.5 min.

## Support in MS-DIAL

The aim is for MS-DIAL to write the claim and the evidence inventory in its standard outputs, so that
the inventory is visible wherever MS-DIAL results are used. MS-DIAL can supply only part of it:
context, other spectroscopy and many in-silico results come from elsewhere. MS-DIAL does not yet write
the claim or the inventory itself. The middle column describes the released MS-DIAL 5.5.261105, checked
against the source at commit `92cbf0d`; the one item that exists only in source is marked.

| Tag | In MS-DIAL 5.5.261105 | Not yet available |
| --- | --- | --- |
| `RS` | none | A declaration that a library holds authentic standards measured on the same method |
| `FM` | Precursor m/z agreement, adduct type, MS1 isotope pattern, the reference's formula | Isotope similarity in the alignment table and mzTab-M; formula assignment, which MS-FINDER provides |
| `SL` | MS/MS similarity scores (cosine, not squared), matched peaks, `Evidence source`; up to three ranked candidates per annotator in `.mdcandidate.tsv` and mzTab-M | Library name, version and checksum per row; reference collision energy |
| `DF` | Rule-based lipid annotation (x = 0–2 above); diagnostic-rule results written by MS-FINDER (in the MS-FINDER source, not yet released) | Which rule, ion or neutral loss supported a lipid annotation |
| `RT` | RT or RI verdict, reference RT, RT similarity; a retention-order check for lipid annotations | Whether the reference RT was measured on the same method or predicted |
| `IM` | CCS, reference CCS, CCS verdict and similarity | Measured versus predicted reference CCS |
| `IS` | Matches to predicted spectra (`450`); MS-FINDER structures accepted in the graphical interface (`Evidence source`) | Tool version, rank and score type per candidate; in-silico tools run from the console |
| `MN` | Molecular networks from the graphical interface; export for GNPS feature-based molecular networking | Network evidence linked to individual annotations |
| `HO` | none | Homologous-series evidence in the exports |
| `CO`, `UN`, `OS` | none | Outside MS-DIAL today |

## How to cite

- For the Level 3 sublevels `L3-SP` and `L3-CP`, cite Matsuda et al. (2026) [1], which introduced them
  as Level 3-SP and Level 3-CP (see [Claim](#claim)). The preprint's Level 3-SP names the most likely
  planar structure without ruling out others; listing the remaining candidates and stating how the list
  was bounded is a proposal of this page.
- For the framework as a whole (the three layers, `L3-SC` and `L3-SI`, the evidence inventory and the
  study-level criteria) and for MS-DIAL's annotation code, cite this page by its URL at the tag below,
  with the version number:

<!-- AUTHOR DECISION: credit wording -->
<!-- The tag annotation-framework-v0.1 will be created on the merge commit of this page; the link below resolves only after that. -->
> systemsomicslab. Annotation evidence framework: claim, evidence inventory and study-level criteria.
> Working proposal, version 0.1 (25 November 2026).
> https://github.com/systemsomicslab/msdial_spectrum_catalog/blob/annotation-framework-v0.1/docs/annotation_evidence_framework.md

## Relation to other frameworks

- **Schymanski et al. (2014).** The levels are theirs. Their distinction between Level 2a (library) and
  2b (diagnostic) is carried here by the `SL` and `DF` tags rather than by sublevels, and the claim tags
  subdivide their Level 3 by how specific the proposal is.
- **Metabolomics Standards Initiative (Sumner et al. 2007).** Its levels separate identified compounds,
  putatively annotated compounds, putatively characterised compound classes and unknowns. `L3-CP` is
  close to its class level.
- **Evidence codes attached to a level.** Other proposals add letters or codes for the kinds of evidence
  used to a confidence level, and address the same need. The arrangement here keeps the inventory apart
  from both the level and the claim, and suggests stating the thresholds once per study, so that
  evidence can be re-tagged when criteria change.
- Lipid examples use the shorthand notation of Liebisch et al. (2020).

## References

1. Matsuda F, Hayasaka R, Torigoe T, Takahashi M, Takahashi Y, Yoshizawa AC, Tsugawa H, Okuda S,
   Neumann S, Schymanski EL, Izumi Y, Hirayama A. MassBank Human: The annotated human metabolome.
   bioRxiv 2026.10.10.758200 (2026). doi:10.64898/2026.10.10.758200.
   <https://www.biorxiv.org/content/10.64898/2026.10.10.758200v1>
2. Schymanski EL, Jeon J, Gulde R, Fenner K, Ruff M, Singer HP, Hollender J. Identifying small molecules
   via high resolution mass spectrometry: communicating confidence. Environ Sci Technol 48, 2097–2098
   (2014). <https://doi.org/10.1021/es5002105>
3. Sumner LW, Amberg A, Barrett D, et al. Proposed minimum reporting standards for chemical analysis.
   Metabolomics 3, 211–221 (2007). <https://doi.org/10.1007/s11306-007-0082-2>
4. Liebisch G, Fahy E, Aoki J, et al. Update on LIPID MAPS classification, nomenclature, and shorthand
   notation for MS-derived lipid structures. J Lipid Res 61, 1539–1555 (2020).
   <https://doi.org/10.1194/jlr.S120001025>

## Version history

- 0.1 (25 November 2026): first public summary, as a working proposal.
