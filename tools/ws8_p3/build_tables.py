"""WS8 P3 shared tables: the draft rule set's elements, adducts, units and flavonoid aglycone list (stdlib only).

The rows below were transcribed for the draft rule set (ruleset-draft-20261118) from the normalised candidate
records (candidates-20261118: the DOI and locator of each record, and the source field of the NPR records) and
from the MS-DIAL adduct resource. Every row is `as_published` with a page-level source. The script writes them in
storage form (rule_format.py) after cross-checking the aglycone formulas against the ions the records print.

    python -B tools/ws8_p3/build_tables.py write <RuleSet> [--rule-set-version 0.1.0-draft]
    python -B tools/ws8_p3/build_tables.py compare <RuleSet>

`write` creates ruleset.json, tables/elements.json, tables/adducts.json, tables/units.json and
scaffolds/flavonoid_aglycones.json under <RuleSet>. A file that already holds other content is never overwritten:
the run stops before writing anything and names the file. `compare` reads those files and reports, per table, the
built rows that are equal, differ or are missing, and the rows the builder does not produce (for example unit and
adduct requests merged into the tables later). The builder reproduces only its own rows.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

try:
    from . import rule_format as rf
except ImportError:  # run as a script
    import rule_format as rf

DEFAULT_RULE_SET_VERSION = "0.1.0-draft"
MSDIAL = "MsdialWorkbench 563cc1c86 src/Common/CommonStandard/Resources/AdductIonResource_{}.txt"


def cite(doi, locator):
    return {"doi": doi, "locator": locator}


def msdial(mode, line, printed):
    return {"reference": MSDIAL.format(mode), "locator": f"line {line} ({printed})"}


ELEMENTS = {"format": "msfinder-diagnostic-elements", "format_version": 1, "electron": 0.000549,
            "elements": {"C": 12.0, "Cl": 34.96885268, "H": 1.007825032, "K": 38.96370668, "N": 14.003074,
                         "Na": 22.98976928, "O": 15.99491462, "P": 30.97376163, "S": 31.972071}}


def adduct(name, charge, add, remove, common, aliases, source):
    return {"adduct": name, "aliases": aliases, "charge": charge, "multiplier": 1, "add": add, "remove": remove,
            "common": common, "origin": "as_published", "source": source}


ADDUCTS = [
    adduct("[M+H]+", 1, "H", None, True, ["M+H"], [msdial("Positive", 2, "[M+H]+")]),
    adduct("[M+NH4]+", 1, "H4N", None, True, ["M+NH4"], [msdial("Positive", 3, "[M+NH4]+")]),
    adduct("[M+Na]+", 1, "Na", None, True, ["M+Na"], [msdial("Positive", 4, "[M+Na]+")]),
    adduct("[M+K]+", 1, "K", None, True, ["M+K"], [msdial("Positive", 6, "[M+K]+")]),
    adduct("[M+H-H2O]+", 1, "H", "H2O", False, ["[M-H2O+H]+", "M+H-H2O", "M-H2O+H"],
           [msdial("Positive", 9, "[M+H-H2O]+")]),
    # intrinsic charge (flavylium and other permanent cations): not in the MS-DIAL adduct resource
    adduct("[M]+", 1, None, None, False, ["[M+]+", "M+"],
           [cite("10.1021/jf981242+", "Results, ES-MS Anthocyanin Profile"),
            cite("10.1016/j.chroma.2005.07.036", "Sec. 3; Figures 2B-C, 3B-C, 4B-F, 6B-G")]),
    adduct("[M-H]-", -1, None, "H", True, ["M-H"], [msdial("Negative", 2, "[M-H]-")]),
    adduct("[M+Cl]-", -1, "Cl", None, True, ["M+Cl"], [msdial("Negative", 5, "[M+Cl]-")]),
    adduct("[M+HCOO]-", -1, "CHO2", None, True,
           ["[M+FA-H]-", "[M-H+HCOOH]-", "[M+HCOOH-H]-", "M+FA-H", "M+HCOO"],
           [msdial("Negative", 7, "[M+FA-H]-")]),
    adduct("[M+CH3COO]-", -1, "C2H3O2", None, True,
           ["[M+Hac-H]-", "[M-H+CH3COOH]-", "[M+CH3COOH-H]-", "M+Hac-H", "M+CH3COO"],
           [msdial("Negative", 8, "[M+Hac-H]-")]),
]

JMS585 = cite("10.1002/jms.585", "Flavonoid O-glycosides; Fig. 4a")
CHROMA2005 = "10.1016/j.chroma.2005.07.036"
RCM1697 = "10.1002/rcm.1697"
JF026187Q = cite("10.1021/jf026187q", "Results, General LC-MS Fragmentation; Table 2")


def unit(name, cls, residue, source):
    return {"name": name, "class": cls, "residue": residue, "origin": "as_published", "source": source}


UNITS = [
    unit("Hex", "sugar", "C6H10O5", [cite("10.1016/j.jasms.2003.07.002", "Sec. Negative Ion Mode; Tables 1-2; Figs 2, 4"),
                                     JMS585]),
    unit("dHex", "sugar", "C6H10O4", [cite(CHROMA2005, "Sec. 3; Figures 2D-F, 3D, 3F, 3H-I"), JMS585]),
    unit("Pen", "sugar", "C5H8O4", [cite(CHROMA2005, "Sec. 3; Figure 4G-H; Figure 5"), JMS585]),
    unit("HexA", "sugar", "C6H8O6", [cite("10.1002/pca.1189", "Results, fragmentation of flavone di- and "
                                                             "triglycosides; Table 2; Figure 2"), JMS585]),
    unit("malonyl", "acyl", "C3H2O3", [cite("10.1002/jms.1344", "Results, malonylated flavonoid glycosides; "
                                                               "Table 2; Fig. 3")]),
    unit("acetyl", "acyl", "C2H2O", [cite("10.1002/jms.3413", "Results (acetylated flavone C-glycosides); Table 4")]),
    unit("coumaroyl", "acyl", "C9H6O2", [cite("10.1021/jf048068b", "Table 1 (Concord grape peaks 17, 19, 21, 23, 24)"),
                                         cite("10.1002/jms.3803", "Results")]),
    unit("caffeoyl", "acyl", "C9H6O3", [cite(RCM1697, "Table 2; Fig. 3a, 3c"), JF026187Q]),
    unit("feruloyl", "acyl", "C10H8O3", [JF026187Q]),
    unit("sinapoyl", "acyl", "C11H10O4", [cite("10.1002/rcm.4639", "Table 2; Results; Fig. 5")]),
    unit("galloyl", "acyl", "C7H4O4", [cite(RCM1697, "Table 2; Table 1; Fig. 4")]),
    unit("sulfate", "conjugate", "O3S", [cite("10.1021/acs.analchem.6b03671", "Results (CID); Table 1; Fig. 3"),
                                         cite("10.1002/jms.635", "Results, CID conditions")]),
]

# Record sources (DOI + locator as in candidates-20261118/records.tsv, or the NPR record's source field).
NPR = "10.1039/d1np00014d"
NPR_FLAVANONE = cite(NPR, "p.1740, Fig. 5 Flavanone panel")
NPR_FLAVONOL = cite(NPR, "p.1740, Fig. 5 Flavonol panel")
NPR_FLAVONE = cite(NPR, "p.1740, Fig. 5 Flavone panel")
NPR_C_GLYC = cite(NPR, "p.1740, Fig. 5 Flavone C-glycosides panel")
M12030593_T2 = cite("10.3390/12030593", "Table 2")
M12030593_T23 = cite("10.3390/12030593", "Table 2, Table 3")
M29225246_153 = cite("10.3390/molecules29225246", "Table 2; Tables A5, A6, A8; Section 2 text")
M29225246_167 = cite("10.3390/molecules29225246", "Table 2; Table A4")
M29225246_OME = cite("10.3390/molecules29225246", "Section 2 text; Tables 1, 4, A4, A5, A10, A11")
M29225246_DAI = cite("10.3390/molecules29225246", "Section 2 text; Tables A8, A9")
M29225246_FIS = cite("10.3390/molecules29225246", "Tables A6, A7; Section 2 text")
M29225246_F3OL = cite("10.3390/molecules29225246", "Table 2; Table 3 (structures as images); Section 2 text")
PLANTS12193388 = cite("10.3390/plants12193388", "Section 2")
JCHROMB_153 = cite("10.1016/j.jchromb.2026.125348", "Sec. 3.3.1; Table 2; Fig. 3A")
JCHROMB_WOG = cite("10.1016/j.jchromb.2026.125348", "Sec. 3.4.2; Fig. 3B")
JCHROMB_ACA = cite("10.1016/j.jchromb.2026.125348", "Sec. 3.5; Table 3; Fig. 5A")
JCHROMB_IRI = cite("10.1016/j.jchromb.2026.125348", "Sec. 3.5; Table 3; Fig. 5B")
JCHROMB_LIQ = cite("10.1016/j.jchromb.2026.125348", "Sec. 3.6.3; Sec. 3.9; Fig. 7B")
JCHROMB_TAX = cite("10.1016/j.jchromb.2026.125348", "Sec. 3.8; Fig. 8")
RCM983 = cite("10.1002/(sici)1097-0231(199708)11:12<1357::aid-rcm983>3.0.co;2-9", "Tables 1-2; Scheme 1")
IJMS_KAE = cite("10.1016/j.ijms.2005.09.011", "Sec. 4.1.2, 4.2.2.2; Scheme 5; Fig. 2c; Table 6")
IJMS_NAR = cite("10.1016/j.ijms.2005.09.011", "Sec. 4.2.2.1-4.2.2.3; Table 3; Table 6")
JF062431S = cite("10.1021/jf062431s", "Tables 2-6")
METAB0287 = cite("10.1007/s11306-011-0287-2", "Table 1 footnotes a-c; Methods")
TPC058040 = cite("10.1105/tpc.108.058040", "Table 1; Methods")
JASMS2007 = cite("10.1016/j.jasms.2007.11.004", "Results, negative-ion MSn section; Figure 1 and Figure 3")
EJMS_FLAVANONE = cite("10.1255/ejms.1206", "Table 5; Scheme 1; CID of flavanones and chalcones")
EJMS_FLAVANONOL = cite("10.1255/ejms.1206", "CID of flavanonol; Scheme 1")
EJMS_ISO = cite("10.1255/ejms.1206", "Table 4; CID of isoflavones")
JASMS2001_FLAVANONE = cite("10.1016/s1044-0305(01)00226-4", "Table 3; Scheme 4; Results, flavanones")
AC102546X = cite("10.1021/ac102546x", "Results, hydroxylation patterns; Supplemental Fig. 3")
MOL24213929 = cite("10.3390/molecules24213929", "Sec. 2.2.1; Table 2")
CHROMA2024 = cite("10.1016/j.chroma.2024.465548", "Sec. 3.2.1; Table 1 (no. 12); Fig. 2")
JASMS2016 = cite("10.1007/s13361-015-1240-9", "Results, basic fragmentation patterns")
JMS4793 = cite("10.1002/jms.4793", "Results, Study on fragmentation; Fig. 2; Table S2")


def rda(pos=(), neg=()):
    """1,3A ion (A-ring retro-Diels-Alder) as the record gives it for this aglycone."""
    return {"rda_13a": {"pos": list(pos), "neg": list(neg)}}


A_DIOH_POS, A_DIOH_NEG = "C7H5O4", "C7H3O4"  # 5,7-dihydroxy A ring with C-4 carbonyl: 153.0182 / 151.0037
A_OH_POS, A_OH_NEG = "C7H5O3", "C7H3O3"  # 7-hydroxy A ring: 137.0233 / 135.0088
A_OME_POS = "C8H7O4"  # 5-hydroxy-7-methoxy (or 6-methoxy-7-hydroxy) A ring: 167.0339
A_F3OL_POS = "C7H7O3"  # flavan-3-ol, no C-4 carbonyl: 139.0390


def entry(sid, formula, cls, source, ions=None, name=None):
    return {"scaffold_id": sid, "name": name or sid.replace("_", " "), "formula": formula, "class": cls,
            "ions": ions or {}, "origin": "as_published", "source": source}


ENTRIES = [
    # flavones
    entry("apigenin", "C15H10O5", "flavone", [NPR_FLAVONE, M12030593_T2, PLANTS12193388],
          rda([A_DIOH_POS], [A_DIOH_NEG])),
    entry("luteolin", "C15H10O6", "flavone", [NPR_FLAVONE, M12030593_T2, JCHROMB_153], rda([A_DIOH_POS])),
    entry("chrysin", "C15H10O4", "flavone", [M12030593_T2, M29225246_153], rda([A_DIOH_POS])),
    entry("acacetin", "C16H12O5", "flavone", [NPR_C_GLYC, JCHROMB_ACA, RCM983], rda([A_DIOH_POS])),
    entry("diosmetin", "C16H12O6", "flavone", [NPR_FLAVONE, M29225246_OME]),
    entry("chrysoeriol", "C16H12O6", "flavone", [METAB0287, RCM983], rda([A_DIOH_POS])),
    entry("tricin", "C17H14O7", "flavone", [NPR_FLAVONE]),
    entry("baicalein", "C15H10O5", "flavone", [JCHROMB_WOG]),
    entry("wogonin", "C16H12O5", "flavone", [JCHROMB_WOG, M29225246_OME]),
    entry("hispidulin", "C16H12O6", "flavone", [M29225246_OME]),
    entry("genkwanin", "C16H12O5", "flavone", [M29225246_167, M29225246_OME], rda([A_OME_POS])),
    entry("dimethoxyluteolin", "C17H14O6", "flavone", [NPR_C_GLYC], name="dimethoxyluteolin (as printed)"),
    # flavonols
    entry("kaempferol", "C15H10O6", "flavonol", [NPR_FLAVONOL, M29225246_153, JCHROMB_153, IJMS_KAE],
          rda([A_DIOH_POS], [A_DIOH_NEG])),
    entry("quercetin", "C15H10O7", "flavonol", [NPR_FLAVONOL, M29225246_153, JCHROMB_153], rda([A_DIOH_POS])),
    entry("myricetin", "C15H10O8", "flavonol", [JF062431S, RCM983], rda([A_DIOH_POS])),
    entry("isorhamnetin", "C16H12O7", "flavonol", [M29225246_153, TPC058040], rda([A_DIOH_POS])),
    entry("rhamnetin", "C16H12O7", "flavonol", [M29225246_167, M29225246_OME], rda([A_OME_POS])),
    entry("galangin", "C15H10O5", "flavonol", [M29225246_153], rda([A_DIOH_POS])),
    entry("fisetin", "C15H10O6", "flavonol", [M29225246_FIS, cite("10.3390/12030593", "Table 2, Results")]),
    entry("morin", "C15H10O7", "flavonol", [M29225246_153, AC102546X], rda([A_DIOH_POS])),
    entry("kaempferide", "C16H12O6", "flavonol", [NPR_FLAVONOL, M29225246_153], rda([A_DIOH_POS])),
    entry("syringetin", "C17H14O8", "flavonol", [NPR_FLAVONOL]),
    entry("quercetagetin", "C15H10O8", "flavonol", [NPR_FLAVONOL]),
    entry("patuletin", "C16H12O8", "flavonol", [NPR_FLAVONOL]),
    entry("dimethoxyflavonol", "C17H14O5", "flavonol", [NPR_FLAVONOL], name="dimethoxyflavonol (as printed)"),
    # flavanones
    entry("naringenin", "C15H12O5", "flavanone",
          [NPR_FLAVANONE, M12030593_T23, PLANTS12193388, IJMS_NAR, JASMS2001_FLAVANONE],
          rda([A_DIOH_POS], [A_DIOH_NEG])),
    entry("eriodictyol", "C15H12O6", "flavanone", [NPR_FLAVANONE, M12030593_T23, JASMS2001_FLAVANONE],
          rda([A_DIOH_POS], [A_DIOH_NEG])),
    entry("hesperetin", "C16H14O6", "flavanone", [NPR_FLAVANONE, M12030593_T23], rda([A_DIOH_POS])),
    entry("liquiritigenin", "C15H12O4", "flavanone", [NPR_FLAVANONE, JCHROMB_LIQ, EJMS_FLAVANONE],
          rda([A_OH_POS], [A_OH_NEG])),
    entry("sakuranetin", "C16H14O5", "flavanone", [NPR_FLAVANONE, M29225246_167], rda([A_OME_POS])),
    entry("isosakuranetin", "C16H14O5", "flavanone", [JF062431S, JASMS2001_FLAVANONE]),
    entry("homoeriodictyol", "C16H14O6", "flavanone", [JASMS2007]),
    # dihydroflavonols
    entry("taxifolin", "C15H12O7", "dihydroflavonol", [M12030593_T23, JCHROMB_TAX], rda([A_DIOH_POS])),
    entry("aromadendrin", "C15H12O6", "dihydroflavonol", [PLANTS12193388, M29225246_153],
          rda([A_DIOH_POS], [A_DIOH_NEG])),
    entry("dihydromyricetin", "C15H12O8", "dihydroflavonol", [EJMS_FLAVANONOL]),
    # isoflavones
    entry("daidzein", "C15H10O4", "isoflavone", [NPR_FLAVONE, M29225246_DAI], rda([A_OH_POS])),
    entry("genistein", "C15H10O5", "isoflavone", [NPR_FLAVONE, M29225246_153, JASMS2016], rda([A_DIOH_POS])),
    entry("glycitein", "C16H12O5", "isoflavone", [NPR_FLAVONE, JF062431S]),
    entry("formononetin", "C16H12O4", "isoflavone", [NPR_FLAVONE, MOL24213929], rda([A_OH_POS])),
    entry("biochanin_a", "C16H12O5", "isoflavone", [METAB0287, EJMS_ISO], name="biochanin A"),
    entry("afrormosin", "C17H14O5", "isoflavone", [NPR_FLAVONE, MOL24213929], rda([A_OME_POS])),
    entry("calycosin", "C16H12O5", "isoflavone", [MOL24213929, CHROMA2024], rda([A_OH_POS])),
    entry("orobol", "C15H10O6", "isoflavone", [JMS4793]),
    entry("irisolidone", "C17H14O6", "isoflavone", [JCHROMB_IRI, METAB0287]),
    # flavan-3-ols
    entry("catechin", "C15H14O6", "flavan_3_ol", [M12030593_T2, M29225246_F3OL], rda([A_F3OL_POS])),
    entry("epicatechin", "C15H14O6", "flavan_3_ol", [M12030593_T2, M29225246_F3OL], rda([A_F3OL_POS])),
    entry("gallocatechin", "C15H14O7", "flavan_3_ol", [M29225246_F3OL]),
    entry("epigallocatechin", "C15H14O7", "flavan_3_ol", [M29225246_F3OL]),
    # chalcones and dihydrochalcones
    entry("isoliquiritigenin", "C15H12O4", "chalcone", [EJMS_FLAVANONE], rda([], [A_OH_NEG])),
    entry("butein", "C15H12O5", "chalcone", [EJMS_FLAVANONE], rda([], [A_OH_NEG])),
    entry("phloretin", "C15H14O5", "dihydrochalcone", [JF062431S]),
]
SCAFFOLD_LIST = "flavonoid_aglycones"

# Cross-checks against the records: NPR Y0 [Agly+H]+ formulas (exact) and nominal aglycone ions printed elsewhere.
NPR_Y0 = {"naringenin": "C15H13O5", "liquiritigenin": "C15H13O4", "sakuranetin": "C16H15O5",
          "eriodictyol": "C15H13O6", "hesperetin": "C16H15O6", "kaempferol": "C15H11O6",
          "dimethoxyflavonol": "C17H15O5", "kaempferide": "C16H13O6", "quercetin": "C15H11O7",
          "quercetagetin": "C15H11O8", "patuletin": "C16H13O8", "syringetin": "C17H15O8", "tricin": "C17H15O7",
          "daidzein": "C15H11O4", "formononetin": "C16H13O4", "apigenin": "C15H11O5", "genistein": "C15H11O5",
          "glycitein": "C16H13O5", "luteolin": "C15H11O6", "afrormosin": "C17H15O5", "diosmetin": "C16H13O6"}
NPR_PLUS_C = {"acacetin": "C17H13O5", "dimethoxyluteolin": "C18H15O6", "apigenin": "C16H11O5",
              "luteolin": "C16H11O6", "diosmetin": "C17H13O6"}
NOMINAL_POS = {"isorhamnetin": 317, "myricetin": 319, "isosakuranetin": 287, "phloretin": 275, "biochanin_a": 285,
               "chrysoeriol": 301, "irisolidone": 315, "acacetin": 285, "calycosin": 285, "taxifolin": 305,
               "rhamnetin": 317, "catechin": 291, "chrysin": 255, "galangin": 271, "morin": 303, "fisetin": 287}
NOMINAL_NEG = {"homoeriodictyol": 301, "dihydromyricetin": 319}
# Ring-cleavage ions with the accurate m/z the records print.
PRINTED_IONS = {A_DIOH_POS: (1, 153.0182), A_DIOH_NEG: (-1, 151.0037), A_OH_POS: (1, 137.0233),
                A_OH_NEG: (-1, 135.0088), A_OME_POS: (1, 167.0339), A_F3OL_POS: (1, 139.0390)}


class CheckError(ValueError):
    """A built row disagrees with an ion the records print."""


def load_elements(doc=None):
    errors = []
    el = rf.load_elements(rf.Ctx(errors, rf.ELEMENTS_FILE), ELEMENTS if doc is None else doc)
    if el is None:
        raise CheckError("elements: " + "; ".join(str(e) for e in errors))
    return el


def check(el, entries=None, y0=None, plus_c=None, nominal_pos=None, nominal_neg=None, printed=None):
    """Cross-check the aglycone entries against the ions the records give; raise CheckError on the first
    disagreement. The defaults are this module's tables; tests pass synthetic ones."""
    entries = ENTRIES if entries is None else entries
    y0 = NPR_Y0 if y0 is None else y0
    plus_c = NPR_PLUS_C if plus_c is None else plus_c
    nominal_pos = NOMINAL_POS if nominal_pos is None else nominal_pos
    nominal_neg = NOMINAL_NEG if nominal_neg is None else nominal_neg
    printed = PRINTED_IONS if printed is None else printed
    pda = el.pda
    seen = set()
    for e in entries:
        sid = e["scaffold_id"]
        seen.add(sid)
        a = rf.parse_formula(e["formula"], pda)
        if sid in y0 and rf.hill_format(rf.comp_add(a, {"H": 1})) != y0[sid]:
            raise CheckError(f"{sid}: [A+H]+ is not {y0[sid]}")
        if sid in plus_c and rf.hill_format(rf.comp_add(a, {"H": 1, "C": 1})) != plus_c[sid]:
            raise CheckError(f"{sid}: [A+C+H]+ is not {plus_c[sid]}")
        if sid in nominal_pos and round(rf.mz_of(el.ion_pda(rf.comp_add(a, {"H": 1}), 1))) != nominal_pos[sid]:
            raise CheckError(f"{sid}: [A+H]+ is not nominal m/z {nominal_pos[sid]}")
        if sid in nominal_neg and round(rf.mz_of(el.ion_pda(rf.comp_add(a, {"H": 1}, -1), -1))) != nominal_neg[sid]:
            raise CheckError(f"{sid}: [A-H]- is not nominal m/z {nominal_neg[sid]}")
        for key, modes in e["ions"].items():
            for mode, formulas in modes.items():
                for f in formulas:
                    rf.parse_formula(f, pda)
    unknown = sorted((set(y0) | set(plus_c) | set(nominal_pos) | set(nominal_neg)) - seen)
    if unknown:
        raise CheckError(f"cross-checks name entries that do not exist: {unknown}")
    for f, (z, mz) in printed.items():
        got = rf.mz_of(el.ion_pda(rf.parse_formula(f, pda), z))
        if abs(got - mz) >= 1e-4:
            raise CheckError(f"{f}: computed m/z {got:.6f}, printed {mz}")


def built_files(rule_set_version=DEFAULT_RULE_SET_VERSION):
    """{posix path relative to RuleSet/: object} for every file the builder writes."""
    return {
        rf.RULESET_FILE: {"format": "msfinder-diagnostic-ruleset", "format_version": 1,
                          "rule_set_version": rule_set_version},
        rf.ELEMENTS_FILE: ELEMENTS,
        rf.ADDUCTS_FILE: {"format": "msfinder-diagnostic-adducts", "format_version": 1, "adducts": ADDUCTS},
        rf.UNITS_FILE: {"format": "msfinder-diagnostic-units", "format_version": 1, "units": UNITS},
        f"scaffolds/{SCAFFOLD_LIST}.json": {"list_id": SCAFFOLD_LIST, "entries": ENTRIES},
    }


def write(root, rule_set_version=DEFAULT_RULE_SET_VERSION):
    """-> ([(status, path)], ok). Writes nothing when any target holds other content."""
    root = Path(root)
    check(load_elements())
    plan = []
    for rel, obj in built_files(rule_set_version).items():
        data = rf.storage_bytes(obj)
        rf.load_strict(data)
        path = root / rel
        if not path.exists():
            plan.append(("write", rel, path, data))
        elif path.read_bytes() == data:
            plan.append(("unchanged", rel, path, data))
        else:
            plan.append(("differs", rel, path, data))
    if any(p[0] == "differs" for p in plan):
        return [(s, rel) for s, rel, _, _ in plan], False
    for status, rel, path, data in plan:
        if status == "write":
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
    return [(s, rel) for s, rel, _, _ in plan], True


ROW_KEYS = {rf.ADDUCTS_FILE: ("adducts", "adduct"), rf.UNITS_FILE: ("units", "name"),
            f"scaffolds/{SCAFFOLD_LIST}.json": ("entries", "scaffold_id")}


def _same(a, b):
    return rf.canonical_bytes(a) == rf.canonical_bytes(b)


def compare(root, rule_set_version=DEFAULT_RULE_SET_VERSION):
    """-> (report lines, ok). ok: every built row is present and equal, and the built rows lead each table."""
    root = Path(root)
    lines, ok = [], True
    for rel, built in built_files(rule_set_version).items():
        path = root / rel
        if not path.exists():
            lines.append(f"missing_file\t{rel}")
            ok = False
            continue
        try:
            have = rf.load_strict(path.read_bytes())
        except rf.RuleFormatError as e:
            lines.append(f"unreadable\t{rel}\t{e.code}")
            ok = False
            continue
        if rel not in ROW_KEYS:
            same = _same(have, built)
            lines.append(f"file\t{rel}\t{'equal' if same else 'differs'}")
            ok = ok and same
            continue
        list_key, id_key = ROW_KEYS[rel]
        head_same = _same({k: v for k, v in have.items() if k != list_key},
                          {k: v for k, v in built.items() if k != list_key})
        have_rows = have.get(list_key) if isinstance(have.get(list_key), list) else []
        have_by = {r.get(id_key): r for r in have_rows if isinstance(r, dict)}
        built_ids = [r[id_key] for r in built[list_key]]
        equal = [i for i, r in zip(built_ids, built[list_key]) if i in have_by and _same(have_by[i], r)]
        differ = [i for i, r in zip(built_ids, built[list_key]) if i in have_by and not _same(have_by[i], r)]
        missing = [i for i in built_ids if i not in have_by]
        extra = [r.get(id_key) for r in have_rows if not isinstance(r, dict) or r.get(id_key) not in built_ids]
        lead = [r.get(id_key) if isinstance(r, dict) else None for r in have_rows[:len(built_ids)]] == built_ids
        lines.append(f"file\t{rel}\tbuilt={len(built_ids)}\tequal={len(equal)}\tdiffer={len(differ)}"
                     f"\tmissing={len(missing)}\textra={len(extra)}\tbuilt_rows_lead={str(lead).lower()}"
                     f"\theader={'equal' if head_same else 'differs'}")
        lines += [f"differ\t{rel}\t{i}" for i in differ] + [f"missing\t{rel}\t{i}" for i in missing]
        lines += [f"extra\t{rel}\t{i}" for i in extra]
        ok = ok and head_same and not differ and not missing and lead
    return lines, ok


def main(argv=None):
    ap = argparse.ArgumentParser(prog="build_tables.py", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, text in (("write", "write the tables into a RuleSet directory (never overwrites other content)"),
                       ("compare", "compare a RuleSet directory's tables with the built rows")):
        p = sub.add_parser(name, help=text)
        p.add_argument("ruleset_dir")
        p.add_argument("--rule-set-version", default=DEFAULT_RULE_SET_VERSION)
    args = ap.parse_args(argv)
    root = Path(args.ruleset_dir)
    try:
        if args.cmd == "write":
            if root.exists() and not root.is_dir():
                print(f"not a directory: {root}", file=sys.stderr)
                return 2
            report, ok = write(root, args.rule_set_version)
            for status, rel in report:
                print(f"{status}\t{rel}")
            if not ok:
                print("refused: a file holds other content; nothing was written (see compare)", file=sys.stderr)
            return 0 if ok else 1
        if not root.is_dir():
            print(f"not a directory: {root}", file=sys.stderr)
            return 2
        lines, ok = compare(root, args.rule_set_version)
        for line in lines:
            print(line)
        print("OK built rows reproduced" if ok else "DIFFERS")
        return 0 if ok else 1
    except CheckError as e:
        print(f"cross-check failed: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
