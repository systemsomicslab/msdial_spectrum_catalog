"""WS8 P3 candidate normaliser, on synthetic records only (never the literature files)."""
import contextlib
import csv
import importlib.util
import io
import json
from pathlib import Path
import re
import tempfile
import unittest

from tools.ws8_p3 import normalise_candidates as nc

HAS_RDKIT = importlib.util.find_spec("rdkit") is not None
TOOL_DIR = Path(nc.__file__).resolve().parent


def record(rid, **kw):
    base = {"id": rid, "family": "flavonoid_O_glycoside", "scope": "synthetic", "kind": "fragment_ion",
            "value": "", "ion_mode": "positive", "adduct": "[M+H]+", "instrument": "Agilent 6530 Q-TOF",
            "standards_based": False, "author_reviewed": False}
    base.update(kw)
    return base


def stage2(rid, paper, **kw):
    return record(rid, paper_id=paper, citation=f"Author A 2020, J Synth, doi:10.1000/{paper.upper()}", **kw)


def quiet_main(argv):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return nc.main(argv)


def build(inputs):
    return nc.build({label: (Path(f"{label}.json"), recs) for label, recs in inputs.items()}, nc.load_maps())


class KindMapTests(unittest.TestCase):
    def setUp(self):
        self.kinds = nc.load_kind_map()

    def test_vocabularies_map_onto_one(self):
        self.assertEqual(self.kinds[("stage2", "fragment_ion")]["kind"], "product_ion")
        self.assertEqual(self.kinds[("stage3", "fragment_ion")]["kind"], "product_ion")
        self.assertEqual(self.kinds[("npr", "fragment_ion")]["kind"], "product_ion")
        self.assertEqual(self.kinds[("stage3", "adduct_behaviour")]["kind"], "adduct_behaviour")
        self.assertEqual(self.kinds[("stage2", "precondition")]["kind"], "precondition")

    def test_negative_evidence_kinds(self):
        self.assertEqual(self.kinds[("stage3", "diagnostic_absence")]["kind"], "absence")
        self.assertEqual(self.kinds[("stage2", "veto")]["kind"], "veto")
        for key in (("stage3", "diagnostic_absence"), ("stage2", "veto")):
            self.assertEqual(self.kinds[key]["evidence_sense"], "negative")

    def test_other_goes_to_manual_queue_and_is_never_guessed(self):
        self.assertEqual(self.kinds[("stage3", "other")]["kind"], "manual")
        result = build({"stage3": [stage2("s3-x-001", "p3x", kind="other", formula="C9H7O3+",
                                          value="a caveat that is not one ion")]})
        row = result["records"][0]
        self.assertEqual((row["kind"], row["kind_original"], row["manual_queue"]), ("manual", "other", True))
        self.assertEqual([m["id"] for m in result["manual"]], ["s3-x-001"])
        self.assertEqual(result["manual"][0]["reason"], "kind_other")

    def test_unknown_kind_is_unmapped_and_queued(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", kind="mystery")]})
        self.assertIn(("kind", "stage2:mystery"), result["unmapped"])
        self.assertEqual(result["manual"][0]["reason"], "kind_unmapped")


class FamilyMapTests(unittest.TestCase):
    def setUp(self):
        self.fam = nc.load_family_map()

    def test_synonyms_unify(self):
        for label in ("flavone_C_glycoside", "flavonoid_C_glycoside", "flavone_c_glycoside"):
            self.assertEqual(self.fam[label]["family"], "flavonoid_c_glycoside")
        self.assertEqual(self.fam["triterpene_saponin"]["family"], "triterpenoid_saponin")
        self.assertEqual(self.fam["alkaloid_MIA"]["family"], self.fam["alkaloid_monoterpene_indole"]["family"])
        self.assertEqual(self.fam["sulfur_metabolites"]["family"], self.fam["sulfur_metabolite"]["family"])
        self.assertEqual(self.fam["acyl_phenylpropanoid"]["family"], self.fam["phenylpropanoid_acyl"]["family"])
        self.assertEqual(self.fam["flavonoid_O_glycoside"]["family"], "flavonoid_o_glycoside")

    def test_distinct_classes_keep_labels_and_share_a_group(self):
        self.assertNotEqual(self.fam["flavone_CC_glycoside"]["family"], "flavonoid_c_glycoside")
        self.assertEqual(self.fam["flavone_CC_glycoside"]["family_group"], "flavonoid_c_glycoside")
        labels = {self.fam[k]["family"] for k in ("flavonol", "flavanone", "flavonoid_aglycone",
                                                  "flavonoid_O_glycoside")}
        self.assertEqual(len(labels), 4)
        groups = {self.fam[k]["family_group"] for k in ("flavonol", "flavanone", "flavonoid_aglycone",
                                                        "flavonoid_O_glycoside")}
        self.assertEqual(groups, {"flavonoid_glycoside_aglycone"})

    def test_scopes(self):
        for label in ("lipid", "lipid_sn_position", "lipid_double_bond_position"):
            self.assertEqual(self.fam[label]["scope"], "out_of_scope")
        for label in ("tool_note", "general_heuristic", "generic_substructure", "general_substructure"):
            self.assertEqual(self.fam[label]["scope"], "generic")
        self.assertEqual(self.fam["glucosinolate"]["scope"], "in_scope")
        self.assertTrue(all(re.fullmatch(r"[a-z0-9_]+", v["family"]) for v in self.fam.values()))

    def test_original_kept_and_unknown_family_unmapped(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", family="flavone_C_glycoside"),
                                   stage2("s2-x-002", "p3x", family="not_a_family")]})
        rows = {r["id"]: r for r in result["records"]}
        self.assertEqual(rows["s2-x-001"]["family_original"], "flavone_C_glycoside")
        self.assertEqual(rows["s2-x-001"]["family"], "flavonoid_c_glycoside")
        self.assertEqual(rows["s2-x-002"]["scope"], "unmapped")
        self.assertIn(("family", "not_a_family"), result["unmapped"])


class IonModeAndAdductTests(unittest.TestCase):
    def test_ion_mode(self):
        self.assertEqual(nc.normalise_ion_mode("positive"), ("positive", False))
        self.assertEqual(nc.normalise_ion_mode("negative (inferred)"), ("negative", True))
        self.assertEqual(nc.normalise_ion_mode("both"), ("both", False))
        self.assertEqual(nc.normalise_ion_mode("positive and negative (Table 1 +/+)")[0], "both")
        self.assertEqual(nc.normalise_ion_mode("negative (341/297); positive (299)")[0], "both")
        self.assertEqual(nc.normalise_ion_mode("negative (one search; positive mode added)")[0], "negative")
        self.assertEqual(nc.normalise_ion_mode("positive only (Table 1 pos +, neg -)")[0], "positive")
        self.assertEqual(nc.normalise_ion_mode("not stated (positive inferred)"), ("unknown", False))
        self.assertEqual(nc.normalise_ion_mode(None), ("unknown", False))
        self.assertIsNone(nc.normalise_ion_mode("sideways"))

    def test_adduct_canonical_strings(self):
        self.assertEqual(nc.parse_adducts("[M−H]−")[0], ["[M-H]-"])
        self.assertEqual(nc.parse_adducts("[M–H]–")[0], ["[M-H]-"])
        self.assertEqual(nc.parse_adducts("M+")[0], ["[M]+"])
        self.assertEqual(nc.parse_adducts("M+ (flavylium)")[0], ["[M]+"])
        self.assertEqual(nc.parse_adducts("not stated; fits [M+H]+")[0], ["unstated"])
        self.assertEqual(nc.parse_adducts(None)[0], ["unstated"])
        self.assertEqual(nc.parse_adducts("[M-H]-/[M-2H]2-")[0], ["[M-H]-", "[M-2H]2-"])
        self.assertEqual(nc.parse_adducts("[M+FA-H]-")[0], ["[M+HCOO]-"])
        self.assertEqual(nc.parse_adducts("[M+NH4]+")[0], ["[M+NH4]+"])
        self.assertEqual(nc.parse_adducts("[M+CH3COO]-")[0], ["[M+CH3COO]-"])
        self.assertEqual(nc.parse_adducts("[M+H]+ -> [Agl+H]+")[0], ["[M+H]+"])
        self.assertEqual(nc.parse_adducts("Y0- from [M-H]-")[0], ["[M-H]-"])
        self.assertEqual(nc.parse_adducts("[M]+ (printed as '[M-H]+')")[0], ["[M]+"])
        self.assertEqual(nc.parse_adducts("[(M+Na)-aglycone]+")[0], ["[(M+Na)-aglycone]+"])
        self.assertEqual(nc.parse_adducts("[M-H]- or M-.")[0], ["[M-H]-", "[M]-•"])
        adducts, flags = nc.parse_adducts("see value")
        self.assertIsNone(adducts)

    def test_adduct_original_kept_and_unmapped_logged(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", adduct="[M−H]−", ion_mode="negative"),
                                   stage2("s2-x-002", "p3x", adduct="molecular ion")]})
        rows = {r["id"]: r for r in result["records"]}
        self.assertEqual(rows["s2-x-001"]["adduct_original"], "[M−H]−")
        self.assertEqual(rows["s2-x-001"]["adduct"], ["[M-H]-"])
        self.assertIn(("adduct", "molecular ion"), result["unmapped"])

    def test_both_expands_by_polarity(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", ion_mode="both", adduct="[M+H]+ / [M−H]−",
                                          kind="neutral_loss", formula="C6H10O5 (Glc); C9H7O3+")]})
        terms = {t["term_id"]: t for t in result["terms"]}
        self.assertEqual(set(terms), {"s2-x-001#t1.pos", "s2-x-001#t1.neg", "s2-x-001#t2.pos"})
        self.assertEqual(terms["s2-x-001#t1.pos"]["adduct"], ["[M+H]+"])
        self.assertEqual(terms["s2-x-001#t1.neg"]["adduct"], ["[M-H]-"])
        self.assertEqual(terms["s2-x-001#t1.neg"]["ion_mode"], "negative")
        self.assertTrue(all(t["expanded_from_both"] for t in terms.values()))
        self.assertIn("kind_conflict", terms["s2-x-001#t2.pos"]["flags"])
        self.assertEqual(result["records"][0]["ion_mode"], "both")

    def test_both_without_adduct_of_a_polarity(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", ion_mode="both", adduct="[M+H]+", formula="H2O")]})
        terms = {t["ion_mode"]: t for t in result["terms"]}
        self.assertEqual(terms["negative"]["adduct"], ["unstated"])


class FormulaTests(unittest.TestCase):
    def tokens(self, text, **kw):
        return nc.tokenise_formulas(text, **kw)[0]

    def test_multiple_tokens_labels_and_hill(self):
        toks = self.tokens("C6H8O6 (GlcA loss); C6H10O5 (Glc loss), C9H8O4N- (194.0459)")
        self.assertEqual([t["formula"] for t in toks], ["C6H8O6", "C6H10O5", "C9H8NO4-"])
        self.assertEqual([t["label"] for t in toks], ["GlcA loss", "Glc loss", "194.0459"])
        self.assertEqual([t["formula"] for t in self.tokens("CH3OH; SO3; HSO4-")], ["CH4O", "O3S", "HO4S-"])

    def test_charges(self):
        toks = self.tokens("C9H7O3+; C7H11O6−, C21H18O12²⁻")
        self.assertEqual([(t["formula"], t["charge"]) for t in toks],
                         [("C9H7O3+", 1), ("C7H11O6-", -1), ("[C21H18O12]2-", -2)])
        self.assertEqual(self.tokens("C6H11O8S2-")[0]["formula"], "C6H11O8S2-")
        self.assertEqual(self.tokens("C6H10O5-H2O loss")[0]["charge"], 0)

    def test_radicals(self):
        tok = self.tokens("C3H9N•+ (calc 59.0730)")[0]
        self.assertEqual((tok["formula"], tok["radical_printed"], tok["odd_electron"]), ("C3H9N+", True, True))
        tok = self.tokens("computed C15H8O7-. 300.0276")[0]
        self.assertEqual((tok["formula"], tok["radical_printed"]), ("C15H8O7-", True))
        tok = self.tokens("C15H9O7- 301.0354")[0]
        self.assertFalse(tok["odd_electron"])
        tok = self.tokens("C9H7O3•+")[0]
        self.assertTrue(tok["radical_mark_mismatch"])

    def test_unsupported_element(self):
        tok = self.tokens("C6H10O5Li+")[0]
        self.assertEqual(tok["unsupported_elements"], ["Li"])
        self.assertIsNone(tok["mass"])
        result = build({"stage3": [stage2("s3-x-001", "p3x", formula="C6H10O5Li+")]})
        self.assertIn("unsupported_element", result["terms"][0]["flags"])
        self.assertEqual(result["terms"][0]["mass_check"], "no_mass")

    def test_words_notations_and_ion_notation_contents_are_not_terms(self):
        self.assertEqual(self.tokens("K C18H11O4+"), self.tokens("C18H11O4+"))
        self.assertEqual(self.tokens("PC and SI, CQA, Glc, C2 position"), [])
        self.assertEqual(self.tokens("BH2+ (protonated base)"), [])
        self.assertEqual([t["formula"] for t in self.tokens("read as [M-H-CO2]- and [Agl+H-H2O]+")], [])
        self.assertEqual([t["formula"] for t in self.tokens("[C4H5O2]+ 85.0284")], ["C4H5O2+"])
        self.assertEqual([t["formula"] for t in self.tokens("O")], ["O"])
        self.assertEqual(self.tokens("Ag = aglycone"), [])
        self.assertEqual(nc.tokenise_formulas("CH3R3 loss")[1], ["CH3R3"])

    def test_value_extraction_needs_carbon_and_hydrogen(self):
        toks = self.tokens("m/z 163.0390; C9H7O3+ (caffeoyl); 96.9601 HSO4-; CO2 loss; CH2=CH2", from_value=True)
        self.assertEqual([t["formula"] for t in toks], ["C9H7O3+"])
        result = build({"stage2": [stage2("s2-x-001", "p3x", value="m/z 163.0390 (C9H7O3+)"),
                                   stage2("s2-x-002", "p3x", value="loss of 162 Da")]})
        rows = {r["id"]: r for r in result["records"]}
        self.assertEqual((rows["s2-x-001"]["formula_from"], rows["s2-x-001"]["term_status"]), ("value", "has_terms"))
        self.assertEqual(rows["s2-x-002"]["term_status"], "no_formula")
        self.assertEqual(rows["s2-x-002"]["n_term_rows"], 0)

    def test_unparsed_formula_is_logged_not_dropped(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", formula="not printed", value="C9H7O3+")]})
        self.assertEqual(result["records"][0]["term_status"], "formula_unparsed")
        self.assertEqual(result["terms"], [])
        self.assertIn(("formula", "not printed"), result["unmapped"])

    def test_kind_conflict(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", kind="neutral_loss", formula="C9H7O3+"),
                                   stage2("s2-x-002", "p3x", kind="fragment_ion", formula="C6H10O5"),
                                   stage2("s2-x-003", "p3x", kind="fragment_ion", formula="C9H7O3+")]})
        flags = {t["record_id"]: t["flags"] for t in result["terms"]}
        self.assertIn("kind_conflict", flags["s2-x-001"])
        self.assertIn("kind_conflict", flags["s2-x-002"])
        self.assertNotIn("kind_conflict", flags["s2-x-003"])

    def test_charge_against_ion_mode(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", ion_mode="negative", adduct="[M-H]-",
                                          formula="C9H7O3+")]})
        self.assertIn("charge_mode_conflict", result["terms"][0]["flags"])


class MassTests(unittest.TestCase):
    # Hand-computed from the element table and the electron mass, written out independently of the module.
    IONS = {"C9H7O3+": 163.038970, "C7H11O6-": 191.056112, "C5H15NO4P+": 184.073321, "C6H11O9S-": 259.012927,
            "HSO4-": 96.960104, "C21H18O12 2-": 231.040462, "C15H8O7-": 300.027552, "Cl-": 34.969402,
            "C2H3KO2+": 97.976462, "C18H30NaO14+": 493.152776}
    NEUTRALS = {"H2O": 18.010565, "CO2": 43.989829, "C6H10O5": 162.052823, "C3H9N": 59.073499,
                "C4H8O4": 120.042259}

    def test_element_table_matches_p2(self):
        self.assertEqual(nc.ELEMENT_MASS, {"C": 12.0, "H": 1.007825032, "N": 14.003074, "O": 15.99491462,
                                           "P": 30.97376163, "S": 31.972071, "Na": 22.98976928,
                                           "K": 38.96370668, "Cl": 34.96885268})
        self.assertEqual(nc.ELECTRON_MASS, 0.000549)

    @unittest.skipUnless(HAS_RDKIT, "RDKit not installed in this interpreter")
    def test_element_table_matches_rdkit(self):
        from rdkit import Chem
        table = Chem.GetPeriodicTable()
        for el, mass in nc.ELEMENT_MASS.items():
            self.assertAlmostEqual(mass, table.GetMostCommonIsotopeMass(el), delta=1e-9)

    def test_ion_mz(self):
        for text, expected in self.IONS.items():
            tok = nc.tokenise_formulas(text)[0][0]
            self.assertAlmostEqual(tok["mass"], expected, delta=1e-5, msg=text)
            self.assertNotEqual(tok["charge"], 0)

    def test_neutral_mass(self):
        for text, expected in self.NEUTRALS.items():
            tok = nc.tokenise_formulas(text)[0][0]
            self.assertAlmostEqual(tok["mass"], expected, delta=1e-5, msg=text)
            self.assertEqual(tok["charge"], 0)

    def check(self, mass, *printed, note=()):
        return nc.mass_check(mass, nc.printed_numbers(" ".join(printed)), nc.printed_numbers(" ".join(note)))

    def test_product_tolerance_edge_with_printed_precision(self):
        m = self.IONS["C9H7O3+"]
        self.assertEqual(self.check(m, "163.0390")["status"], "matched")
        self.assertEqual(self.check(m, "163.0440")["status"], "matched")  # 5.03 mDa <= 5 + 0.05 mDa
        self.assertEqual(self.check(m, "163.0441")["status"], "printed_mismatch")  # 5.13 mDa > 5.05 mDa
        self.assertEqual(self.check(m, "163.538")["status"], "printed_mismatch")
        self.assertEqual(self.check(m, "163.6")["status"], "no_printed_value")  # beyond 0.5 Da
        self.assertEqual(self.check(m, "no number")["status"], "no_printed_value")

    def test_lowres_printed(self):
        m = self.IONS["C9H7O3+"]
        self.assertEqual(self.check(m, "163.04")["status"], "matched")  # 1.0 mDa <= 5 + 5 mDa
        self.assertEqual(self.check(m, "163.05")["status"], "lowres_printed")  # 2 decimals cannot resolve 5 mDa
        self.assertEqual(self.check(m, "163.1")["status"], "lowres_printed")
        self.assertEqual(self.check(m, "163.05 and 163.0600")["status"], "printed_mismatch")
        self.assertEqual(self.check(m, "15.2% at 163.05%")["status"], "no_printed_value")  # intensities skipped

    def test_two_decimals_resolve_above_500(self):
        self.assertEqual(nc.mass_check(400.0, nc.printed_numbers("400.02"))["status"], "lowres_printed")
        self.assertEqual(nc.mass_check(600.0, nc.printed_numbers("600.02"))["status"], "printed_mismatch")
        big = 941.511540  # C48H77O18-, tolerance 9.4 mDa
        self.assertEqual(self.check(big, "941.52")["status"], "matched")
        self.assertEqual(self.check(big, "941.53")["status"], "printed_mismatch")
        self.assertEqual(self.check(big, "941.6")["status"], "lowres_printed")
        self.assertAlmostEqual(nc.product_tolerance(big), 9.41511540e-3, delta=1e-12)
        self.assertEqual(nc.product_tolerance(100.0), 0.005)

    def test_value_first_then_note(self):
        m = self.IONS["C9H7O3+"]
        self.assertEqual(self.check(m, "printed 163.0300", note=["calc 163.0390"])["status"], "printed_mismatch")
        result = self.check(m, "no number", note=["calc 163.0390"])
        self.assertEqual((result["status"], result["printed_from"]), ("matched", "note"))

    def test_mass_flags_reach_terms(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", formula="C9H7O3+", value="m/z 163.0300")]})
        self.assertEqual(result["terms"][0]["mass_check"], "printed_mismatch")
        self.assertIn("printed_mismatch", result["terms"][0]["flags"])


class InstrumentTests(unittest.TestCase):
    def setUp(self):
        self.rules = nc.load_instrument_map()

    def classify(self, text):
        return nc.classify_instrument(text, self.rules)

    def test_single_classes(self):
        self.assertEqual(self.classify("Agilent 6530 Q-TOF")[:2], ("qtof", "qtof"))
        self.assertEqual(self.classify("Thermo Q Exactive Orbitrap")[:2], ("orbitrap_ft", "orbitrap_ft"))
        self.assertEqual(self.classify("LTQ Orbitrap XL")[:2], ("orbitrap_ft", "orbitrap_ft"))
        self.assertEqual(self.classify("Bruker solariX 7.0 T FTICR")[:2], ("orbitrap_ft", "orbitrap_ft"))
        self.assertEqual(self.classify("Finnigan LCQ Deca ion trap")[:2], ("ion_trap", "ion_trap"))
        self.assertEqual(self.classify("Thermo LTQ linear ion trap")[:2], ("ion_trap", "ion_trap"))
        self.assertEqual(self.classify("AB Sciex 4000 QTRAP")[:2], ("qqq_or_qtrap", "qqq_or_qtrap"))
        self.assertEqual(self.classify("AB Sciex QTRAP (quadrupole linear ion trap)")[:2],
                         ("qqq_or_qtrap", "qqq_or_qtrap"))
        self.assertEqual(self.classify("Waters Xevo TQ-S")[:2], ("qqq_or_qtrap", "qqq_or_qtrap"))

    def test_mixed_counts_toward_the_fragmenting_instrument(self):
        self.assertEqual(self.classify("Finnigan LCQ ion trap (MSn); Waters Q-TOF Ultima for accurate mass")[:2],
                         ("mixed", "ion_trap"))
        self.assertEqual(self.classify("Thermo LTQ-Orbitrap XL (MSn in the LTQ; accurate mass in the Orbitrap)")[:2],
                         ("mixed", "ion_trap"))
        self.assertEqual(self.classify("Orbitrap Fusion (ion-trap CID)")[:2], ("mixed", "ion_trap"))
        self.assertEqual(self.classify("Bruker HCT Ultra and Impact HD")[:2], ("mixed", "mixed"))
        self.assertEqual(self.classify("Sciex API 4000 Q-TRAP; accurate mass on Waters QTOF Ultima")[:2],
                         ("mixed", "qqq_or_qtrap"))

    def test_compilation_unknown_and_unmapped(self):
        self.assertEqual(self.classify("review across ESI with ion trap and QTOF")[0], "compilation")
        self.assertEqual(self.classify("not stated")[0], "unknown")
        self.assertEqual(self.classify(None), ("unknown", "unknown", "missing"))
        self.assertEqual(self.classify("Micromass AutoSpec-oaTOF")[2], "unknown.sector")
        self.assertIsNone(self.classify("C18 column, formic acid"))


class GroupingTests(unittest.TestCase):
    def test_dedup_across_papers_with_citations_and_counts(self):
        recs = {
            "stage2": [
                stage2("s2-a-001", "p3a", family="flavone_C_glycoside", formula="C4H8O4", kind="neutral_loss",
                       ion_mode="negative", adduct="[M-H]-", standards_based=True),
                stage2("s2-a-002", "p3a", family="flavone_C_glycoside", formula="C4H8O4 (0,2X)",
                       kind="neutral_loss", ion_mode="negative", adduct="[M−H]−", standards_based=True),
                stage2("s2-b-001", "p3b", family="flavonoid_C_glycoside", formula="C4H8O4", kind="neutral_loss",
                       ion_mode="negative", adduct="[M-H]-"),
            ],
            "stage3": [
                stage2("s3-c-001", "p3c", family="flavonoid_C_glycoside", value="loss of C4H8O4 (120.0423)",
                       kind="neutral_loss", ion_mode="negative", adduct="[M-H]-", standards_based=True,
                       instrument="Finnigan LCQ ion trap"),
                stage2("s3-c-002", "p3c", family="flavonoid_C_glycoside", formula="C4H8O4", kind="neutral_loss",
                       ion_mode="positive", adduct="[M+H]+"),
            ],
            "npr": [
                record("npr-a-001", family="flavone_C_glycoside", formula="C4H8O4", kind="neutral_loss",
                       ion_mode="negative", adduct="[M-H]-", standards_based=None,
                       source="NPR 2021 (doi:10.1039/D1NP00014D), p.1", primary_refs=["Ref 1, J Synth 2019"]),
                record("npr-b-002", family="flavone_c_glycoside", formula="C4H8O4", kind="neutral_loss",
                       ion_mode="negative", adduct="[M-H]-", standards_based=None, source="NPR 2021, p.2",
                       primary_refs=["Ref 2, J Synth 2018"]),
            ],
        }
        result = build(recs)
        neg = [g for g in result["groups"] if g["ion_mode"] == "negative"]
        self.assertEqual(len(neg), 1)
        g = neg[0]
        self.assertEqual((g["family"], g["kind"], g["adduct"], g["formula"]),
                         ("flavonoid_c_glycoside", "neutral_loss", ["[M-H]-"], "C4H8O4"))
        self.assertEqual(g["record_ids"], ["npr-a-001", "npr-b-002", "s2-a-001", "s2-a-002", "s2-b-001", "s3-c-001"])
        self.assertEqual(g["n_papers"], 4)  # p3a, p3b, p3c and the NPR review (npr-a and npr-b are one paper)
        self.assertEqual(g["n_standards_papers"], 2)  # p3a and p3c
        self.assertIn("doi:10.1039/d1np00014d", g["citations"])
        self.assertIn("doi:10.1000/p3a", g["citations"])
        self.assertEqual(g["secondary_citations"], ["Ref 1, J Synth 2019", "Ref 2, J Synth 2018"])
        self.assertEqual(g["sources"], ["npr", "stage2", "stage3"])
        self.assertEqual(g["instrument_classes"], {"ion_trap": 1, "qtof": 5})
        self.assertEqual(len(result["groups"]), 2)  # the positive-mode row is its own group
        self.assertTrue(all(t["group_id"] for t in result["terms"]))

    def test_counts_tables(self):
        result = build({"stage2": [stage2("s2-a-001", "p3a", formula="C9H7O3+", standards_based=True)],
                        "stage3": [stage2("s3-b-001", "p3b", formula="C9H7O3+"),
                                   stage2("s3-b-002", "p3b", kind="other")]})
        fields, rows = nc.counts_family_kind(result["records"], ["stage2", "stage3"])
        total = rows[-1]
        self.assertEqual((total["family"], total["stage2"], total["stage3"], total["total"]), ("TOTAL", 1, 2, 3))
        fields, rows = nc.counts_family_group(result["records"], result["terms"], result["groups"])
        fam = next(r for r in rows if r["family"] == "flavonoid_o_glycoside")
        self.assertEqual((fam["n_records"], fam["n_term_rows"], fam["n_groups"], fam["n_standards_papers"]),
                         (3, 2, 1, 1))


class OutputFolderTests(unittest.TestCase):
    def test_folder_name_rules(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp) / "ws8-p3"
            parent.mkdir()
            self.assertEqual(nc.check_out_dir(parent / "candidates-20261118"), parent / "candidates-20261118")
            for bad in (Path(tmp) / "elsewhere" / "candidates-20261118", parent / "candidates-2026111",
                        parent / "candidates-20261118x", parent / "run-20261118"):
                with self.assertRaises(nc.NormaliseError, msg=str(bad)):
                    nc.check_out_dir(bad)

    def test_refuses_non_empty_folder_and_accepts_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "ws8-p3" / "candidates-20261118"
            out.mkdir(parents=True)
            nc.check_out_dir(out)
            (out / "old.tsv").write_text("x", encoding="utf-8")
            with self.assertRaises(nc.NormaliseError):
                nc.check_out_dir(out)


class EndToEndTests(unittest.TestCase):
    def write_inputs(self, tmp):
        inputs = {
            "npr": [record("npr-a-001", formula="C9H7O3+", source="NPR 2021, p.1", primary_refs=["Ref 1"])],
            "stage2": [stage2("s2-a-001", "p3a", formula="C9H7O3+", value="m/z 163.0390"),
                       stage2("s2-a-002", "p3a", kind="veto", value="m/z 137 alone is not specific")],
            "stage3": [stage2("s3-b-001", "p3b", kind="other", value="a caveat"),
                       stage2("s3-b-002", "p3b", formula="C9H7O3+", value="m/z 163.0300")],
        }
        args = []
        for label, recs in inputs.items():
            path = Path(tmp) / f"{label}.json"
            path.write_text(json.dumps(recs), encoding="utf-8")
            args.append(f"{label}={path}")
        return inputs, args

    def test_main_writes_all_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            inputs, args = self.write_inputs(tmp)
            (Path(tmp) / "ws8-p3").mkdir()
            out = Path(tmp) / "ws8-p3" / "candidates-20261118"
            self.assertEqual(quiet_main(["--inputs", *args, "--out", str(out)]), 0)
            names = {p.name for p in out.iterdir()}
            self.assertEqual(names, {"records.tsv", "terms.tsv", "groups.json", "counts_family_kind.tsv",
                                     "counts_family_group.tsv", "manual_queue.tsv", "mass_flags.tsv",
                                     "unmapped.tsv", "run.json"})
            run = json.loads((out / "run.json").read_text(encoding="utf-8"))
            self.assertEqual(run["reconciliation"]["status"], "ok")
            self.assertEqual(run["counts"]["n_records_in"], 5)
            self.assertEqual(run["counts"]["n_records_out"], 5)
            self.assertEqual([i["label"] for i in run["inputs"]], ["npr", "stage2", "stage3"])
            self.assertTrue(all(re.fullmatch(r"[0-9a-f]{64}", i["sha256"]) for i in run["inputs"]))
            self.assertRegex(run["code"]["script_sha256"], r"^[0-9a-f]{64}$")
            with open(out / "records.tsv", encoding="utf-8", newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(sorted(r["id"] for r in rows), sorted(r["id"] for v in inputs.values() for r in v))
            with open(out / "manual_queue.tsv", encoding="utf-8", newline="") as f:
                self.assertEqual([r["id"] for r in csv.DictReader(f, delimiter="\t")], ["s3-b-001"])
            with open(out / "mass_flags.tsv", encoding="utf-8", newline="") as f:
                self.assertEqual([r["record_id"] for r in csv.DictReader(f, delimiter="\t")], ["s3-b-002"])
            self.assertEqual(quiet_main(["--inputs", *args, "--out", str(out)]), 2)  # now non-empty

    def test_reconciliation_failure_is_a_hard_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            inputs, args = self.write_inputs(tmp)
            loaded = nc.load_inputs(args)
            result = nc.build(loaded, nc.load_maps())
            with self.assertRaises(nc.ReconciliationError):
                nc.reconcile(loaded, result["records"][:-1])
            with self.assertRaises(nc.ReconciliationError):
                nc.reconcile(loaded, result["records"] + result["records"][:1])
            extra = dict(result["records"][0], id="not-an-input")
            with self.assertRaises(nc.ReconciliationError):
                nc.reconcile(loaded, result["records"] + [extra])
            (Path(tmp) / "ws8-p3").mkdir()
            out = Path(tmp) / "ws8-p3" / "candidates-20261118"
            broken = dict(result, records=result["records"][:-1])
            with self.assertRaises(nc.ReconciliationError):
                nc.write_outputs(out, loaded, broken, [])
            self.assertFalse((out / "run.json").exists())

    def test_repeated_input_id_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "stage2.json"
            path.write_text(json.dumps([stage2("s2-a-001", "p3a"), stage2("s2-a-001", "p3a")]), encoding="utf-8")
            loaded = nc.load_inputs([f"stage2={path}"])
            result = nc.build(loaded, nc.load_maps())
            with self.assertRaises(nc.ReconciliationError):
                nc.reconcile(loaded, result["records"])

    def test_bad_input_label_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.json"
            path.write_text("[]", encoding="utf-8")
            (Path(tmp) / "ws8-p3").mkdir()
            out = Path(tmp) / "ws8-p3" / "candidates-20261118"
            self.assertEqual(quiet_main(["--inputs", f"stage9={path}", "--out", str(out)]), 2)


class PublicRepositoryTests(unittest.TestCase):
    def test_no_local_paths_in_tool_files(self):
        for path in TOOL_DIR.iterdir():
            if path.suffix in (".py", ".tsv", ".md"):
                text = path.read_text(encoding="utf-8")
                self.assertIsNone(re.search(r"\b[A-Za-z]:[\\/]", text), path.name)
                self.assertNotIn("Users", text, path.name)


if __name__ == "__main__":
    unittest.main()
