"""WS8 P3 candidate normaliser, on synthetic records only (never the literature files)."""
from collections import defaultdict
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
        # its tokens stay in terms.tsv but are not grouped until a reader assigns a kind
        (term,) = result["terms"]
        self.assertEqual((term["group_exclusion"], term["group_id"]), (["manual"], ""))
        self.assertIn("manual", term["flags"])
        self.assertEqual(result["groups"], [])

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

    def test_record_override_moves_di_c_records(self):
        overrides = nc.load_family_override()
        self.assertTrue(overrides)
        unified = {v["family"] for v in self.fam.values()}
        for (source, rid), o in overrides.items():
            self.assertIn(o["family"], unified)
            if o["family_original"] == "flavonoid_C_glycoside":
                self.assertEqual((o["family"], o["family_group"]),
                                 ("flavonoid_di_c_glycoside", "flavonoid_c_glycoside"))
            else:
                # a steryl_glycoside record about steryl esters only: an acyl-lipid class, out of scope
                self.assertEqual((o["family_original"], o["family"], o["scope"]),
                                 ("steryl_glycoside", "steryl_ester", "out_of_scope"))
        self.assertEqual(self.fam["steryl_glycoside"]["scope"], "in_scope")
        c_glycoside = next(o for o in overrides.values() if o["family_original"] == "flavonoid_C_glycoside")
        maps = nc.load_maps()
        maps["family_override"] = {("stage3", "s3-x-001"): dict(c_glycoside)}
        recs = {"stage3": [stage2("s3-x-001", "p3x", family="flavonoid_C_glycoside"),
                           stage2("s3-x-002", "p3x", family="flavonoid_C_glycoside")]}
        result = nc.build({k: (Path(f"{k}.json"), v) for k, v in recs.items()}, maps)
        rows = {r["id"]: r for r in result["records"]}
        self.assertEqual((rows["s3-x-001"]["family"], rows["s3-x-001"]["family_original"]),
                         ("flavonoid_di_c_glycoside", "flavonoid_C_glycoside"))
        self.assertTrue(rows["s3-x-001"]["family_override"])
        self.assertEqual((rows["s3-x-002"]["family"], rows["s3-x-002"]["family_override"]),
                         ("flavonoid_c_glycoside", ""))
        recs["stage3"][0]["family"] = "flavonol"
        with self.assertRaises(nc.NormaliseError):
            nc.build({k: (Path(f"{k}.json"), v) for k, v in recs.items()}, maps)

    def test_override_file_is_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "override.tsv"
            path.write_text("source\tid\tfamily_original\tfamily\treason\n"
                            "stage3\ts3-x-001\tflavonoid_C_glycoside\tnot_a_family\tx\n", encoding="utf-8")
            with self.assertRaises(nc.NormaliseError):
                nc.load_family_override(path)

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
        self.assertEqual(nc.parse_adducts("[M+H]+.")[0], ["[M+H]+•"])  # a radical mark after the charge
        self.assertEqual(nc.parse_adducts("[M+H]+. (radical cation)")[0], ["[M+H]+•"])
        self.assertEqual(nc.parse_adducts("[M+H]+.5")[0], ["[M+H]+"])
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
        self.assertEqual(terms["s2-x-001#t1.pos"]["adduct"], "[M+H]+")
        self.assertEqual(terms["s2-x-001#t1.neg"]["adduct"], "[M-H]-")
        self.assertEqual(terms["s2-x-001#t1.neg"]["ion_mode"], "negative")
        self.assertTrue(all(t["expanded_from_both"] for t in terms.values()))
        self.assertFalse(any(t["expanded_from_adducts"] for t in terms.values()))
        self.assertIn("kind_conflict", terms["s2-x-001#t2.pos"]["flags"])
        self.assertEqual(result["records"][0]["ion_mode"], "both")

    def test_both_takes_the_polarity_the_text_gives_a_neutral_token(self):
        value = ("fragments common to all isomers: C14H9O2 and C13H9O2 (negative), C13H9 and C12H9 (positive); "
                 "every isomer differs; C12H9 occurs for all of them in both modes")
        result = build({"stage3": [stage2("s3-x-001", "p3x", ion_mode="both", adduct="[M-H]- / [M+H]+",
                                          value=value)]})
        modes = defaultdict(set)
        for t in result["terms"]:
            modes[t["formula"]].add(t["ion_mode"])
        self.assertEqual(dict(modes), {"C14H9O2": {"negative"}, "C13H9O2": {"negative"}, "C13H9": {"positive"},
                                       "C12H9": {"positive", "negative"}})
        terms = {t["term_id"]: t for t in result["terms"]}
        self.assertEqual(set(terms), {"s3-x-001#t1.neg", "s3-x-001#t2.neg", "s3-x-001#t3.pos", "s3-x-001#t4.pos",
                                      "s3-x-001#t4.neg"})
        self.assertIn("polarity_from_text", terms["s3-x-001#t1.neg"]["flags"])
        self.assertEqual(terms["s3-x-001#t4.neg"]["polarity_label"], "both")
        self.assertNotIn("polarity_unverified", terms["s3-x-001#t4.neg"]["flags"])
        # a heading of the clause gives the polarity; a polarity word elsewhere leaves the token unverified
        result = build({"stage3": [
            stage2("s3-x-002", "p3x", ion_mode="both", kind="neutral_loss", adduct="[M+H]+ / [M-H]-",
                   value="positive: losses of C2H2O (ketene; 245 -> 203), then CH3 to 188; negative: CO2 loss"),
            stage2("s3-x-003", "p3x", ion_mode="both", kind="neutral_loss", adduct="[M+H]+ / [M-H]-",
                   value="radical loss of CH3O (31) is the base peak: 477->446 (neg); 523->492 (pos)"),
            stage2("s3-x-004", "p3x", ion_mode="both", kind="neutral_loss", adduct="[M+H]+ / [M-H]-",
                   formula="C6H10O5")]})
        modes = defaultdict(set)
        flags = defaultdict(set)
        for t in result["terms"]:
            modes[(t["record_id"], t["formula"])].add(t["ion_mode"])
            flags[(t["record_id"], t["formula"])].update(t["flags"])
        self.assertEqual(modes[("s3-x-002", "C2H2O")], {"positive"})
        self.assertEqual(modes[("s3-x-002", "CH3")], {"positive"})
        self.assertEqual(modes[("s3-x-003", "CH3O")], {"positive", "negative"})
        self.assertIn("polarity_unverified", flags[("s3-x-003", "CH3O")])
        self.assertEqual(modes[("s3-x-004", "C6H10O5")], {"positive", "negative"})
        self.assertNotIn("polarity_unverified", flags[("s3-x-004", "C6H10O5")])

    def test_both_without_adduct_of_a_polarity(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", ion_mode="both", adduct="[M+H]+", formula="H2O")]})
        terms = {t["ion_mode"]: t for t in result["terms"]}
        self.assertEqual(terms["negative"]["adduct"], "unstated")

    def test_one_term_row_per_adduct(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", ion_mode="negative", adduct="[M-H]-/[M-2H]2-",
                                          formula="O3S-; C6H10O5"),
                                   stage2("s2-x-002", "p3x", ion_mode="both", adduct="[M+H]+, [M+Na]+, [M-H]-",
                                          kind="neutral_loss", formula="C6H10O5")]})
        terms = {t["term_id"]: t for t in result["terms"]}
        self.assertEqual(set(terms), {"s2-x-001#t1.a1", "s2-x-001#t1.a2", "s2-x-001#t2.a1", "s2-x-001#t2.a2",
                                      "s2-x-002#t1.pos.a1", "s2-x-002#t1.pos.a2", "s2-x-002#t1.neg"})
        self.assertEqual([terms[f"s2-x-001#t1.a{j}"]["adduct"] for j in (1, 2)], ["[M-H]-", "[M-2H]2-"])
        self.assertEqual([terms[f"s2-x-002#t1.pos.a{j}"]["adduct"] for j in (1, 2)], ["[M+H]+", "[M+Na]+"])
        self.assertTrue(terms["s2-x-001#t1.a2"]["expanded_from_adducts"])
        self.assertFalse(terms["s2-x-002#t1.neg"]["expanded_from_adducts"])
        self.assertEqual(result["records"][0]["n_term_rows"], 4)
        self.assertEqual(result["records"][0]["n_formula_tokens"], 2)


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

    def test_star_and_dot_radical_marks_keep_the_charge(self):
        cases = {"C22H28N2O3*-": ("C22H28N2O3-", -1), "C19H15NO4*+": ("C19H15NO4+", 1),
                 "C23H28O5.+": ("C23H28O5+", 1), "SO3.- (79.9574)": ("O3S-", -1), "C7H6O.-": ("C7H6O-", -1),
                 "C8H15N-*": ("C8H15N-", -1), "C7H4O3.-, C8H7O3-": ("C7H4O3-", -1)}
        for text, (formula, charge) in cases.items():
            tok = self.tokens(text)[0]
            self.assertEqual((tok["formula"], tok["charge"], tok["radical_printed"]), (formula, charge, True), text)
        self.assertEqual([t["formula"] for t in self.tokens("C7H4O3.-, C8H7O3-")], ["C7H4O3-", "C8H7O3-"])
        tok = self.tokens("CH3* (I = C22H28N2O3*-)")[0]
        self.assertEqual((tok["formula"], tok["charge"], tok["radical_printed"]), ("CH3", 0, True))
        self.assertEqual(self.tokens("C9H7O3*+")[0]["radical_mark_mismatch"], True)
        # a full stop after an uncharged formula ends a sentence and is no radical mark
        tok = self.tokens("loss of C6H10O5. Then")[0]
        self.assertEqual((tok["formula"], tok["radical_printed"]), ("C6H10O5", False))
        result = build({"stage2": [stage2("s2-x-001", "p3x", kind="neutral_loss", ion_mode="negative",
                                          adduct="[M-H]-", formula="CH3* (I = C22H28N2O3*-)")]})
        flags = {t["formula"]: t["flags"] for t in result["terms"]}
        self.assertIn("kind_conflict", flags["C22H28N2O3-"])
        self.assertNotIn("kind_conflict", flags["CH3"])

    def test_unsupported_element(self):
        tok = self.tokens("C6H10O5Li+")[0]
        self.assertEqual(tok["unsupported_elements"], ["Li"])
        self.assertIsNone(tok["mass"])
        result = build({"stage3": [stage2("s3-x-001", "p3x", formula="C6H10O5Li+")]})
        self.assertIn("unsupported_element", result["terms"][0]["flags"])
        self.assertEqual(result["terms"][0]["mass_check"], "no_mass")

    def test_words_notations_and_ion_notation_contents_are_not_terms(self):
        self.assertEqual([t["formula"] for t in self.tokens("K C18H11O4+")], ["C18H11O4+"])
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

    def test_composites_and_fragments_are_logged(self):
        tokens, rejected, composites, _ = nc.tokenise_formulas("computed CH2O+2H2O 66.0317; C3H6O3 90.0317")
        self.assertEqual([(t["formula"], t["composite_part"]) for t in tokens],
                         [("CH2O", True), ("C3H6O3", False)])
        self.assertEqual((rejected, composites), (["2H2O"], ["CH2O+2H2O"]))
        tokens, rejected, composites, _ = nc.tokenise_formulas("RCOOH + NH3")
        self.assertEqual(([t["formula"] for t in tokens], rejected, composites), (["H3N"], ["RCOOH"],
                                                                                 ["RCOOH + NH3"]))
        # a formula that also stands alone is no composite part; a charge is no '+' operator; units are no operand
        tokens, _, composites, _ = nc.tokenise_formulas("C8H16O2 (144.1150); C8H20O4 = C8H16O2 + 2H2O")
        self.assertEqual(composites, ["C8H16O2 + 2H2O"])
        self.assertEqual(nc.tokenise_formulas("C9H7O3+ C7H5O2+")[2], [])
        self.assertEqual(nc.tokenise_formulas("72 Da + C2H4 losses")[2], [])
        self.assertEqual(nc.tokenise_formulas("0,2X C17H12O6; 1,3B+ C9H9O+")[1], [])
        result = build({"stage2": [stage2("s2-x-001", "p3x", kind="neutral_loss",
                                          formula="C8H16O2 (144.1150); C8H20O4 = C8H16O2 + 2H2O"),
                                   stage2("s2-x-002", "p3x", kind="neutral_loss", formula="CH2O+2H2O; C3H6O3")]})
        self.assertIn(("formula_composite", "C8H16O2 + 2H2O"), result["unmapped"])
        self.assertEqual(sorted(result["unmapped"][("formula_fragment", "2H2O")]), ["s2-x-001", "s2-x-002"])
        rows = {r["id"]: r for r in result["records"]}
        self.assertEqual(rows["s2-x-002"]["formula_flags"], ["fragment_not_a_term", "composite_expression"])
        terms = {(t["record_id"], t["formula"]): t for t in result["terms"]}
        self.assertEqual(terms[("s2-x-001", "C8H16O2")]["group_exclusion"], [])
        self.assertEqual(terms[("s2-x-002", "CH2O")]["group_exclusion"], ["composite_part"])
        self.assertEqual(terms[("s2-x-002", "CH2O")]["group_id"], "")
        self.assertTrue(terms[("s2-x-002", "C3H6O3")]["group_id"])

    def test_precursor_role(self):
        def roles(text, **kw):
            return [(t["formula"], t["role"]) for t in self.tokens(text, **kw)]
        self.assertEqual(roles("C9H7O4- (179); C9H5O3- (161); precursor C16H15O8-"),
                         [("C9H7O4-", ""), ("C9H5O3-", ""), ("C16H15O8-", "precursor")])
        self.assertEqual(roles("C14H9O3- (from precursor C14H11O4-)"),
                         [("C14H9O3-", ""), ("C14H11O4-", "precursor")])
        self.assertEqual(roles("C27H44NO3+ (precursor); 269 C19H25O+"),
                         [("C27H44NO3+", "precursor"), ("C19H25O+", "")])
        self.assertEqual(roles("C34H29O15- ([M-H]-, theor. 677.1512)"), [("C34H29O15-", "precursor")])
        self.assertEqual(roles("C15H9O6- ([M-H-CO2]-)"), [("C15H9O6-", "")])
        self.assertEqual(roles("precursors m/z 556.363 (C33H50NO6+ 556.3633), 572.358 (C33H50NO7+ 572.3582)"),
                         [("C33H50NO6+", "precursor"), ("C33H50NO7+", "precursor")])
        self.assertEqual(roles("M+ 1609.4442 (C78H81O37+) -> 1271 (C60H63O28+)", from_value=True),
                         [("C78H81O37+", "precursor"), ("C60H63O28+", "")])
        self.assertEqual(roles("[M-H]- 963 (C43H48O25 neutral) -> 801", from_value=True),
                         [("C43H48O25", "precursor")])
        self.assertEqual(roles("not printed; MH+ C10H16N3O5+ 258.1084, BH2+ C5H8N3O+"),
                         [("C10H16N3O5+", "precursor"), ("C5H8N3O+", "")])
        result = build({"stage2": [stage2("s2-x-001", "p3x", kind="neutral_loss", ion_mode="negative",
                                          adduct="[M-H]-", formula="C6H10O5; [M-H]- C21H19O11- 447.0933")]})
        terms = {t["formula"]: t for t in result["terms"]}
        self.assertEqual(terms["C21H19O11-"]["group_exclusion"], ["precursor"])
        self.assertNotIn("kind_conflict", terms["C21H19O11-"]["flags"])
        self.assertEqual(terms["C21H19O11-"]["group_id"], "")
        self.assertEqual([g["formula"] for g in result["groups"]], ["C6H10O5"])

    def test_precursor_heading_and_notation_before_at(self):
        def roles(text, **kw):
            return [(t["formula"], t["role"], t["role_evidence"]) for t in self.tokens(text, **kw)]
        text = ("loss C6H10O5 (162.0528, recomputed; not printed). Precursor formulas from Table 1: C25H32O13, "
                "C25H32O12, C24H30O13")
        self.assertEqual(roles(text), [("C6H10O5", "", ""), ("C25H32O13", "precursor", "clause"),
                                       ("C25H32O12", "precursor", "clause"), ("C24H30O13", "precursor", "clause")])
        # a sentence end closes the heading's clause; 'precursor m/z 353:' is no heading
        self.assertEqual(roles("Precursor ions: C16H17O9-. Products are C9H7O4- and C7H11O6-")[1:],
                         [("C9H7O4-", "", ""), ("C7H11O6-", "", "")])
        self.assertEqual(roles("MS3 of precursor m/z 353: C7H11O6- 191"), [("C7H11O6-", "", "")])
        self.assertEqual(roles("observed as a radical cation [M+H]+. at m/z 369.1686 (C18H27NO7)", from_value=True),
                         [("C18H27NO7", "precursor", "before")])
        self.assertEqual(roles("[M-H]- ion at m/z 447.0933 (C21H19O11-)"), [("C21H19O11-", "precursor", "before")])
        self.assertEqual(roles("loss of CO2 from [M-H]- is C14H9O3-"), [("CO2", "", ""), ("C14H9O3-", "", "")])

    def test_formulas_named_as_not_observed_or_quoted_as_printed(self):
        def roles(text, **kw):
            return [(t["formula"], t["role"]) for t in self.tokens(text, **kw)]
        self.assertEqual(roles("recomputed: C6H3O5- 154.9986; C5H3O3- 111.0088 (the QA-path analogues would be "
                               "C7H7O4- 155.0350 and C6H7O2- 111.0452)"),
                         [("C6H3O5-", ""), ("C5H3O3-", ""), ("C7H7O4-", "not_observed"), ("C6H7O2-", "not_observed")])
        self.assertEqual(roles("(would be C9H7O3+); C9H5O2+"), [("C9H7O3+", "not_observed"), ("C9H5O2+", "")])
        self.assertEqual(roles("C7H11O6- (printed as C7H12O6); C7H9O5-; C16H17O9-"),
                         [("C7H11O6-", ""), ("C7H12O6", "printed_as"), ("C7H9O5-", ""), ("C16H17O9-", "")])
        self.assertEqual(roles("58 u, printed as CH2COOH; C2H2O2 = 58.0055 by calculation"),
                         [("C2H3O2", "printed_as"), ("C2H2O2", "")])
        self.assertEqual(roles("C3H6O3 (90.0317); fragments printed as C18H13O8 (357), C17H11O7 (327)"),
                         [("C3H6O3", ""), ("C18H13O8", "printed_as"), ("C17H11O7", "printed_as")])
        self.assertEqual(roles("(printed as C7H12O6), C7H9O5-"), [("C7H12O6", "printed_as"), ("C7H9O5-", "")])
        text = ("C6H10O2 (114.0681); 297 printed as C21H27O+, but m/z 297 fits C21H29O+ (297.2213; C21H27O+ would be "
                "295.2056)")
        self.assertEqual(roles(text), [("C6H10O2", ""), ("C21H27O+", "printed_as"), ("C21H29O+", ""),
                                       ("C21H27O+", "not_observed")])
        result = build({"stage2": [
            stage2("s2-x-001", "p3x", ion_mode="negative", adduct="[M-H]-",
                   formula="C6H3O5- 154.9986 (the analogues would be C7H7O4- 155.0350)"),
            stage2("s2-x-002", "p3y", ion_mode="negative", adduct="[M-H]-", formula="C7H7O4-"),
            stage2("s2-x-003", "p3z", kind="neutral_loss", formula=text, adduct="[M+H]+")]})
        terms = {(t["record_id"], t["formula"]): t for t in result["terms"]}
        self.assertEqual(terms[("s2-x-001", "C7H7O4-")]["group_exclusion"], ["not_observed"])
        self.assertEqual(terms[("s2-x-001", "C7H7O4-")]["group_id"], "")
        self.assertEqual(terms[("s2-x-003", "C21H27O+")]["group_exclusion"], ["printed_as"])  # both occurrences
        self.assertEqual(terms[("s2-x-003", "C21H29O+")]["group_exclusion"], [])
        group = next(g for g in result["groups"] if g["formula"] == "C7H7O4-")
        self.assertEqual((group["record_ids"], group["n_papers"]), (["s2-x-002"], 1))
        self.assertEqual(result["unmapped"][("formula_not_observed", "C7H7O4-")], ["s2-x-001"])
        self.assertEqual(result["unmapped"][("formula_printed_as", "C21H27O+")], ["s2-x-003"])

    def test_isotope_labels_are_not_terms(self):
        tokens, rejected, composites, isotopes = nc.tokenise_formulas(
            "C10H10N+ (144.0808); C10H10(15N)+ in labelled plants")
        self.assertEqual([t["formula"] for t in tokens], ["C10H10N+"])
        self.assertEqual(isotopes, ["C10H10(15N)+"])
        self.assertEqual(nc.tokenise_formulas("(13C)C5H10O5 and C6H10O5 (162)")[3], ["(13C)C5H10O5"])
        self.assertEqual(nc.tokenise_formulas("C6H10O5 (162.0528)")[3], [])
        result = build({"stage3": [stage2("s3-x-001", "p3x", formula="C10H10N+ (144.0808); C10H10(15N)+ labelled")]})
        self.assertEqual([t["formula"] for t in result["terms"]], ["C10H10N+"])
        self.assertEqual(result["unmapped"][("formula_isotope_label", "C10H10(15N)+")], ["s3-x-001"])
        self.assertIn("isotope_label_not_a_term", result["records"][0]["formula_flags"])

    def test_ion_before_an_arrow_is_the_precursor_of_a_loss(self):
        text = "C24H39O5- 407.2803 -> C23H35O2- 343.2643 (calculated)"
        result = build({"stage2": [
            stage2("s2-x-001", "p3x", kind="neutral_loss", ion_mode="negative", adduct="[M-H]-", formula=text),
            stage2("s2-x-002", "p3x", kind="fragment_ion", ion_mode="negative", adduct="[M-H]-", formula=text)]})
        terms = {(t["record_id"], t["formula"]): t for t in result["terms"]}
        self.assertEqual((terms[("s2-x-001", "C24H39O5-")]["role"], terms[("s2-x-001", "C24H39O5-")]["role_evidence"]),
                         ("precursor", "arrow"))
        self.assertNotIn("kind_conflict", terms[("s2-x-001", "C24H39O5-")]["flags"])
        self.assertIn("kind_conflict", terms[("s2-x-001", "C23H35O2-")]["flags"])
        # in a product-ion record the first ion of a pathway may itself be the diagnostic product
        self.assertEqual(terms[("s2-x-002", "C24H39O5-")]["role"], "")

    def test_term_override(self):
        maps = nc.load_maps()
        maps["term_override"] = {("stage2", "s2-x-001", "C16H17O9-"): {"role": "precursor", "reason": "test"}}
        recs = {"stage2": [stage2("s2-x-001", "p3x", ion_mode="negative", adduct="[M-H]-",
                                  formula="C7H11O6-; C16H17O9-")]}
        result = nc.build({k: (Path(f"{k}.json"), v) for k, v in recs.items()}, maps)
        terms = {t["formula"]: t for t in result["terms"]}
        self.assertEqual((terms["C16H17O9-"]["role"], terms["C16H17O9-"]["role_evidence"]), ("precursor", "override"))
        self.assertEqual(terms["C16H17O9-"]["group_id"], "")
        self.assertEqual(result["overrides"]["term_override_rows_used"], 1)
        maps["term_override"] = {("stage2", "s2-x-001", "C9H7O4-"): {"role": "precursor", "reason": "stale"}}
        with self.assertRaises(nc.NormaliseError):
            nc.build({k: (Path(f"{k}.json"), v) for k, v in recs.items()}, maps)
        self.assertTrue(nc.load_term_override())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "term.tsv"
            path.write_text("source\tid\tformula\trole\treason\nstage2\ts2-x-001\tC9H7O4-\tproduct\tx\n",
                            encoding="utf-8")
            with self.assertRaises(nc.NormaliseError):
                nc.load_term_override(path)

    def test_negative_evidence_terms_are_not_grouped(self):
        result = build({"stage3": [stage2("s3-x-001", "p3x", kind="diagnostic_absence", ion_mode="positive",
                                          formula="C8H12NO2+ (154.0863); C8H10NO+"),
                                   stage2("s3-x-002", "p3x", formula="C8H12NO2+")],
                        "stage2": [stage2("s2-x-001", "p3y", kind="veto", formula="C8H12NO2+")]})
        for t in result["terms"]:
            if t["record_id"] in ("s3-x-001", "s2-x-001"):
                self.assertEqual(t["group_exclusion"], ["role_unverified"])
                self.assertIn("role_unverified", t["flags"])
                self.assertEqual(t["group_id"], "")
        self.assertEqual([(g["kind"], g["record_ids"]) for g in result["groups"]], [("product_ion", ["s3-x-002"])])

    def test_odd_electron_neutral_loss(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", kind="neutral_loss",
                                          formula="CH3; C27H29O16; C6H10O5")]})
        flags = {t["formula"]: t["flags"] for t in result["terms"]}
        self.assertIn("odd_electron_neutral", flags["CH3"])
        self.assertIn("odd_electron_neutral", flags["C27H29O16"])
        self.assertNotIn("odd_electron_neutral", flags["C6H10O5"])

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

    def check(self, mass, *printed, note=(), hill=None):
        return nc.mass_check(mass, nc.printed_numbers(" ".join(printed)), nc.printed_numbers(" ".join(note)), hill)

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

    def test_note_mismatch_is_advisory(self):
        m = self.IONS["C9H7O3+"]
        result = self.check(m, "m/z 163 (nominal)", note=["recomputed 163.0300"])
        self.assertEqual((result["status"], result["printed_from"]), ("note_mismatch", "note"))
        self.assertEqual(self.check(m, "m/z 163", note=["about 163.1"])["status"], "note_lowres")
        result = build({"stage2": [stage2("s2-x-001", "p3x", formula="C9H7O3+", value="m/z 163",
                                          note="recomputed 163.0300")]})
        self.assertIn("note_mismatch", result["terms"][0]["flags"])
        self.assertNotIn("printed_mismatch", result["terms"][0]["flags"])

    def test_number_next_to_another_formula_is_not_compared(self):
        sodiated = self.IONS["C18H30NaO14+"]  # 493.1528; the aglycone ion C30H46NaO4+ is 493.3288
        note = ["493 is also the sodiated aglycone (computed C30H46NaO4+ 493.3288)"]
        hill = "C18H30NaO14"
        self.assertEqual(self.check(sodiated, "m/z 493", note=note, hill=hill)["status"], "no_printed_value")
        self.assertEqual(self.check(sodiated, "m/z 493", note=note)["status"], "note_mismatch")  # no formula given
        self.assertEqual(self.check(sodiated, "m/z 493.3288 (C30H46NaO4+)", hill=hill)["status"], "no_printed_value")
        self.assertEqual(self.check(sodiated, "C18H30NaO14+ 493.1600", hill=hill)["status"], "printed_mismatch")
        self.assertEqual(self.check(sodiated, "493.1600", hill=hill)["status"], "printed_mismatch")
        result = build({"stage2": [stage2("s2-x-001", "p3x", formula="C18H30NaO14+", value="m/z 493 (nominal)",
                                          adduct="[M+Na]+", note=note[0])]})
        self.assertEqual(result["terms"][0]["mass_check"], "no_printed_value")
        numbers = nc.printed_numbers("C9H7O3+ 163.0390; 145.0284 (C9H5O2+), 117.0335")
        self.assertEqual([(n[2], n[3]) for n in numbers],
                         [("163.0390", "C9H7O3"), ("145.0284", "C9H5O2"), ("117.0335", None)])

    def test_pairing_does_not_hide_a_match(self):
        # a number next to a formula is that formula's mass only when it lies near it; one-element tokens, ':' and
        # composite parts never pair
        cases = [
            ("C21H31+", "I(m/z 283.2420) > I(m/z 289.2162) for the 12-OH aglycone", "matched"),
            ("C27H45+", "sterol cation after loss of hexose + NH3: 369.3505 (cholesterol)", "matched"),
            ("C6H15NO6", "NL 197.09 (C6H12O6 + NH3)", "matched"),
            ("C28H48NO2-", "474.3589 -> 430.3670 (CO2 loss, as printed)", "matched"),
            ("C13H18NO2+", "losses: C5H8O4 (220.1342, C13H18NO2+)", "matched"),
            ("C13H18NO2+", "losses: C5H8O4 (220.1442, C13H18NO2+)", "printed_mismatch"),
        ]
        for text, value, status in cases:
            tok = nc.tokenise_formulas(text)[0][0]
            got = nc.mass_check(tok["mass"], nc.printed_numbers(value), (), tok["hill"])["status"]
            self.assertEqual(got, status, f"{text} | {value}")
        self.assertEqual([n[3] for n in nc.printed_numbers("I(m/z 283.2420)")], [None])
        self.assertEqual([n[3] for n in nc.printed_numbers("C9H7O3+: 163.0390")], [None])
        self.assertEqual([n[3] for n in nc.printed_numbers("C5H8O4 (132.0423)")], ["C5H8O4"])

    def test_sentence_final_numbers_are_read(self):
        self.assertEqual([n[2] for n in nc.printed_numbers("ions at 303.0510, 305.0667.")],
                         ["303.0510", "305.0667"])
        self.assertEqual([n[2] for n in nc.printed_numbers("1.2.3 and 163.0390.5")], [])
        m = 305.066677  # C15H13O7-
        self.assertAlmostEqual(nc.tokenise_formulas("C15H13O7-")[0][0]["mass"], m, delta=1e-5)
        self.assertEqual(self.check(m, "no number", note=["computed 303.0510, 305.0667."])["status"], "matched")

    def test_mass_flags_reach_terms(self):
        result = build({"stage2": [stage2("s2-x-001", "p3x", formula="C9H7O3+", value="m/z 163.0300")]})
        self.assertEqual(result["terms"][0]["mass_check"], "printed_mismatch")
        self.assertIn("printed_mismatch", result["terms"][0]["flags"])


class InstrumentTests(unittest.TestCase):
    def setUp(self):
        self.rules = nc.load_instrument_map()

    def classify(self, text, collision=None, condition=None):
        return nc.classify_instrument(text, self.rules, collision, condition)

    def test_single_classes(self):
        self.assertEqual(self.classify("Agilent 6530 Q-TOF")[:2], ("qtof", "qtof"))
        self.assertEqual(self.classify("Thermo Q Exactive Orbitrap")[:2], ("orbitrap_ft", "orbitrap_ft"))
        self.assertEqual(self.classify("Thermo Q Exactive HF", "CID 35 eV")[:2], ("orbitrap_ft", "orbitrap_ft"))
        self.assertEqual(self.classify("Orbitrap Exploris 240", "stepped NCE 20/40/60")[:2],
                         ("orbitrap_ft", "orbitrap_ft"))
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

    def test_ion_trap_hybrids_take_the_stratum_of_their_activation(self):
        cases = [
            # (instrument, collision, condition, class, stratum)
            ("Thermo LTQ-Orbitrap XL", "ion trap CID 35-40%", None, "mixed", "ion_trap"),
            ("Thermo LTQ Orbitrap", "linear ion trap CID, 35%", None, "mixed", "ion_trap"),
            ("Thermo LTQ-FTICR 7 T", "LIT CID MSn", None, "mixed", "ion_trap"),
            ("Thermo LTQ FT Ultra", "CID 35%", None, "mixed", "ion_trap"),
            ("LTQ Orbitrap XL", "NCE 35%", None, "mixed", "ion_trap"),
            ("LTQ FT Ultra (MSn)", "35%", None, "mixed", "ion_trap"),
            ("Orbitrap Elite", "CID ramp 0-50 eV", None, "mixed", "ion_trap"),
            ("Thermo LTQ Orbitrap Velos Pro (HESI), CID in the linear ion trap", "NCE 25%", None, "mixed", "ion_trap"),
            ("Thermo LTQ-Orbitrap XL", "not stated (ion trap MS/MS)", None, "mixed", "ion_trap"),
            ("Orbitrap Elite (ion trap-Orbitrap)", "HCD 30, 70, 110 NCE", None, "mixed", "orbitrap_ft"),
            ("Thermo LTQ Orbitrap XL", "HCD", None, "mixed", "orbitrap_ft"),
            ("Orbitrap Fusion Lumos", "HCD CE 80", None, "mixed", "orbitrap_ft"),
            ("LIT-FTICR", "IRMPD", None, "mixed", "orbitrap_ft"),
            ("Thermo LTQ-FTICR 7 T", "IRMPD (CO2 laser, FTICR), not CID", "activation outside CID/HCD scope",
             "mixed", "orbitrap_ft"),
            ("Thermo LTQ Orbitrap XL", "HCD RCE 60% (Orbitrap) and LIT CID MSn", None, "mixed", "mixed"),
            ("LTQ Orbitrap XL", "CID 35% NCE; HCD for accurate mass", None, "mixed", "mixed"),
            ("Thermo LTQ-FTICR 7 T", "LIT CID MS2-MS5 (accurate masses from IRMPD-FTICR)", None, "mixed", "mixed"),
            ("Orbitrap Fusion Lumos (HESI)", "conventional CID or low-energy HCD", None, "mixed", "mixed"),
            ("Thermo LTQ-Orbitrap XL", "source/pre-trap SID 70%", None, "mixed", "mixed"),
            ("Thermo LTQ-Orbitrap XL", None, None, "mixed", "mixed"),
            ("Orbitrap Fusion Lumos", "abundant at low NCE, gone at high NCE", None, "mixed", "mixed"),
            ("Thermo LTQ Orbitrap XL", None, "counted among the top-10 HCD peaks", "mixed", "orbitrap_ft"),
            ("Thermo LTQ FT Ultra", "CID 35%; in-source", None, "mixed", "ion_trap"),
            ("XCT Ultra 6340 ion trap; LTQ Orbitrap for accurate mass", "ion trap MS2-MS4", None, "mixed", "ion_trap"),
        ]
        for inst, collision, condition, klass, strat in cases:
            got = self.classify(inst, collision, condition)
            self.assertEqual(got[:2], (klass, strat), f"{inst} | {collision} | {got}")
            self.assertIn(got[2], ("mixed.ltq_hybrid", "mixed.tribrid"))
        self.assertEqual(self.classify("Thermo LTQ FT Ultra", "CID 35%")[3], "trap (collision)")
        self.assertEqual(self.classify("Orbitrap Elite", "CID 25 eV; HCD 35")[3], "beam+trap (collision)")
        self.assertEqual(self.classify("Thermo LTQ-Orbitrap XL")[3], "none")
        self.assertEqual(self.classify("Agilent 6530 Q-TOF", "CID 20 eV")[3], "")
        # a standalone linear ion trap stays an ion trap
        self.assertEqual(self.classify("Thermo LTQ XL linear ion trap (ESI)", "HCD")[:2], ("ion_trap", "ion_trap"))
        self.assertEqual(self.classify("Thermo Velos Pro linear ion trap (MSn); accurate mass on Q Exactive Focus")[:2],
                         ("mixed", "ion_trap"))

    def test_tribrid_names_without_orbitrap(self):
        self.assertEqual(self.classify("Orbitrap Exploris 120/240, Ascend")[:3], ("mixed", "mixed", "mixed.tribrid"))
        self.assertEqual(self.classify("Orbitrap Exploris 480, Ascend", "HCD NCE 30")[:3],
                         ("mixed", "orbitrap_ft", "mixed.tribrid"))
        self.assertEqual(self.classify("Thermo Eclipse", "CID 35%")[:3], ("mixed", "ion_trap", "mixed.tribrid"))
        # the Zorbax Eclipse column is no instrument
        self.assertEqual(self.classify("Agilent 6545 Q-TOF; Zorbax Eclipse Plus C18")[:2], ("qtof", "qtof"))
        self.assertEqual(self.classify("Agilent 6545 Q-TOF; Eclipse XDB-C18")[:2], ("qtof", "qtof"))

    def test_n2_collision_gas_marks_beam_cid_on_trap_hybrids(self):
        got = self.classify("Orbitrap Elite (FTMS), infusion", "CID (N2), ramp 0-50 eV; 2 m/z isolation")
        self.assertEqual(got, ("mixed", "orbitrap_ft", "mixed.ltq_hybrid", "beam (collision gas)"))
        self.assertEqual(self.classify("Orbitrap Elite", "CID with nitrogen, 30 eV")[:2], ("mixed", "orbitrap_ft"))
        # without a beam-cell gas, a bare CID stays trap CID on these hybrids
        self.assertEqual(self.classify("Orbitrap Elite", "CID ramp 0-50 eV")[:2], ("mixed", "ion_trap"))
        self.assertEqual(self.classify("LTQ Orbitrap XL", "CID (He), 35%")[:2], ("mixed", "ion_trap"))

    def test_collision_regime(self):
        cases = [
            # (instrument, collision, regime, basis, in_source_mentioned)
            ("Exactive Orbitrap", "in-source CID 25-100 eV", "in_source", "collision: in-source only", True),
            ("Agilent 6220 TOF", "in-source CID, fragmentor 170-330 V", "in_source", "collision: in-source only",
             True),
            ("Q-TOF Premier", "cone voltage 30 V", "in_source", "collision: source voltage only", False),
            ("Waters Xevo TQ", "cone 30, 110, 140 V", "in_source", "collision: source voltage only", False),
            ("Waters Xevo TQ", "cone 30 V, collision energy 20 eV", "beam", "instrument", False),
            ("Q Exactive", "in-source CID 35 eV, then HCD NCE 27", "beam", "instrument", True),
            ("Agilent 6530 Q-TOF", "CID 20 eV", "beam", "instrument", False),
            ("Agilent 6530 Q-TOF", "fragmentor 175 V; CE 20 eV", "beam", "instrument", False),
            ("Finnigan LCQ Deca ion trap", "NCE 35%", "trap", "instrument", False),
            ("Thermo LTQ FT Ultra", "CID 35%", "trap", "activation", False),
            ("Orbitrap Elite", "CID (N2), ramp 0-50 eV", "beam", "activation", False),
            ("Thermo LTQ Orbitrap XL", "HCD RCE 60% (Orbitrap) and LIT CID MSn", "mixed", "activation", False),
            ("Thermo LTQ-Orbitrap XL", None, "unknown", "activation not named", False),
            ("Bruker solariX 7.0 T FTICR", "IRMPD", "ft_cell", "collision", False),
            ("review across ESI with ion trap and QTOF", None, "unknown", "instrument", False),
            ("review across ESI with ion trap and QTOF", "HCD", "beam", "collision", False),
        ]
        for inst, collision, regime, basis, in_src in cases:
            got = nc.collision_regime(self.classify(inst, collision), inst, collision, None, self.rules)
            self.assertEqual(got, (regime, basis, in_src), f"{inst} | {collision}")
        result = build({"stage3": [stage2("s3-x-001", "p3x", instrument="Exactive Orbitrap",
                                          collision="in-source CID 50 eV", formula="C9H7O3+")]})
        row = result["records"][0]
        self.assertEqual((row["instrument_strat"], row["collision_regime"], row["in_source_mentioned"]),
                         ("orbitrap_ft", "in_source", True))
        self.assertEqual(result["groups"][0]["collision_regimes"], {"in_source": 1})
        self.assertEqual(result["support"][0]["collision_regimes"], {"in_source": 1})

    def test_stratum_override(self):
        maps = nc.load_maps()
        maps["stratum_override"] = [("stage2", re.compile(r"^s2-p3x-"), "orbitrap_ft", "beam", "test"),
                                    ("stage2", re.compile(r"^s2-p3z-"), "qtof", "beam", "unused")]
        recs = {"stage2": [stage2("s2-p3x-001", "p3x", instrument="Orbitrap Elite", collision="CID ramp 0-50 eV",
                                  formula="C9H7O3+"),
                           stage2("s2-p3y-001", "p3y", instrument="Orbitrap Elite", collision="CID ramp 0-50 eV",
                                  formula="C9H7O3+")]}
        result = nc.build({k: (Path(f"{k}.json"), v) for k, v in recs.items()}, maps)
        rows = {r["id"]: r for r in result["records"]}
        self.assertEqual((rows["s2-p3x-001"]["instrument_strat"], rows["s2-p3x-001"]["instrument_strat_computed"],
                          rows["s2-p3x-001"]["collision_regime"], rows["s2-p3x-001"]["collision_regime_basis"],
                          rows["s2-p3x-001"]["instrument_override"]),
                         ("orbitrap_ft", "ion_trap", "beam", "override", "test"))
        self.assertEqual((rows["s2-p3y-001"]["instrument_strat"], rows["s2-p3y-001"]["instrument_override"]),
                         ("ion_trap", ""))
        self.assertEqual(result["overrides"]["stratum_override_rows_used"], ["^s2-p3x-"])
        self.assertEqual(result["overrides"]["stratum_override_rows_unused"], ["^s2-p3z-"])
        self.assertEqual(result["groups"][0]["instrument_strata"], {"ion_trap": 1, "orbitrap_ft": 1})
        self.assertTrue(nc.load_stratum_override())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "stratum.tsv"
            path.write_text("source\tid_pattern\tinstrument_strat\tcollision_regime\treason\n"
                            "stage2\t^s2-x-\torbitrap\tbeam\tx\n", encoding="utf-8")
            with self.assertRaises(nc.NormaliseError):
                nc.load_stratum_override(path)

    def test_hybrid_stratum_reaches_records(self):
        result = build({"stage3": [stage2("s3-x-001", "p3x", instrument="Thermo LTQ FT Ultra",
                                          collision="CID 35%", formula="C9H7O3+")]})
        row = result["records"][0]
        self.assertEqual((row["instrument_class"], row["instrument_strat"], row["instrument_rule"],
                          row["instrument_activation"]), ("mixed", "ion_trap", "mixed.ltq_hybrid", "trap (collision)"))
        self.assertEqual(result["groups"][0]["instrument_strata"], {"ion_trap": 1})

    def test_compilation_unknown_and_unmapped(self):
        self.assertEqual(self.classify("review across ESI with ion trap and QTOF")[0], "compilation")
        self.assertEqual(self.classify("not stated")[0], "unknown")
        self.assertEqual(self.classify(None), ("unknown", "unknown", "missing", ""))
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
                         ("flavonoid_c_glycoside", "neutral_loss", "[M-H]-", "C4H8O4"))
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

    def test_support_rollup_counts_papers_across_adducts(self):
        common = dict(family="acylcarnitine", formula="[C4H5O2]+ 85.0284")
        result = build({"stage2": [stage2("s2-a-001", "p3a", adduct="[M+H]+", standards_based=True, **common),
                                   stage2("s2-b-001", "p3b", adduct="[M]+", **common),
                                   stage2("s2-c-001", "p3c", adduct="not stated", **common),
                                   stage2("s2-d-001", "p3d", adduct="[M+H]+, [M+Na]+", **common),
                                   stage2("s2-e-001", "p3e", kind="veto", adduct="[M+H]+", **common)]})
        groups = {g["adduct"]: g for g in result["groups"]}
        self.assertEqual(sorted(groups), ["[M+H]+", "[M+Na]+", "[M]+", "unstated"])
        self.assertEqual(groups["[M+H]+"]["paper_keys"], ["p3a", "p3d"])  # the veto row is not grouped
        self.assertEqual({g["support_n_papers"] for g in groups.values()}, {4})
        self.assertEqual({g["support_n_standards_papers"] for g in groups.values()}, {1})
        self.assertEqual(len({g["support_id"] for g in groups.values()}), 1)
        (support,) = result["support"]
        self.assertEqual((support["family"], support["kind"], support["ion_mode"], support["formula"]),
                         ("acylcarnitine", "product_ion", "positive", "C4H5O2+"))
        self.assertEqual((support["n_groups"], support["n_papers"], support["n_records"]), (4, 4, 4))
        self.assertEqual(support["adducts"], ["[M+H]+", "[M+Na]+", "[M]+", "unstated"])
        self.assertTrue(all(t["support_id"] == support["support_id"] for t in result["terms"] if t["group_id"]))

    def test_kind_consistent(self):
        result = build({"stage2": [
            stage2("s2-a-001", "p3a", kind="neutral_loss", formula="C6H10O5"),
            stage2("s2-b-001", "p3b", kind="neutral_loss", formula="C9H7O3+"),  # a charged neutral loss
            stage2("s2-c-001", "p3c", kind="fragment_ion", formula="C9H7O3+"),
            stage2("s2-d-001", "p3d", kind="fragment_ion", formula="C9H7O3")]})  # an uncharged product ion
        groups = {(g["kind"], g["formula"]): g for g in result["groups"]}
        self.assertEqual({k: (g["kind_consistent"], g["n_kind_conflict_rows"]) for k, g in groups.items()},
                         {("neutral_loss", "C6H10O5"): (True, 0), ("neutral_loss", "C9H7O3+"): (False, 1),
                          ("product_ion", "C9H7O3+"): (True, 0), ("product_ion", "C9H7O3"): (False, 1)})
        self.assertEqual({(s["kind"], s["formula"]): s["kind_consistent"] for s in result["support"]},
                         {k: g["kind_consistent"] for k, g in groups.items()})

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
                       stage2("s2-a-002", "p3a", kind="veto", formula="C7H5O3+",
                              value="m/z 137 alone is not specific")],
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
            self.assertEqual(names, {"records.tsv", "terms.tsv", "groups.json", "term_support.tsv",
                                     "counts_family_kind.tsv", "counts_family_group.tsv", "manual_queue.tsv",
                                     "mass_flags.tsv", "unmapped.tsv", "run.json"})
            run = json.loads((out / "run.json").read_text(encoding="utf-8"))
            self.assertEqual(run["reconciliation"]["status"], "ok")
            self.assertEqual(run["counts"]["n_records_in"], 5)
            self.assertEqual(run["counts"]["n_records_out"], 5)
            self.assertEqual([i["label"] for i in run["inputs"]], ["npr", "stage2", "stage3"])
            self.assertTrue(all(re.fullmatch(r"[0-9a-f]{64}", i["sha256"]) for i in run["inputs"]))
            self.assertRegex(run["code"]["script_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(set(run["code"]["mapping_files"]), {p.name for p in nc.MAP_FILES})
            self.assertIn("term-override-p3.tsv", run["code"]["mapping_files"])
            self.assertIn("stratum-override-p3.tsv", run["code"]["mapping_files"])
            self.assertIn("collision_regime", run["counts"])
            self.assertEqual(run["counts"]["terms_excluded_from_groups"], {"role_unverified": 1})
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
