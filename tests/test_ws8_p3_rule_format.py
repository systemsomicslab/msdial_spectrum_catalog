"""WS8 P3 rule-file format reference (tools/ws8_p3/rule_format.py), on synthetic rule sets only."""
import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from tools.ws8_p3 import rule_format as rf

ELEMENTS = {"format": "msfinder-diagnostic-elements", "format_version": 1, "electron": 0.000549,
            "elements": {"C": 12.0, "Cl": 34.96885268, "H": 1.007825032, "K": 38.96370668, "N": 14.003074,
                         "Na": 22.98976928, "O": 15.99491462, "P": 30.97376163, "S": 31.972071}}
CITE = [{"doi": "10.1000/synthetic.1", "locator": "Table 1"}]
REF = [{"reference": "synthetic adduct table", "locator": "line 2"}]


def adduct(name, charge, add=None, remove=None, common=True, aliases=(), multiplier=1, origin="as_published"):
    return {"adduct": name, "aliases": list(aliases), "charge": charge, "multiplier": multiplier, "add": add,
            "remove": remove, "common": common, "origin": origin, "source": REF}


def unit(name, cls, residue, origin="as_published"):
    return {"name": name, "class": cls, "residue": residue, "origin": origin, "source": CITE}


def motif(mid, smarts="c1ccccc1", origin="source_defined"):
    return {"motif_id": mid, "smarts": smarts, "origin": origin, "source": CITE,
            "chemont": {"compatible": [{"id": "CHEMONTID:0000334", "name": "Flavonoids"}], "broader": []}}


ADDUCTS = [adduct("[M+H]+", 1, add="H", aliases=["M+H"]), adduct("[M+Na]+", 1, add="Na"),
           adduct("[M]+", 1, common=False, aliases=["[M+]+"]), adduct("[M-H]-", -1, remove="H"),
           adduct("[M+HCOO]-", -1, add="CHO2", aliases=["[M+FA-H]-"]), adduct("[M+Cl]-", -1, add="Cl", common=False)]
UNITS = [unit("Hex", "sugar", "C6H10O5"), unit("dHex", "sugar", "C6H10O4"), unit("Pen", "sugar", "C5H8O4"),
         unit("HexA", "sugar", "C6H8O6"), unit("caffeoyl", "acyl", "C9H6O3")]
MOTIFS = [motif("flavonoid_o_glycoside"), motif("flavonoid"), motif("hydroxycinnamoyl_hexose")]
NESTING = [{"child": "flavonoid_o_glycoside", "parent": "flavonoid", "source": CITE}]
SCAFFOLDS = {"list_id": "aglycones", "entries": [
    {"scaffold_id": "apigenin", "name": "apigenin", "formula": "C15H10O5", "class": "flavone",
     "ions": {"rda": {"pos": ["C7H5O4"], "neg": ["C7H3O4"]}}, "origin": "as_published", "source": CITE},
    {"scaffold_id": "quercetin", "name": "quercetin", "formula": "C15H10O7", "class": "flavonol",
     "ions": {"rda": {"pos": ["C7H5O4"], "neg": ["C7H3O4"]}}, "origin": "as_published", "source": CITE}]}


def term(tid, kind, role, **fields):
    t = {"term_id": tid, "kind": kind, "role": role, "term_source": "src:0", "min_rel_intensity": 0.05,
         "floor_source": "src:0", "than": None, "ratio": None, "ratio_source": None, "note": None}
    t.update(fields)
    return t


def base_rule():
    return {
        "format": "msfinder-diagnostic-rule", "format_version": 1, "rule_id": "flavonoid.o_glycoside.neg",
        "version": 1, "source": copy.deepcopy(CITE), "transcribed_by": "synthetic", "source_checked": False,
        "author_reviewed": False, "measure": True, "review": None, "parameters_origin": "as_published",
        "primary_precondition_mode": "adduct_ambiguous", "calibrated_on": ["QTOF-type"], "convention": "p2",
        "product_tolerance": "p2", "loss_reference": "hypothesis", "veto_excludes_motif": False,
        "claim": {"motif_id": "flavonoid_o_glycoside", "claim_term_ids": ["Y0", "rda"], "claim_level": "class",
                  "non_claims": ["sugar position unresolved"],
                  "positive_examples": [{"pubchem_cid": 1, "smiles": "C"}, {"pubchem_cid": 2, "smiles": "CC"}],
                  "negative_examples": [{"pubchem_cid": 3, "smiles": "CCC"}]},
        "preconditions": {"ion_mode": "neg", "adducts": ["[M-H]-"], "precursor_tolerance_da": 0.01,
                          "composition": {"core": None, "scaffold_list": "aglycones",
                                          "units": {"Hex": [0, 2], "dHex": [0, 1]},
                                          "unit_totals": {"sugar": [1, 2]}},
                          "priors": [], "competitors": []},
        "terms": [term("Y0", "aglycone_ion", "required", form="even", charge_carrier="proton", plus=None, minus=None),
                  term("loss_sugar", "loss_from_precursor", "counted",
                       unit={"class": "sugar", "plus": None, "minus": None}),
                  term("rda", "product_ion", "required", scaffold_ions="rda"),
                  term("caffeate", "product_ion", "veto_present", formula="C9H7O4")],
        "complementary": [],
        "logic": {"all": [{"term": "Y0"}, {"term": "rda"},
                          {"k_of": {"k": 1, "k_source": "src:0", "of": [{"term": "loss_sugar"}]}},
                          {"not": {"term": "caffeate"}}]},
        "notes": "synthetic",
    }


def core_rule(rule_id, logic, terms, claim, core="C15H10O5", mode="neg", adducts=("[M-H]-",), measure=True):
    r = base_rule()
    r["rule_id"] = rule_id
    r["measure"] = measure
    r["preconditions"]["ion_mode"] = mode
    r["preconditions"]["adducts"] = list(adducts)
    r["preconditions"]["composition"] = {"core": core, "scaffold_list": None, "units": {}, "unit_totals": {}}
    r["terms"] = terms
    r["logic"] = logic
    r["claim"]["claim_term_ids"] = claim
    return r


def dumps(o):
    return rf.storage_bytes(o)


def wrap(fmt, key, rows):
    return {"format": fmt, "format_version": 1, key: rows}


class Set:
    """A synthetic rule set as plain objects; files() renders it."""

    def __init__(self, rules=None):
        self.ruleset = {"format": "msfinder-diagnostic-ruleset", "format_version": 1, "rule_set_version": "0.0.1-test"}
        self.elements = copy.deepcopy(ELEMENTS)
        self.adducts = copy.deepcopy(ADDUCTS)
        self.units = copy.deepcopy(UNITS)
        self.motifs = copy.deepcopy(MOTIFS)
        self.nesting = copy.deepcopy(NESTING)
        self.scaffolds = {"aglycones": copy.deepcopy(SCAFFOLDS)}
        self.rules = [base_rule()] if rules is None else rules
        self.extra = {}
        self.drop = set()

    @property
    def rule(self):
        return self.rules[0]

    def files(self):
        f = {"ruleset.json": dumps(self.ruleset), "tables/elements.json": dumps(self.elements),
             "tables/adducts.json": dumps(wrap("msfinder-diagnostic-adducts", "adducts", self.adducts)),
             "tables/units.json": dumps(wrap("msfinder-diagnostic-units", "units", self.units)),
             "tables/motifs.json": dumps(wrap("msfinder-diagnostic-motifs", "motifs", self.motifs)),
             "tables/motif_nesting.json": dumps(wrap("msfinder-diagnostic-motif-nesting", "nesting", self.nesting))}
        for lid, doc in self.scaffolds.items():
            f[f"scaffolds/{lid}.json"] = dumps(doc)
        for r in self.rules:
            f[f"rules/{r['rule_id']}.json"] = dumps(r)
        f.update(self.extra)
        for p in self.drop:
            f.pop(p, None)
        return f

    def load(self, draft=True, **kw):
        return rf.load_files(self.files(), draft=draft, **kw)

    def codes(self, draft=True):
        return {e.code for e in self.load(draft).errors}

    def errors(self, draft=True):
        return [(e.code, e.file, e.pointer) for e in self.load(draft).errors]


RULE_FILE = "rules/flavonoid.o_glycoside.neg.json"


def code_of(b):
    try:
        rf.load_strict(b)
    except rf.RuleFormatError as e:
        return e.code
    return None


class StrictReaderTests(unittest.TestCase):
    def test_reader_codes(self):
        cases = {
            b"\xef\xbb\xbf{}": "bom", b'{"a":"\xff"}': "invalid_utf8", b'{"a":1,}': "syntax",
            b'{"a":1 /* c */}': "syntax", b"{'a':1}": "syntax", b'{"a":"x\ny"}': "syntax", b"[]": "top_level_not_object",
            b'"x"': "top_level_not_object", b'{"a":1,"a":2}': "duplicate_key", b'{"a":NaN}': "nonfinite",
            b'{"a":-Infinity}': "nonfinite", b'{"a":12345678901234567}': "number_range",
            b'{"a":' + b"9" * 4301 + b"}": "number_range", b'{"a":9007199254740992}': "number_range",
            b'{"a":-0}': "negative_zero", b'{"a":-0.0}': "negative_zero", b'{"a":-0.000}': "negative_zero",
            b'{"a":1e5}': "number_form", b'{"a":1.5E3}': "number_form", b'{"a":0.1234567890123456}': "float_precision",
            b'{"a":0.00001}': "number_range", b'{"a":1000000000000000.0}': "number_range",
            b'{"a":"\\u007f"}': "control_char", b'{"a":"\x7f"}': "control_char", b'{"a":"\\u0080"}': "non_ascii",
            '{"a":"é"}'.encode(): "non_ascii", b'{"a":"\\ud83d\\ude00"}': "non_ascii",
            b'{"a":"\\ud800"}': "non_ascii", b'{"\\u00e9":1}': "non_ascii",
        }
        for b, code in cases.items():
            with self.subTest(b=b[:40]):
                self.assertEqual(code_of(b), code)

    def test_nesting_depth_boundary(self):
        self.assertIsNone(code_of(b'{"a":' + b"[" * 127 + b"]" * 127 + b"}"))
        self.assertEqual(code_of(b'{"a":' + b"[" * 128 + b"]" * 128 + b"}"), "nesting_depth")
        self.assertEqual(code_of(b'{"a":' + b"[" * 5000 + b"]" * 5000 + b"}"), "nesting_depth")

    def test_accepted_numbers(self):
        o = rf.load_strict(b'{"a":9007199254740991,"b":0.0001,"c":123456789012345.0,"d":0.0,"e":-5}')
        self.assertEqual(o, {"a": 2 ** 53 - 1, "b": 0.0001, "c": 123456789012345.0, "d": 0.0, "e": -5})

    def test_canonical_vectors(self):
        o = rf.load_strict(b'{"b":1,"B":[0.050,5.00,10.0,0.0001,123456789012345.0],"a":"q\\"\\\\"}')
        self.assertEqual(rf.canonical_bytes(o),
                         b'{"B":[0.05,5.0,10.0,0.0001,123456789012345.0],"a":"q\\"\\\\","b":1}')

    def test_crlf_and_whitespace_do_not_change_hashes(self):
        a = b'{\r\n\t"x": [1, 2.50],\r\n "y": null\r\n}\r\n'
        b = b'{"y":null,"x":[1,2.5]}'
        self.assertEqual(rf.file_sha256(a), rf.file_sha256(b))
        self.assertEqual(rf.file_sha256(b), hashlib.sha256(b'{"x":[1,2.5],"y":null}').hexdigest())

    def test_ruleset_sha256_definition(self):
        files = {"ruleset.json": b'{"a":1}', "rules/x.json": b'{"b": true}'}
        m = {"files": {p: hashlib.sha256(rf.canonical_bytes(json.loads(b))).hexdigest() for p, b in files.items()},
             "format": "msfinder-diagnostic-ruleset-hash", "format_version": 1}
        expected = hashlib.sha256(json.dumps(m, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        self.assertEqual(rf.ruleset_sha256(files), expected)

    def test_rule_semantic_hash_ignores_review_fields(self):
        r = base_rule()
        h = rf.rule_semantic_sha256(r, {"x": 1})
        r2 = dict(r, author_reviewed=True, review={"by": "x"})
        self.assertEqual(rf.rule_semantic_sha256(r2, {"x": 1}), h)
        self.assertNotEqual(rf.rule_semantic_sha256(r, {"x": 2}), h)
        self.assertNotEqual(rf.rule_semantic_sha256(r, {"x": 1}, "msfinder-diagnostic-rules/2"), h)

    def test_storage_form(self):
        self.assertEqual(rf.storage_bytes({"b": [1, 2.50], "a": "x"}),
                         b'{\n "a": "x",\n "b": [\n  1,\n  2.5\n ]\n}\n')

    def test_write_storage_refuses_what_the_reader_refuses(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.json"
            for obj, code in (({"a": "é"}, "non_ascii"), ({"a": 0.00001}, "number_form"),
                              ({"a": "x\ty"}, "control_char")):
                with self.assertRaises(rf.RuleFormatError) as cm:
                    rf.write_storage(p, obj)
                self.assertEqual(cm.exception.code, code)
            self.assertTrue(rf.write_storage(p, {"b": 1, "a": 0.050}))
            self.assertEqual(p.read_bytes(), b'{\n "a": 0.05,\n "b": 1\n}\n')
            self.assertFalse(rf.write_storage(p, {"a": 0.05, "b": 1}))


class FormulaAndMassTests(unittest.TestCase):
    el = rf.Elements(0.000549, 549_000_000, ELEMENTS["elements"],
                     {s: rf.lexeme_pda(m) for s, m in ELEMENTS["elements"].items()})

    def mass(self, f):
        return rf.mz_of(self.el.mass_pda(rf.parse_formula(f, self.el.pda)))

    def ion(self, f, z):
        return rf.mz_of(self.el.ion_pda(rf.parse_formula(f, self.el.pda), z))

    def test_formula_grammar(self):
        ok = {"C6H10O5": {"C": 6, "H": 10, "O": 5}, "H2O": {"H": 2, "O": 1}, "ClH": {"Cl": 1, "H": 1},
              "H4N": {"H": 4, "N": 1}, "C6H5Cl": {"C": 6, "H": 5, "Cl": 1}, "C7H12NO9S2": None, "O3S": None,
              "CHO2": None, "C999": None, "Na": None}
        for f, comp in ok.items():
            got = rf.parse_formula(f, self.el.pda)
            if comp:
                self.assertEqual(got, comp)
        bad = {"C1H4": "bad_formula", "C0": "bad_formula", "c6": "bad_formula", "C1000": "bad_formula",
               "C06": "bad_formula", "": "bad_formula", "C6 H6": "bad_formula", "C6H6+": "bad_formula",
               "HCl": "formula_not_hill", "NH4": "formula_not_hill", "H2OC": "formula_not_hill",
               "CH4C": "formula_not_hill", "OH": "formula_not_hill", "C6H10O5Br": "unknown_element", "Xx": "unknown_element"}
        for f, code in bad.items():
            with self.subTest(f=f):
                with self.assertRaises(rf.RuleFormatError) as cm:
                    rf.parse_formula(f, self.el.pda)
                self.assertEqual(cm.exception.code, code)
        self.assertEqual(rf.hill_format({"O": 2, "C": 1, "H": 0}), "CO2")
        self.assertEqual(rf.hill_format({"H": 1, "Cl": 1}), "ClH")

    def test_lexeme_pda(self):
        self.assertEqual(rf.lexeme_pda(1.007825032), 1_007_825_032_000)
        self.assertEqual(rf.lexeme_pda(0.000549), 549_000_000)
        self.assertEqual(rf.lexeme_pda(12.0), 12 * 10 ** 12)
        self.assertEqual(rf.lexeme_pda(34.96885268), 34_968_852_680_000)
        with self.assertRaises(rf.RuleFormatError) as cm:
            rf.lexeme_pda(1.0078250321234)
        self.assertEqual(cm.exception.code, "mass_precision")

    def test_masses_against_hand_values(self):
        hand = {  # each summed by hand from the element literals and the electron 0.000549
            ("C6H10O5", 0): 162.05282342, ("C6H10O4", 0): 146.0579088, ("C5H8O4", 0): 132.04225848,
            ("C6H8O6", 0): 176.03208806, ("C9H6O3", 0): 162.031694, ("C9H7O4", -1): 179.034982704,
            ("C15H17O9", -1): 341.087806124, ("C21H21O11", 1): 449.107837492, ("C15H11O5", 1): 271.060099452,
            ("C7H5O4", 1): 153.01823464, ("C7H3O4", -1): 151.003682576, ("O3S", 0): 79.95681486,
        }
        for (f, z), value in hand.items():
            with self.subTest(f=f):
                got = self.mass(f) if z == 0 else self.ion(f, z)
                self.assertAlmostEqual(got, value, delta=1e-5)

    def test_hypothesis_precursors_and_targets(self):
        res = Set().load()
        self.assertEqual(res.errors, [])
        rule = res.rules["flavonoid.o_glycoside.neg"]
        hyps = {h.hypothesis_id: h for h, _ in rf.iter_hypotheses(rule, res.tables)}
        h = hyps["apigenin|Hex1|[M-H]-"]
        self.assertAlmostEqual(rf.mz_of(h.neutral_pda), 432.10564684, delta=1e-5)
        self.assertAlmostEqual(rf.mz_of(h.precursor_pda), 431.098370808, delta=1e-5)
        self.assertAlmostEqual(rf.mz_of(hyps["apigenin|Hex1|[M+HCOO]-"].precursor_pda), 477.103850112, delta=1e-5)
        targets = {t.tid: [rf.mz_of(m[0]) for m in rf.term_members(t, h, rule, res.tables)] for t in rule.terms}
        self.assertAlmostEqual(targets["Y0"][0], 269.045547388, delta=1e-5)
        self.assertAlmostEqual(targets["loss_sugar"][0], 431.098370808 - 162.05282342, delta=1e-5)
        self.assertAlmostEqual(targets["rda"][0], 151.003682576, delta=1e-5)
        self.assertAlmostEqual(targets["caffeate"][0], 179.034982704, delta=1e-5)

    def test_aglycone_forms_and_adduct_carrier(self):
        s = Set()
        r = s.rule
        r["terms"] += [
            term("Yr", "aglycone_ion", "upgrade", form="radical", charge_carrier="proton", plus=None, minus=None),
            term("Ya", "aglycone_ion", "upgrade", form="even", charge_carrier="adduct", plus="C", minus=None)]
        res = s.load()
        self.assertEqual(res.errors, [])
        rule = res.rules[r["rule_id"]]
        hyps = {h.hypothesis_id: h for h, _ in rf.iter_hypotheses(rule, res.tables)}
        h = hyps["apigenin|Hex1|[M+HCOO]-"]
        got = {t.tid: rf.mz_of(rf.term_members(t, h, rule, res.tables)[0][0]) for t in rule.terms}
        self.assertAlmostEqual(got["Yr"], 268.037722356, delta=1e-5)
        self.assertAlmostEqual(got["Ya"], 270.05282342 + 44.997654272 + 12.0 + 0.000549, delta=1e-5)

    def test_flavylium_intrinsic_charge(self):
        r = core_rule("anthocyanin.glucoside.pos", {"all": [{"term": "a"}, {"term": "b"}]},
                      [term("a", "product_ion", "required", formula="C15H11O6"),
                       term("b", "loss_from_precursor", "required", formula="C6H10O5")], ["a", "b"],
                      core="C21H21O11", mode="pos", adducts=("[M]+",))
        res = Set([r]).load()
        self.assertEqual(res.errors, [])
        rule = res.rules[r["rule_id"]]
        hyps = {h.hypothesis_id: h for h, _ in rf.iter_hypotheses(rule, res.tables)}
        self.assertEqual(list(hyps), ["core||[M]+", "core||[M+H]+", "core||[M+Na]+"])
        self.assertAlmostEqual(rf.mz_of(hyps["core||[M]+"].precursor_pda), 449.107837492, delta=1e-5)

    def test_enumeration_order_and_ion_series(self):
        s = Set()
        s.rule["terms"].append(term("ys", "ion_series", "upgrade", unit_classes=["sugar"], charge_carrier="proton",
                                    form="even", plus=None, minus=None))
        res = s.load()
        self.assertEqual(res.errors, [])
        rule = res.rules[s.rule["rule_id"]]
        ids = [h.hypothesis_id for h, _ in rf.iter_hypotheses(rule, res.tables)]
        self.assertEqual(ids[:8], ["apigenin|dHex1|[M-H]-", "apigenin|dHex1|[M+HCOO]-", "apigenin|Hex1|[M-H]-",
                                   "apigenin|Hex1|[M+HCOO]-", "apigenin|Hex1,dHex1|[M-H]-",
                                   "apigenin|Hex1,dHex1|[M+HCOO]-", "apigenin|Hex2|[M-H]-", "apigenin|Hex2|[M+HCOO]-"])
        self.assertEqual(len(ids), 16)
        self.assertEqual(rf.count_series(rule, res.tables), 8)
        hyps = {h.hypothesis_id: h for h, _ in rf.iter_hypotheses(rule, res.tables)}
        ys = rule.term_by_id["ys"]
        keys = [m[1] for m in rf.term_members(ys, hyps["quercetin|Hex1,dHex1|[M-H]-"], rule, res.tables)]
        self.assertEqual(keys, ["dHex1", "Hex1"])
        self.assertEqual([m[1] for m in rf.term_members(ys, hyps["quercetin|Hex2|[M-H]-"], rule, res.tables)], ["Hex1"])
        loss = rule.term_by_id["loss_sugar"]
        self.assertEqual([m[1] for m in rf.term_members(loss, hyps["quercetin|Hex1,dHex1|[M-H]-"], rule, res.tables)],
                         ["Hex", "dHex"])

    def test_caffeoyl_hexose_loss_and_caffeate_coincide_in_pda(self):
        r = core_rule("hydroxycinnamoyl.hexose_ester.neg", {"all": [{"term": "loss_hex"}, {"term": "caffeate"}]},
                      [term("loss_hex", "loss_from_precursor", "required", formula="C6H10O5"),
                       term("caffeate", "product_ion", "required", formula="C9H7O4")],
                      ["loss_hex", "caffeate"], core="C15H18O9", measure=False)
        r["claim"]["motif_id"] = "hydroxycinnamoyl_hexose"
        r["complementary"] = [["loss_hex", "caffeate"]]
        res = Set([r]).load()
        self.assertEqual(res.errors, [])
        rule = res.rules[r["rule_id"]]
        h = next(h for h, _ in rf.iter_hypotheses(rule, res.tables))
        self.assertAlmostEqual(rf.mz_of(h.precursor_pda), 341.087806124, delta=1e-5)
        a, b = (rf.term_members(t, h, rule, res.tables)[0][0] for t in rule.terms)
        self.assertEqual(a, b)


def defect(mutate, draft=True):
    s = Set()
    mutate(s)
    return s


def set_rule(path, value):
    def m(s):
        o = s.rule
        for k in path[:-1]:
            o = o[k]
        o[path[-1]] = value
    return m


def t_of(s, tid):
    return next(t for t in s.rule["terms"] if t["term_id"] == tid)


def raw_rule(b):
    def m(s):
        s.extra[RULE_FILE] = b
    return m


def deep(n):
    node = {"term": "Y0"}
    for _ in range(n):
        node = {"all": [node]}
    return node


DEFECTS = {
    "bom": raw_rule(b"\xef\xbb\xbf{}"),
    "invalid_utf8": raw_rule(b'{"a":"\xff"}'),
    "syntax": raw_rule(b'{"a":1,}'),
    "top_level_not_object": raw_rule(b"[]"),
    "duplicate_key": raw_rule(b'{"a":1,"a":2}'),
    "nonfinite": raw_rule(b'{"a":NaN}'),
    "nesting_depth": raw_rule(b'{"a":' + b"[" * 129 + b"]" * 129 + b"}"),
    "number_range": raw_rule(b'{"a":12345678901234567}'),
    "negative_zero": raw_rule(b'{"a":-0}'),
    "number_form": raw_rule(b'{"a":1e5}'),
    "float_precision": raw_rule(b'{"a":0.1234567890123456}'),
    "control_char": raw_rule(b'{"a":"\\u007f"}'),
    "non_ascii": raw_rule(b'{"a":"\\u00e9"}'),
    "missing_file": lambda s: s.drop.add("tables/units.json"),
    "unexpected_file": lambda s: s.extra.update({"tables/extra.json": b"{}"}),
    "unknown_key": set_rule(["extra"], 1),
    "missing_key": lambda s: s.rule.pop("notes"),
    "bad_type": set_rule(["version"], "1"),
    "bad_value": set_rule(["convention"], "other"),
    "kind_field_choice": lambda s: t_of(s, "rda").update(formula="C7H3O4"),
    "formula_dependent_key": set_rule(["preconditions", "oracle_formula"], "C6H12O6"),
    "typed_mass": lambda s: t_of(s, "caffeate").update(mz=179.035),
    "out_of_scope_kind": lambda s: t_of(s, "caffeate").update(kind="chain_set"),
    "out_of_scope_variable": set_rule(["logic", "all", 2, "k_of", "k"], {"var": "n_chains", "offset": 0, "min": 1}),
    "bad_formula": lambda s: t_of(s, "caffeate").update(formula="C9h7O4"),
    "formula_not_hill": lambda s: t_of(s, "caffeate").update(formula="H7C9O4"),
    "unknown_element": lambda s: t_of(s, "caffeate").update(formula="C9H7BrO4"),
    "mass_precision": lambda s: s.elements["elements"].update(H=1.0078250321234),
    "missing_source": lambda s: s.units[4].update(source=[]),
    "duplicate_id": lambda s: s.units.append(unit("Hex", "sugar", "C6H10O5")),
    "adduct_polarity": lambda s: s.adducts[0].update(charge=-1),
    "adduct_name_mismatch": lambda s: s.adducts[1].update(add="K"),
    "missing_required_unit": lambda s: s.units.pop(2),
    "nesting_cycle": lambda s: s.nesting.append({"child": "flavonoid", "parent": "flavonoid_o_glycoside",
                                                 "source": CITE}),
    "file_name_mismatch": lambda s: s.extra.update({RULE_FILE: dumps(dict(base_rule(), rule_id="flavonoid.other.neg"))}),
    "mode_mismatch": set_rule(["preconditions", "ion_mode"], "pos"),
    "unknown_term": set_rule(["logic", "all", 0, "term"], "nope"),
    "unknown_unit": set_rule(["preconditions", "composition", "units", "Xyl"], [0, 1]),
    "unknown_scaffold_list": set_rule(["preconditions", "composition", "scaffold_list"], "nope"),
    "unknown_ion_key": lambda s: t_of(s, "rda").update(scaffold_ions="nokey"),
    "unknown_adduct": set_rule(["preconditions", "adducts"], ["[M+Br]-"]),
    "unknown_motif": set_rule(["claim", "motif_id"], "nope"),
    "unknown_source_ref": lambda s: t_of(s, "Y0").update(term_source="src:5"),
    "unknown_variable": set_rule(["logic", "all", 2, "k_of", "k"], {"var": "n_foo", "offset": 0, "min": 1}),
    "duplicate_term": lambda s: s.rule["terms"].append(term("Y0", "product_ion", "upgrade", formula="C7H3O4")),
    "unused_term": lambda s: s.rule["terms"].append(term("x", "product_ion", "required", formula="C7H3O4")),
    "claim_term_role": set_rule(["claim", "claim_term_ids"], ["Y0", "caffeate"]),
    "claim_terms_too_few": set_rule(["claim", "claim_term_ids"], ["Y0"]),
    "role_placement": set_rule(["logic", "all", 3], {"not": {"term": "Y0"}}),
    "logic_node": set_rule(["logic", "all", 0], {"xor": []}),
    "logic_depth": set_rule(["logic", "all", 0], deep(16)),
    "k_range": set_rule(["logic", "all", 2, "k_of", "k"], 0),
    "series_too_large": set_rule(["preconditions", "composition"], {
        "core": None, "scaffold_list": "aglycones",
        "units": {"Hex": [0, 100], "Pen": [0, 100], "HexA": [0, 100], "dHex": [0, 1]},
        "unit_totals": {"sugar": [1, 400]}}),
    "series_empty": set_rule(["preconditions", "composition", "unit_totals"], {"sugar": [1, 2], "acyl": [1, 1]}),
    "negative_count": set_rule(["preconditions", "composition", "units", "Hex"], [-1, 2]),
    "mass_range": lambda s: s.rule["terms"].append(term("u", "product_ion", "upgrade",
                                                        unit={"name": "Hex", "plus": None, "minus": "C6H10O5"})),
    "complementary_pair": set_rule(["complementary"], [["Y0", "Y0"]]),
    "complementary_multi_member": set_rule(["complementary"], [["loss_sugar", "rda"]]),
    "single_ion_claim": set_rule(["logic"], {"any": [{"term": "Y0"}, {"term": "rda"},
                                                     {"k_of": {"k": 1, "k_source": "src:0",
                                                               "of": [{"term": "loss_sugar"}]}},
                                                     {"not": {"term": "caffeate"}}]}),
    "claim_single_ion": set_rule(["logic"], {"all": [{"term": "Y0"},
                                                     {"k_of": {"k": 1, "k_source": "src:0",
                                                               "of": [{"term": "loss_sugar"}]}},
                                                     {"any": [{"term": "rda"}, {"not": {"term": "caffeate"}}]}]}),
    "support_target_is_precursor": set_rule(["preconditions", "composition", "unit_totals"], {"sugar": [0, 2]}),
    "origin_inconsistent": set_rule(["parameters_origin"], "tuned_vs20"),
    "review_missing": lambda s: s.rule.update(author_reviewed=True, source_checked=True),
    "review_not_allowed": lambda s: reviewed_with(s, source_checked=False),
    "review_stale": lambda s: s.rule.update(author_reviewed=True, source_checked=True,
                                            review=review_block("0" * 64)),
}


def reviewed_with(s, **fields):
    """Give s.rule a review whose rule_semantic_sha256 is current, after applying fields."""
    s.rule.update({"author_reviewed": True, "source_checked": True, **fields})
    s.rule["review"] = review_block("0" * 64)
    s.rule["review"]["rule_semantic_sha256"] = s.load().semantic_sha256(s.rule["rule_id"])
NON_DRAFT_DEFECTS = {
    "history_format": lambda s: s.extra.update({rf.HISTORY_FILE: b"0.0.1-test\tnot-a-sha\t2026-11-18\n"}),
    "ruleset_unpinned": lambda s: s.extra.update({rf.HISTORY_FILE: b""}),
    "ruleset_sha_mismatch": lambda s: s.extra.update({rf.HISTORY_FILE: b"0.0.1-test\t" + b"0" * 64 + b"\t2026-11-18\n"}),
}


def review_block(sha):
    return {"by": "author", "at": "2026-12-01", "basis": "benchmark", "run_id": "p3cs-20261201-abc",
            "rules_precision_sha256": "a" * 64, "review_proposal_sha256": "b" * 64, "rule_semantic_sha256": sha,
            "instrument_classes": ["QTOF-type"], "product_tolerance": {"abs_da": 0.005, "ppm": 10.0, "source": "p2"},
            "precursor_tolerance_da": 0.01}


class ErrorCodeTests(unittest.TestCase):
    def test_base_set_is_clean(self):
        res = Set().load()
        self.assertEqual(res.errors, [])
        self.assertEqual(res.n_rules, 1)

    def test_every_code_has_a_single_defect_case(self):
        self.assertEqual(set(DEFECTS) | set(NON_DRAFT_DEFECTS), set(rf.ERROR_CODES))

    def test_each_defect_reports_its_code(self):
        for code, mutate in DEFECTS.items():
            with self.subTest(code=code):
                s = Set()
                mutate(s)
                self.assertIn(code, s.codes())
        for code in ("review_missing", "review_not_allowed", "review_stale", "missing_source", "claim_single_ion",
                     "support_target_is_precursor", "complementary_multi_member", "origin_inconsistent"):
            with self.subTest(single=code):
                s = Set()
                DEFECTS[code](s)
                self.assertEqual(s.codes(), {code})
        for code, mutate in NON_DRAFT_DEFECTS.items():
            with self.subTest(code=code):
                s = Set()
                mutate(s)
                self.assertIn(code, s.codes(draft=False))

    def test_errors_are_sorted_with_pointers(self):
        s = Set()
        s.rule["version"] = "1"
        s.rule["convention"] = "other"
        t_of(s, "caffeate")["formula"] = "H7C9O4"
        errs = s.errors()
        self.assertEqual(errs, sorted(errs, key=lambda e: (e[1], e[2], e[0])))
        self.assertIn(("formula_not_hill", RULE_FILE, "/terms/3/formula"), errs)
        self.assertIn(("bad_value", RULE_FILE, "/convention"), errs)
        self.assertIn(("bad_type", RULE_FILE, "/version"), errs)

    def test_refused_keys_anywhere(self):
        for key, code in (("inchikey", "formula_dependent_key"), ("predicted_formula", "formula_dependent_key"),
                          ("exact_mass", "typed_mass"), ("max_rel_intensity", "out_of_scope_kind")):
            with self.subTest(key=key):
                s = Set()
                s.rule["claim"]["positive_examples"][0][key] = "x"
                self.assertEqual(s.errors(), [(code, RULE_FILE, "/claim/positive_examples/0/" + key)])

    def test_pointer_escaping(self):
        self.assertEqual(rf.jp("", "a/b", "c~d", 0), "/a~1b/c~0d/0")

    def test_roles(self):
        cases = {
            "counted outside k_of": set_rule(["logic", "all", 2], {"term": "loss_sugar"}),
            "required under k_of": set_rule(["logic", "all", 2, "k_of", "of"], [{"term": "loss_sugar"}, {"term": "Y0"}]),
            "veto outside not": set_rule(["logic", "all", 3], {"term": "caffeate"}),
            "nested not": set_rule(["logic", "all", 3], {"not": {"not": {"term": "caffeate"}}}),
            "k_of inside not": set_rule(["logic", "all", 3], {"not": {"k_of": {"k": 1, "k_source": "src:0",
                                                                              "of": [{"term": "caffeate"}]}}}),
            "upgrade in logic": lambda s: (s.rule["terms"].append(term("up", "product_ion", "upgrade",
                                                                       formula="C7H3O4")),
                                           s.rule["logic"]["all"].append({"term": "up"})),
            "window as support": lambda s: t_of(s, "caffeate").update(
                kind="window", role="required", **{"as": "ion", "half_width_da": 0.5}),
        }
        for name, mutate in cases.items():
            with self.subTest(name=name):
                s = Set()
                mutate(s)
                self.assertIn("role_placement", s.codes())

    def test_window_veto_and_ordering_are_accepted(self):
        s = Set()
        t_of(s, "caffeate").update(kind="window", role="veto_stronger_than", than=["Y0"], ratio=0.1,
                                   ratio_source="src:0", **{"as": "ion", "half_width_da": 0.5})
        t_of(s, "rda").update(role="ordering", than=["Y0"])
        self.assertEqual(s.errors(), [])

    def test_than_and_ratio_consistency(self):
        s = Set()
        t_of(s, "rda").update(ratio=0.5, ratio_source="src:0")
        self.assertIn(("bad_value", RULE_FILE, "/terms/2/ratio"), s.errors())
        s = Set()
        t_of(s, "caffeate").update(role="veto_stronger_than")
        self.assertIn(("bad_value", RULE_FILE, "/terms/3/than"), s.errors())
        s = Set()
        t_of(s, "caffeate").update(role="veto_stronger_than", than=["Y0"], ratio=0.1)
        self.assertIn(("bad_value", RULE_FILE, "/terms/3/ratio_source"), s.errors())

    def test_adduct_polarity_in_rule(self):
        s = Set()
        s.rule["preconditions"]["adducts"] = ["[M+H]+"]
        self.assertEqual(s.errors(), [("adduct_polarity", RULE_FILE, "/preconditions/adducts/0")])

    def test_duplicate_alias(self):
        s = Set()
        s.adducts[1]["aliases"] = ["M + H"]
        self.assertEqual(s.errors(), [("duplicate_id", "tables/adducts.json", "/adducts/1/aliases/0")])

    def test_min_rel_intensity_scale_follows_convention(self):
        s = Set()
        t_of(s, "Y0")["min_rel_intensity"] = 5.0
        self.assertIn(("bad_value", RULE_FILE, "/terms/0/min_rel_intensity"), s.errors())
        s.rule["convention"] = "msdial"
        self.assertEqual(s.errors(), [])
        t_of(s, "Y0")["min_rel_intensity"] = 5
        self.assertEqual(s.errors(), [("bad_type", RULE_FILE, "/terms/0/min_rel_intensity")])


class TwoIonTests(unittest.TestCase):
    X = term("X", "product_ion", "required", formula="C7H3O4")
    Y = term("Y", "product_ion", "required", formula="C8H5O3")
    Z = term("Z", "product_ion", "required", formula="C6H5O")
    W = term("W", "loss_from_precursor", "required", formula="H2O")

    def run_rule(self, logic, terms, claim, measure=True, **kw):
        r = core_rule("flavonoid.synthetic.neg", logic, copy.deepcopy(terms), claim, measure=measure, **kw)
        return {e.code for e in Set([r]).load().errors}

    def test_design_example_one(self):
        # all(X, loss H2O) with claim [X]: passes pass 1, refused with claim_terms_too_few and claim_single_ion
        got = self.run_rule({"all": [{"term": "X"}, {"term": "W"}]}, [self.X, self.W], ["X"])
        self.assertEqual(got, {"claim_terms_too_few", "claim_single_ion"})

    def test_design_example_two(self):
        # all(X, loss H2O, any(Y, Z)) with [X, Y]: passes the static check and pass 1, fails pass 2
        got = self.run_rule({"all": [{"term": "X"}, {"term": "W"}, {"any": [{"term": "Y"}, {"term": "Z"}]}]},
                            [self.X, self.W, self.Y, self.Z], ["X", "Y"])
        self.assertEqual(got, {"claim_single_ion"})

    def test_complementary_group_is_one_unit(self):
        terms = [term("loss_hex", "loss_from_precursor", "required", formula="C6H10O5"),
                 term("caffeate", "product_ion", "required", formula="C9H7O4")]
        logic = {"all": [{"term": "loss_hex"}, {"term": "caffeate"}]}

        def run(measure, pairs):
            r = core_rule("hydroxycinnamoyl.hexose_ester.neg", logic, copy.deepcopy(terms), ["loss_hex", "caffeate"],
                          core="C15H18O9", measure=measure)
            r["claim"]["motif_id"] = "hydroxycinnamoyl_hexose"
            r["complementary"] = pairs
            return {e.code for e in Set([r]).load().errors}

        self.assertIn("single_ion_claim", run(True, [["loss_hex", "caffeate"]]))
        self.assertEqual(run(False, [["loss_hex", "caffeate"]]), set())
        self.assertEqual(run(True, []), set())  # without the pair, two units (overlap is the evaluator's check)

    def test_two_distinct_ions_pass(self):
        got = self.run_rule({"all": [{"term": "X"}, {"term": "Y"}]}, [self.X, self.Y], ["X", "Y"])
        self.assertEqual(got, set())

    def test_ion_series_counted_with_k_two(self):
        s = Set()
        r = s.rule
        r["terms"] = [term("ys", "ion_series", "counted", unit_classes=["sugar"], charge_carrier="proton", form="even",
                           plus=None, minus=None),
                      term("Y0", "aglycone_ion", "required", form="even", charge_carrier="proton", plus=None,
                           minus=None)]
        r["preconditions"]["composition"]["units"] = {"Hex": [0, 3]}
        r["preconditions"]["composition"]["unit_totals"] = {"sugar": [3, 3]}
        r["claim"]["claim_term_ids"] = ["ys"]
        r["logic"] = {"all": [{"term": "Y0"}, {"k_of": {"k": 2, "k_source": "src:0", "of": [{"term": "ys"}]}}]}
        self.assertEqual(s.errors(), [])
        r["logic"]["all"][1]["k_of"]["k"] = {"var": "n_hex", "offset": -1, "min": 2}
        self.assertEqual(s.errors(), [])
        r["logic"]["all"][1]["k_of"]["k"] = 1
        self.assertEqual({e[0] for e in s.errors()}, {"claim_terms_too_few", "claim_single_ion"})

    def test_member_in_precursor_region_is_never_present(self):
        # C16H13O4- is 0.036 Da above the apigenin [M-H]- precursor: never an ion there, an ion under [M+HCOO]-.
        # Monotone logic cannot gain a violation from a member that is never present, so the rule passes.
        r = core_rule("flavonoid.synthetic.neg", {"all": [{"term": "X"}, {"term": "N"}]},
                      [copy.deepcopy(self.X), term("N", "product_ion", "required", formula="C16H13O4")], ["X", "N"])
        s = Set([r])
        res = s.load()
        self.assertEqual(res.errors, [])
        rule = res.rules[r["rule_id"]]
        hyps = list(rf.iter_hypotheses(rule, res.tables))
        support = [t for t in rule.terms]
        zero = {h.adduct.name: tuple(tuple(abs(m[0] - h.precursor_pda) <= rf.PRECURSOR_REGION_PDA
                                           for m in rf.term_members(t, h, rule, res.tables)) for t in support)
                for h, _ in hyps}
        self.assertEqual(zero, {"[M-H]-": ((False,), (True,)), "[M+HCOO]-": ((False,), (False,))})
        # with N never present, X alone would be one ion: both passes hold under [M-H]-
        self.assertEqual(rf.two_ion_passes(rule, support, zero["[M-H]-"], (), {}, {0, 1}), [])
        lines = rf.describe_masses(res, r["rule_id"])
        self.assertTrue(any(x.startswith("  T\tN\t") and x.endswith("precursor_region") for x in lines))

    def test_unmeasured_rule_skips_the_claim_checks(self):
        got = self.run_rule({"all": [{"term": "X"}, {"term": "W"}]}, [self.X, self.W], ["X"], measure=False)
        self.assertEqual(got, set())


class PrecursorTargetTests(unittest.TestCase):
    def test_y0_at_zero_units_is_the_m_minus_h_precursor(self):
        s = Set()
        s.rule["preconditions"]["composition"]["unit_totals"] = {"sugar": [0, 2]}
        self.assertEqual(s.errors(), [("support_target_is_precursor", RULE_FILE, "/terms/0")])
        s.rule["preconditions"]["composition"]["unit_totals"] = {"sugar": [1, 2]}
        self.assertEqual(s.errors(), [])

    def test_upgrade_targets_are_checked_too(self):
        s = Set()
        s.rule["terms"].append(term("up", "aglycone_ion", "upgrade", form="even", charge_carrier="adduct",
                                    plus=None, minus=None))
        s.rule["preconditions"]["composition"]["units"] = {"Hex": [0, 2]}
        s.rule["preconditions"]["composition"]["unit_totals"] = {"sugar": [0, 2]}
        self.assertIn(("support_target_is_precursor", RULE_FILE, "/terms/4"), s.errors())


class ProvenanceTests(unittest.TestCase):
    def derived_error(self, mutate, origin="as_published"):
        s = Set()
        mutate(s)
        s.rule["parameters_origin"] = origin
        return ("origin_inconsistent", RULE_FILE, "/parameters_origin") in s.errors()

    def test_each_source_raises_the_derived_origin(self):
        cases = {
            "term_source": lambda s: t_of(s, "Y0").update(term_source="tuned_vs20"),
            "floor_source": lambda s: t_of(s, "rda").update(floor_source="tuned_vs20"),
            "ratio_source": lambda s: t_of(s, "rda").update(role="ordering", than=["Y0"], ratio=2.0,
                                                            ratio_source="tuned_vs20"),
            "k_source": lambda s: s.rule["logic"]["all"][2]["k_of"].update(k_source="tuned_vs20"),
            "scaffold entry": lambda s: s.scaffolds["aglycones"]["entries"][1].update(origin="tuned_vs20"),
            "unit row": lambda s: s.units[1].update(origin="tuned_vs20"),
            "tried common adduct": lambda s: s.adducts[4].update(origin="tuned_vs20"),
            "motif": lambda s: s.motifs[0].update(origin="tuned_vs20"),
        }
        for name, mutate in cases.items():
            with self.subTest(name=name):
                self.assertTrue(self.derived_error(mutate))
                self.assertFalse(self.derived_error(mutate, "tuned_vs20"))

    def test_pooled_wins(self):
        def both(s):
            t_of(s, "Y0").update(term_source="tuned_vs20")
            s.rule["logic"]["all"][2]["k_of"].update(k_source="tuned_pooled")
        self.assertTrue(self.derived_error(both, "tuned_vs20"))
        self.assertFalse(self.derived_error(both, "tuned_pooled"))

    def test_unrelated_rows_do_not_count(self):
        for mutate in (lambda s: s.units[2].update(origin="tuned_vs20"),  # Pen: not in composition.units
                       lambda s: s.adducts[5].update(origin="tuned_vs20"),  # [M+Cl]-: not common, not listed
                       lambda s: s.adducts[0].update(origin="tuned_vs20"),  # positive adduct
                       lambda s: s.motifs[2].update(origin="tuned_vs20")):
            self.assertFalse(self.derived_error(mutate))


class ReviewTests(unittest.TestCase):
    def reviewed(self):
        s = Set()
        s.rule.update(author_reviewed=True, source_checked=True)
        res = s.load()
        self.assertEqual({e.code for e in res.errors}, {"review_missing"})
        s.rule["review"] = review_block(res.semantic_sha256(s.rule["rule_id"]))
        self.assertEqual(s.errors(), [])
        return s

    def test_stale_for_each_resolved_dependency(self):
        stale = {
            "residue": lambda s: s.units[0].update(residue="C6H12O6"),
            "scaffold entry added": lambda s: s.scaffolds["aglycones"]["entries"].append(dict(
                SCAFFOLDS["entries"][0], scaffold_id="genistein", name="genistein")),
            "tried adduct common flag": lambda s: s.adducts[5].update(common=True),
            "motif smarts": lambda s: s.motifs[0].update(smarts="c1ccccc1O"),
            "element mass": lambda s: s.elements["elements"].update(O=15.994914619),
            "rule term": lambda s: t_of(s, "Y0").update(min_rel_intensity=0.1),
        }
        for name, mutate in stale.items():
            with self.subTest(name=name):
                s = self.reviewed()
                mutate(s)
                self.assertIn(("review_stale", RULE_FILE, "/review/rule_semantic_sha256"), s.errors())

    def test_not_stale_for_unrelated_edits(self):
        fresh = {
            "another unit": lambda s: s.units.append(unit("feruloyl", "acyl", "C10H8O3")),
            "chemont node": lambda s: s.motifs[0]["chemont"]["broader"].append({"id": "CHEMONTID:0000261",
                                                                                "name": "Phenylpropanoids"}),
            "nesting": lambda s: s.nesting.append({"child": "hydroxycinnamoyl_hexose", "parent": "flavonoid",
                                                   "source": CITE}),
            "non-common adduct alias": lambda s: s.adducts[2]["aliases"].append("M+"),
            "unused element": lambda s: s.elements["elements"].update(P=30.973761629),
            "review fields": lambda s: s.rule["review"].update(run_id="p3cs-other"),
        }
        for name, mutate in fresh.items():
            with self.subTest(name=name):
                s = self.reviewed()
                mutate(s)
                self.assertEqual(s.errors(), [])

    def test_ppm_only_msfinder_tolerance_records_abs_da_zero(self):
        # The evaluator applies a ppm-only MS-FINDER tolerance as {0, X, msfinder}; a review measured at it says so.
        pointer = "/review/product_tolerance/abs_da"
        cases = {
            "msfinder ppm only": ({"abs_da": 0.0, "ppm": 10.0, "source": "msfinder"}, []),
            "msfinder Da": ({"abs_da": 0.01, "ppm": 0.0, "source": "msfinder"}, []),
            "msfinder zero": ({"abs_da": 0.0, "ppm": 0.0, "source": "msfinder"}, [("bad_value", RULE_FILE, pointer)]),
            "p2 abs_da zero": ({"abs_da": 0.0, "ppm": 10.0, "source": "p2"}, [("bad_value", RULE_FILE, pointer)]),
            "negative abs_da": ({"abs_da": -0.01, "ppm": 10.0, "source": "msfinder"}, [("bad_value", RULE_FILE, pointer)]),
            "integer ppm": ({"abs_da": 0.0, "ppm": 10, "source": "msfinder"},
                            [("bad_value", RULE_FILE, pointer),
                             ("bad_type", RULE_FILE, "/review/product_tolerance/ppm")]),
        }
        for name, (tolerance, errors) in cases.items():
            with self.subTest(name=name):
                s = self.reviewed()
                s.rule["review"]["product_tolerance"] = tolerance
                self.assertEqual(s.errors(), sorted(errors, key=lambda e: (e[1], e[2], e[0])))

    def test_engine_semantics_is_hashed(self):
        s = self.reviewed()
        res = s.load()
        rule = res.rules[s.rule["rule_id"]]
        resolved = rf.rule_resolved(rule, res.tables)
        self.assertNotEqual(rf.rule_semantic_sha256(rule.doc, resolved, "msfinder-diagnostic-rules/2"),
                            s.rule["review"]["rule_semantic_sha256"])
        self.assertEqual(sorted(resolved["elements"]["masses"]), ["C", "H", "O"])
        self.assertEqual([a["adduct"] for a in resolved["adducts"]], ["[M-H]-", "[M+HCOO]-"])
        self.assertNotIn("chemont", resolved["motif"])


class DraftAndPinTests(unittest.TestCase):
    def staged(self):
        s = Set()
        s.drop |= {"tables/motifs.json", "tables/motif_nesting.json"}
        s.extra["_motifs/a.json"] = dumps(MOTIFS[:2])
        s.extra["_motifs/b.json"] = dumps(MOTIFS[2:])
        s.extra["_nesting/a.json"] = dumps(NESTING)
        s.extra["_requests/a.json"] = dumps({"units": [unit("feruloyl", "acyl", "C10H8O3")],
                                             "adducts": [adduct("[M+K]+", 1, add="K")]})
        s.extra["_audit/notes.json"] = b"not json at all"
        return s

    def test_draft_merges_staging_files(self):
        s = self.staged()
        s.rule["preconditions"]["composition"]["units"]["feruloyl"] = [0, 1]
        res = s.load()
        self.assertEqual(res.errors, [])
        self.assertEqual(sorted(res.tables.motifs), ["flavonoid", "flavonoid_o_glycoside", "hydroxycinnamoyl_hexose"])
        self.assertIn("feruloyl", res.tables.unit_by_name)
        self.assertIn("[M+K]+", res.tables.adduct_by_name)
        self.assertEqual(res.staging_files, ["_motifs/a.json", "_motifs/b.json", "_nesting/a.json", "_requests/a.json"])

    def test_staging_files_are_not_hashed(self):
        s = self.staged()
        with_staging = s.load().ruleset_sha256
        for p in [p for p in s.extra if p.startswith("_")]:
            del s.extra[p]
        self.assertEqual(s.load().ruleset_sha256, with_staging)

    def test_duplicate_ids_across_files(self):
        s = self.staged()
        s.extra["_motifs/b.json"] = dumps(MOTIFS[1:])
        s.extra["_requests/b.json"] = dumps({"units": [unit("feruloyl", "acyl", "C10H8O3")], "adducts": []})
        errs = s.errors()
        self.assertIn(("duplicate_id", "_motifs/b.json", "/0/motif_id"), errs)
        self.assertIn(("duplicate_id", "_requests/b.json", "/units/0/name"), errs)

    def test_staging_without_draft_is_unexpected(self):
        s = self.staged()
        errs = s.errors(draft=False)
        self.assertIn(("unexpected_file", "_motifs/a.json", ""), errs)
        self.assertIn(("missing_file", "tables/motifs.json", ""), errs)

    def test_staging_shapes(self):
        s = Set()
        s.extra["_motifs/a.json"] = dumps({"motifs": MOTIFS})
        s.extra["_requests/a.json"] = dumps({"units": []})
        errs = s.errors()
        self.assertIn(("bad_type", "_motifs/a.json", ""), errs)
        self.assertIn(("missing_key", "_requests/a.json", "/adducts"), errs)

    def test_version_pin(self):
        s = Set()
        self.assertEqual(s.errors(draft=False), [("ruleset_unpinned", "ruleset.json", "/rule_set_version")])
        sha = s.load().ruleset_sha256
        s.extra[rf.HISTORY_FILE] = f"0.0.0\t{'1' * 64}\t2026-11-01\r\n0.0.1-test\t{sha}\t2026-11-18\n".encode()
        self.assertEqual(s.errors(draft=False), [])
        s.rule["notes"] = "changed"
        self.assertEqual(s.errors(draft=False), [("ruleset_sha_mismatch", "ruleset.json", "/rule_set_version")])
        s.extra[rf.HISTORY_FILE] += f"0.0.1-test\t{sha}\t2026-11-19\n".encode()
        self.assertIn(("history_format", rf.HISTORY_FILE, ""), s.errors(draft=False))

    def test_rules_glob(self):
        s = Set()
        other = core_rule("flavonoid.synthetic.neg", {"all": [{"term": "X"}]},
                          [term("X", "product_ion", "required", formula="C7H3O4")], ["X"])
        s.rules.append(other)
        self.assertEqual(s.load(rules_glob="flavonoid.o_*").errors, [])
        self.assertEqual(s.load(rules_glob="flavonoid.o_*").n_rules, 1)
        self.assertTrue(s.load().errors)


class CliTests(unittest.TestCase):
    def run_cli(self, argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            code = rf.main(argv)
        return code, out.getvalue().splitlines()

    def write(self, root, files):
        for rel, b in files.items():
            p = Path(root) / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b)

    def test_validate_sha_masses_fmt(self):
        with tempfile.TemporaryDirectory() as d:
            s = Set()
            self.write(d, s.files())
            code, lines = self.run_cli(["validate", d, "--draft"])
            self.assertEqual((code, lines), (0, ["OK n_rules=1"]))
            code, lines = self.run_cli(["validate", d])
            self.assertEqual((code, lines[-1]), (1, "ERRORS n=1"))
            self.assertEqual(lines[0], "ruleset_unpinned\truleset.json\t/rule_set_version")
            code, lines = self.run_cli(["sha", d, "--draft"])
            self.assertEqual(code, 0)
            self.assertEqual(lines[0], "ruleset_sha256\t" + s.load().ruleset_sha256)
            self.assertIn(rf.file_sha256(s.files()[RULE_FILE]) + "\t" + RULE_FILE, lines)
            code, lines = self.run_cli(["masses", d, "--rule", "flavonoid.o_glycoside.neg", "--draft"])
            self.assertEqual(code, 0)
            h = next(x for x in lines if x.startswith("H\tapigenin|Hex1|[M-H]-\t"))
            self.assertIn("precursor_mz=431.098371", h)
            self.assertTrue(any("ion=C15H9O5\ttarget_mz=269.045547" in x for x in lines))
            p = Path(d) / "loose.json"
            p.write_bytes(b'{"b": 1, "a": [0.050]}')
            self.assertEqual(self.run_cli(["fmt", str(p)])[0], 0)
            self.assertEqual(p.read_bytes(), b'{\n "a": [\n  0.05\n ],\n "b": 1\n}\n')
            p.write_bytes(b'{"a": 1e5}')
            self.assertEqual(self.run_cli(["fmt", str(p)]), (1, [f"number_form\t{p.as_posix()}\t"]))


if __name__ == "__main__":
    unittest.main()
