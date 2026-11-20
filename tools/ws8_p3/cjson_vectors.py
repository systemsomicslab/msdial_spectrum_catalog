"""WS8 P3 parity vectors for the C# port of the rule-file format (stdlib only; CPython is the oracle).

Writes cjson-vectors.json (WS8-P3 C# design, sections 2.5, 3.1-3.2 and 7.2): the strict reader and canonical JSON
(key order, nesting, strings, integers, every accepted float form and boundary, every refused form), float lexemes
to IEEE bits, fixed-decimal and result-double text, string escaping, NaN canonicalisation, the instrument class,
per-file canonical sha256 and the rule-set (Merkle) sha256 of small synthetic sets, rule_semantic_sha256, formula
to integer pico-dalton masses (elements, ions with the electron, adduct precursors), and loader errors: a synthetic
file set with one defect for every code of rule_format.ERROR_CODES, with the exact sorted error list.

    python -B tools/ws8_p3/cjson_vectors.py write <cjson-vectors.json> [--source-md]
    python -B tools/ws8_p3/cjson_vectors.py check <cjson-vectors.json>

`write` prints the file's sha256; `--source-md` also writes SOURCE.md beside it (generator, catalog commit,
interpreter, sha256; repository-relative paths only). `check` regenerates and compares bytes. The output is
deterministic: the same code gives the same bytes under any hash seed. Everything is synthetic.
"""
from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import json
import math
from pathlib import Path
import platform
import struct
import subprocess
import sys

try:
    from . import rule_format as rf
except ImportError:  # run as a script
    import rule_format as rf

FORMAT = "msfinder-diagnostic-cjson-vectors"
GENERATOR = "tools/ws8_p3/cjson_vectors.py"
HISTORY_DATE = "2026-11-19"


def b64(b: bytes) -> str:
    return base64.b64encode(b).decode("ascii")


def ieee_hex(x: float) -> str:
    return struct.pack(">d", x).hex()


def from_hex(h: str) -> float:
    return struct.unpack(">d", bytes.fromhex(h))[0]


def bits(x: float) -> str:
    """Design 4.11: the raw IEEE-754 pattern, big-endian, 16 lowercase hex digits, every NaN canonicalised."""
    return "7ff8000000000000" if x != x else struct.pack(">d", x).hex()


def fixed(x: float, d: int) -> str:
    """Design 4.11: exact half-even fixed decimals; text whose digits are all zero has no sign."""
    s = f"{x:.{d}f}"
    return s[1:] if s.startswith("-") and set(s[1:]) <= set("0.") else s


def result_double(x: float, d: int):
    return None if not math.isfinite(x) else fixed(x, d)


_UP = {c: c - 32 for c in range(0x61, 0x7B)}  # a-z only


def instrument_class(value):
    """Design 4.9, the reference for both languages."""
    if value is None or value.strip(" \t\r\n") == "":
        return "unknown"
    s = value.translate(_UP)
    if "TOF" in s:
        return "QTOF-type"
    if any(t in s for t in ("ORBITRAP", "FT", "Q EXACTIVE", "Q-EXACTIVE")):
        return "Orbitrap/FT"
    return "other"


def utf16_hex(s: str) -> str:
    return s.encode("utf-16-be", "surrogatepass").hex()


def from_utf16_hex(h: str) -> str:
    return bytes.fromhex(h).decode("utf-16-be", "surrogatepass")


# --------------------------------------------------------------------------------------- strict reader cases
def deep(n, inner=b"1"):
    """An object holding n nested arrays: depth n + 1."""
    return b'{"a":' + b"[" * n + inner + b"]" * n + b"}"


def deep_objects(n):
    """n nested objects: depth n."""
    return b'{"a":' * n + b"1" + b"}" * n


# (name, input bytes). The expectation is computed by rule_format.load_strict.
STRICT_CASES = [
    # accepted: structure, key order, whitespace
    ("empty_object", b"{}"),
    ("key_order_code_points", b'{"b":1,"a":2,"B":3,"_":4,"a0":5,"a ":6,"":7}'),
    ("nested_key_order", b'{"x":{"z":[1,{"b":null,"a":true}],"y":false},"w":[]}'),
    ("whitespace_crlf_tab", b'{\r\n\t"a" : [ 1 , 2 ] ,\r\n "b":{ }\r\n}\r\n'),
    ("trailing_whitespace", b'{"a":1}\n\n  \t'),
    ("depth_128", deep(127)),
    ("depth_128_objects", deep_objects(128)),
    ("duplicate_key_in_sibling_objects_ok", b'{"x":{"a":1},"y":{"a":1}}'),
    # strings
    ("string_escapes_quote_backslash_solidus", b'{"s":"q\\"b\\\\s\\/"}'),
    ("string_unicode_escape_ascii", b'{"s":"\\u0041\\u007e\\u0020"}'),
    ("string_printable_ascii", b'{"s":"' + bytes(c for c in range(0x20, 0x7F) if c not in (0x22, 0x5C))
     + b'\\"\\\\"}'),
    ("key_with_escape_sorts_unescaped", b'{"b":1,"\\u0042":2,"a":3}'),
    ("empty_string_key_and_value", b'{"":""}'),
    # integers
    ("int_zero", b'{"a":0}'),
    ("int_negative", b'{"a":-5}'),
    ("int_max_safe", b'{"a":9007199254740991}'),
    ("int_min_safe", b'{"a":-9007199254740991}'),
    ("int_16_digits", b'{"a":1000000000000000}'),
    ("int_list", b'{"a":[1,-1,10,123456789]}'),
    # floats: every accepted form and boundary
    ("float_zero", b'{"a":0.0}'),
    ("float_zero_three_places", b'{"a":0.000}'),
    ("float_lower_bound", b'{"a":0.0001}'),
    ("float_lower_bound_trailing_zero", b'{"a":0.00010}'),
    ("float_negative_lower_bound", b'{"a":-0.0001}'),
    ("float_trailing_zeros", b'{"a":[0.050,5.00,10.0,1.50,100.000]}'),
    ("float_one_with_zeros", b'{"a":1.00000000000000}'),
    ("float_many_trailing_zeros", b'{"a":0.100000000000000000000}'),
    ("float_15_digits_k18", b'{"a":0.000123456789012345}'),
    ("float_15_digits", b'{"a":[1.23456789012345,12345678901234.5,123456789012345.0]}'),
    ("float_upper_range", b'{"a":[999999999999999.0,99999999999999.9,100000000000000.0]}'),
    ("float_masses", b'{"a":[1.007825032,0.000549,34.96885268,12.0,15.99491462,2.675,0.1,0.3]}'),
    ("float_negative", b'{"a":[-1.5,-0.05,-123.456]}'),
    ("mixed_document", b'{"rule":{"terms":[{"min_rel_intensity":0.05,"than":null,"role":"required"}],'
                       b'"version":1,"measure":true},"format":"x"}'),
    # refused: bytes and grammar
    ("bom", b"\xef\xbb\xbf{}"),
    ("invalid_utf8_ff", b'{"a":"\xff"}'),
    ("invalid_utf8_overlong", b'{"a":"\xc0\xaf"}'),
    ("invalid_utf8_lone_continuation", b'{"a":"\x80"}'),
    ("invalid_utf8_truncated", b'{"a":"\xe3\x81"}'),
    ("invalid_utf8_encoded_surrogate", b'{"a":"\xed\xa0\x80"}'),
    ("syntax_trailing_comma", b'{"a":1,}'),
    ("syntax_trailing_comma_array", b'{"a":[1,]}'),
    ("syntax_comment", b'{"a":1 /* c */}'),
    ("syntax_line_comment", b'{"a":1 // c\n}'),
    ("syntax_single_quotes", b"{'a':1}"),
    ("syntax_unquoted_key", b"{a:1}"),
    ("syntax_hex", b'{"a":0x10}'),
    ("syntax_leading_zero_int", b'{"a":01}'),
    ("syntax_leading_zero_float", b'{"a":01.5}'),
    ("syntax_plus_sign", b'{"a":+1}'),
    ("syntax_bare_fraction", b'{"a":.5}'),
    ("syntax_trailing_point", b'{"a":1.}'),
    ("syntax_minus_bare_fraction", b'{"a":-.5}'),
    ("syntax_raw_newline_in_string", b'{"a":"x\ny"}'),
    ("syntax_raw_tab_in_string", b'{"a":"x\ty"}'),
    ("syntax_raw_nul_in_string", b'{"a":"x\x00y"}'),
    ("syntax_form_feed_whitespace", b'{\x0c"a":1}'),
    ("syntax_vertical_tab_whitespace", b'{"a":\x0b1}'),
    ("syntax_nbsp_whitespace", b'{"a":\xc2\xa01}'),
    ("syntax_empty_input", b""),
    ("syntax_whitespace_only", b" \r\n"),
    ("syntax_unclosed", b'{"a":1'),
    ("syntax_two_documents", b"{}{}"),
    ("syntax_bad_escape", b'{"a":"\\x41"}'),
    ("syntax_short_unicode_escape", b'{"a":"\\u12"}'),
    ("syntax_true_case", b'{"a":True}'),
    ("top_level_array", b"[]"),
    ("top_level_string", b'"x"'),
    ("top_level_number", b"1"),
    ("top_level_null", b"null"),
    ("duplicate_key", b'{"a":1,"a":2}'),
    ("duplicate_key_nested", b'{"x":{"a":1,"a":1}}'),
    ("duplicate_key_by_escape", b'{"a":1,"\\u0061":2}'),
    ("nonfinite_nan", b'{"a":NaN}'),
    ("nonfinite_infinity", b'{"a":Infinity}'),
    ("nonfinite_negative_infinity", b'{"a":-Infinity}'),
    ("depth_129", deep(128)),
    ("depth_129_objects", deep_objects(129)),
    ("depth_5000", deep(5000)),
    # refused: integers
    ("int_17_digits", b'{"a":12345678901234567}'),
    ("int_4301_digits", b'{"a":' + b"9" * 4301 + b"}"),
    ("int_above_max_safe", b'{"a":9007199254740992}'),
    ("int_below_min_safe", b'{"a":-9007199254740992}'),
    ("int_16_digits_above_max_safe", b'{"a":9999999999999999}'),
    ("int_negative_zero", b'{"a":-0}'),
    # refused: floats, in the order the checks run
    ("float_exponent", b'{"a":1e5}'),
    ("float_exponent_upper", b'{"a":1E5}'),
    ("float_fraction_exponent", b'{"a":1.5e3}'),
    ("float_negative_exponent", b'{"a":1.0e-4}'),
    ("float_signed_exponent", b'{"a":1.5E+3}'),
    ("float_huge_exponent", b'{"a":1e400}'),
    ("float_16_digits", b'{"a":0.1234567890123456}'),
    ("float_17_digits", b'{"a":0.30000000000000004}'),
    ("float_16_digits_large", b'{"a":999999999999999.5}'),
    ("float_negative_zero", b'{"a":-0.0}'),
    ("float_negative_zero_places", b'{"a":-0.000}'),
    ("float_below_range", b'{"a":0.00009999}'),
    ("float_below_range_k18", b'{"a":0.000099999999999999}'),
    ("float_tiny", b'{"a":0.00001}'),
    ("float_negative_tiny", b'{"a":-0.00001}'),
    ("float_at_upper_bound", b'{"a":1000000000000000.0}'),
    ("float_above_upper_bound", b'{"a":10000000000000000.0}'),
    ("float_17_digits_above_upper_bound", b'{"a":12345678901234567.0}'),
    # refused: strings after unescaping
    ("control_char_escaped_7f", b'{"a":"\\u007f"}'),
    ("control_char_raw_7f", b'{"a":"\x7f"}'),
    ("control_char_escaped_nul", b'{"a":"\\u0000"}'),
    ("control_char_escaped_tab", b'{"a":"\\t"}'),
    ("control_char_escaped_newline", b'{"a":"x\\ny"}'),
    ("control_char_in_key", b'{"\\u001f":1}'),
    ("non_ascii_escaped_80", b'{"a":"\\u0080"}'),
    ("non_ascii_raw_e_acute", '{"a":"é"}'.encode("utf-8")),
    ("non_ascii_surrogate_pair", b'{"a":"\\ud83d\\ude00"}'),
    ("non_ascii_lone_high_surrogate", b'{"a":"\\ud800"}'),
    ("non_ascii_lone_low_surrogate", b'{"a":"\\udfff"}'),
    ("non_ascii_raw_emoji", b'{"a":"\xf0\x9f\x98\x80"}'),
    ("non_ascii_in_key", b'{"\\u00e9":1}'),
    ("first_offending_code_point_decides", b'{"a":"\\u00e9\\u0001"}'),
    ("key_checked_before_value", b'{"\\u0001":"\\u00e9"}'),
]


def strict_vectors():
    out, names = [], set()
    for name, raw in STRICT_CASES:
        assert name not in names, name
        names.add(name)
        case = {"name": name, "input_b64": b64(raw)}
        try:
            obj = rf.load_strict(raw)
        except rf.RuleFormatError as e:
            case.update(expect="error", code=e.code)
        else:
            canon = rf.canonical_bytes(obj)
            case.update(expect="ok", canonical=canon.decode("ascii"), sha256=hashlib.sha256(canon).hexdigest())
            assert case["sha256"] == rf.file_sha256(raw)
        out.append(case)
    return out


# --------------------------------------------------------------------------------------------- float lexemes
FLOAT_LEXEMES = [
    "0.0", "0.000", "0.0001", "0.00010", "-0.0001", "0.05", "0.050", "5.0", "5.00", "10.0", "1.5", "-1.5", "0.1",
    "0.2", "0.3", "0.7", "2.675", "1.007825032", "0.000549", "12.0", "14.003074", "15.99491462", "22.98976928",
    "30.97376163", "31.972071", "34.96885268", "38.96370668", "162.05282342", "449.107837492", "0.01", "0.005",
    "10.0", "100.0", "1.00000000000000", "0.100000000000000000000", "0.000123456789012345", "1.23456789012345",
    "12345678901234.5", "123456789012345.0", "999999999999999.0", "99999999999999.9", "100000000000000.0",
    "0.999999999999999", "9.99999999999999", "0.000100000000000001", "450359962737.049", "0.0009765625",
    "-123.456", "-999999999999999.0", "-0.000549",
]


def lexeme_parts(lex):
    """Design 3.1: strip trailing fraction zeros; D is the remaining digit integer (signed), k the fraction digits."""
    neg = lex.startswith("-")
    ip, fp = lex.lstrip("-").split(".")
    fp = fp.rstrip("0")
    d = int(ip + fp)
    return (-d if neg else d), len(fp)


def float_vectors():
    out, seen = [], set()
    for lex in FLOAT_LEXEMES:
        if lex in seen:
            continue
        seen.add(lex)
        v = rf.load_strict(b'{"a":' + lex.encode() + b"}")["a"]
        assert isinstance(v, float), lex
        mantissa, scale = lexeme_parts(lex)
        assert abs(mantissa) < 10 ** 15 and scale <= 18, lex
        assert mantissa / 10 ** scale == v or (mantissa == 0 and v == 0.0), lex
        canonical = rf.canonical_bytes([v]).decode("ascii")[1:-1]
        out.append({"lexeme": lex, "ieee_hex": ieee_hex(v), "mantissa": mantissa, "scale": scale,
                    "canonical": canonical})
    return out


# ------------------------------------------------------------------------------- fixed decimals and doubles
FIXED_CASES = [
    (0.0, 6), (-0.0, 6), (1e-9, 6), (-1e-9, 6), (-4e-7, 6), (-6e-7, 6), (5e-7, 6), (1.5e-6, 6), (2.5e-6, 6),
    (1e-7, 6), (0.5, 0), (1.5, 0), (2.5, 0), (-0.5, 0), (-1.5, 0), (0.125, 2), (0.375, 2), (-0.125, 2),
    (0.0625, 3), (0.1875, 3), (-0.0625, 3), (1234.5678, 3), (12.3455, 3), (-12.3455, 3), (449.107837492, 6),
    (341.087806124, 6), (153.01823464, 6), (0.000549, 6), (1.007825032, 6), (123.0000005, 6), (2.0000005, 6),
    (1e15, 6), (9007199254740991.0, 6), (0.1, 6), (0.3, 3), (1e-4, 6), (100.0, 1), (-100.0, 1), (5e-324, 6),
    (1.7976931348623157e308, 0), (0.1, 20), (1.0 / 3.0, 6), (-2.0 / 3.0, 3),
]
RESULT_CASES = [
    ("7ff8000000000000", 6), ("fff8000000000000", 6), ("7ff0000000000001", 6), ("7ff0000000000000", 6),
    ("fff0000000000000", 3), ("8000000000000000", 6), ("0000000000000000", 3), ("3ff0000000000000", 6),
    ("bff8000000000000", 0), ("4076c1b9f8c0d3a5", 6),
]
BITS_CASES = [
    "7ff8000000000000", "fff8000000000000", "7ff0000000000001", "7fffffffffffffff", "ffffffffffffffff",
    "fff0000000000001", "7ff4000000000000", "8000000000000000", "0000000000000000", "3ff0000000000000",
    "bff0000000000000", "7ff0000000000000", "fff0000000000000", "0000000000000001", "800fffffffffffff",
    "7fefffffffffffff", "4076c1b9f8c0d3a5",
]


def fixed_vectors():
    return [{"ieee_hex": ieee_hex(x), "d": d, "text": fixed(x, d)} for x, d in FIXED_CASES]


def result_vectors():
    out = []
    for h, d in RESULT_CASES:
        x = from_hex(h)
        out.append({"ieee_hex": h, "d": d, "text": result_double(x, d)})
    return out


def bits_vectors():
    return [{"ieee_hex_in": h, "bits": bits(from_hex(h))} for h in BITS_CASES]


# ------------------------------------------------------------------------------------------ string escaping
STRING_CASES = [
    '"', "\\", "/", "\b", "\f", "\n", "\r", "\t", "\x00", "\x01", "\x1f", "\x7f", "\x80", "é", " ",
    " ", "﻿", "￿", "\U0001f600", "\ud800", "\udfff", "\udc00\ud800", " ~", "plain ASCII text",
    'a"b\\c\ndé\U0001f600\x7f',
    "".join(chr(c) for c in range(0x20, 0x7F)),
]


def string_vectors():
    out = []
    for s in STRING_CASES:
        j = json.dumps(s, ensure_ascii=True)
        assert all(0x20 <= ord(c) < 0x7F for c in j)
        out.append({"text_utf16_hex": utf16_hex(s), "json": j})
    return out


INSTRUMENT_CASES = [
    None, "", " ", " \t\r\n", "\x0b", " ", "QTOF", "qtof", "Q-TOF", "Agilent 6530 Q-TOF", "LC-ESI-QTOF",
    "tOf", "maXis Impact", "Orbitrap", "LTQ Orbitrap XL", "orbitrap fusion", "LC-ESI-QFT", "FT-ICR", "Q Exactive",
    "Q Exactive Plus", "q-exactive", "Q Exactive HF", "QExactive", "LTQ", "QqQ", "Ion Trap", "LC-ESI-IT-TOF",
    "Orbitrap ID-X with TOF", "Microsoft", "LIFT", "IT-FT", "LTQ FT Ultra", "ｑｔｏｆ",
    "toﬀ", "orbitrapé", "éTOF", "Q Exactive",
]


def instrument_vectors():
    return [{"text_utf16_hex": None if t is None else utf16_hex(t), "class": instrument_class(t)}
            for t in INSTRUMENT_CASES]


# ---------------------------------------------------------------------------------------- synthetic rule sets
ELEMENTS = {"format": "msfinder-diagnostic-elements", "format_version": 1, "electron": 0.000549,
            "elements": {"C": 12.0, "Cl": 34.96885268, "H": 1.007825032, "K": 38.96370668, "N": 14.003074,
                         "Na": 22.98976928, "O": 15.99491462, "P": 30.97376163, "S": 31.972071}}
CITE = [{"doi": "10.1000/synthetic.1", "locator": "Table 1"}]
REF = [{"reference": "synthetic adduct table", "locator": "line 2"}]
VERSION = "0.0.1-test"


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


def core_rule(rule_id, logic, terms, claim, core, mode, adducts, measure=True, motif_id="flavonoid_o_glycoside"):
    r = base_rule()
    r.update(rule_id=rule_id, measure=measure, terms=terms, logic=logic)
    r["preconditions"].update(ion_mode=mode, adducts=list(adducts))
    r["preconditions"]["composition"] = {"core": core, "scaffold_list": None, "units": {}, "unit_totals": {}}
    r["claim"].update(claim_term_ids=claim, motif_id=motif_id)
    return r


def wrap(fmt, key, rows):
    return {"format": fmt, "format_version": 1, key: rows}


class SynSet:
    """A synthetic rule set as plain objects; files() renders it in storage form, with a valid version pin."""

    def __init__(self, rules=None):
        self.ruleset = {"format": "msfinder-diagnostic-ruleset", "format_version": 1, "rule_set_version": VERSION}
        self.elements = copy.deepcopy(ELEMENTS)
        self.adducts = copy.deepcopy(ADDUCTS)
        self.units = copy.deepcopy(UNITS)
        self.motifs = copy.deepcopy(MOTIFS)
        self.nesting = copy.deepcopy(NESTING)
        self.scaffolds = {"aglycones": copy.deepcopy(SCAFFOLDS)}
        self.rules = [base_rule()] if rules is None else rules
        self.extra = {}
        self.drop = set()
        self.history = None  # None: a valid row for the set's own sha256

    @property
    def rule(self):
        return self.rules[0]

    def json_files(self):
        d = rf.storage_bytes
        f = {rf.RULESET_FILE: d(self.ruleset), rf.ELEMENTS_FILE: d(self.elements),
             rf.ADDUCTS_FILE: d(wrap("msfinder-diagnostic-adducts", "adducts", self.adducts)),
             rf.UNITS_FILE: d(wrap("msfinder-diagnostic-units", "units", self.units)),
             rf.MOTIFS_FILE: d(wrap("msfinder-diagnostic-motifs", "motifs", self.motifs)),
             rf.NESTING_FILE: d(wrap("msfinder-diagnostic-motif-nesting", "nesting", self.nesting))}
        for lid, doc in self.scaffolds.items():
            f[f"scaffolds/{lid}.json"] = d(doc)
        for r in self.rules:
            f[f"rules/{r['rule_id']}.json"] = d(r)
        f.update(self.extra)
        for p in self.drop:
            f.pop(p, None)
        return f

    def files(self):
        f = self.json_files()
        if self.history is not None:
            f[rf.HISTORY_FILE] = self.history
        else:
            sha = rf.load_files(f, draft=True).ruleset_sha256
            if sha is not None:
                f[rf.HISTORY_FILE] = f"{VERSION}\t{sha}\t{HISTORY_DATE}\n".encode("ascii")
        return dict(sorted(f.items()))

    def load(self, draft=False):
        return rf.load_files(self.files(), draft=draft)


RULE_FILE = "rules/flavonoid.o_glycoside.neg.json"


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


def deep_logic(n):
    node = {"term": "Y0"}
    for _ in range(n):
        node = {"all": [node]}
    return node


def review_block(sha):
    return {"by": "author", "at": "2026-12-01", "basis": "benchmark", "run_id": "p3cs-20261201-abc",
            "rules_precision_sha256": "a" * 64, "review_proposal_sha256": "b" * 64, "rule_semantic_sha256": sha,
            "instrument_classes": ["QTOF-type"], "product_tolerance": {"abs_da": 0.005, "ppm": 10.0, "source": "p2"},
            "precursor_tolerance_da": 0.01}


def reviewed_with(s, **fields):
    """Give s.rule a review whose rule_semantic_sha256 is current, after applying fields."""
    s.rule.update({"author_reviewed": True, "source_checked": True, **fields})
    s.rule["review"] = review_block("0" * 64)
    s.rule["review"]["rule_semantic_sha256"] = rf.load_files(s.json_files(), draft=True).semantic_sha256(
        s.rule["rule_id"])


# One synthetic defect per code of rule_format.ERROR_CODES (non-draft load, version pin checked).
DEFECTS = {
    "bom": raw_rule(b"\xef\xbb\xbf{}"),
    "invalid_utf8": raw_rule(b'{"a":"\xff"}'),
    "syntax": raw_rule(b'{"a":1,}'),
    "top_level_not_object": raw_rule(b"[]"),
    "duplicate_key": raw_rule(b'{"a":1,"a":2}'),
    "nonfinite": raw_rule(b'{"a":NaN}'),
    "nesting_depth": raw_rule(deep(128)),
    "number_range": raw_rule(b'{"a":12345678901234567}'),
    "negative_zero": raw_rule(b'{"a":-0}'),
    "number_form": raw_rule(b'{"a":1e5}'),
    "float_precision": raw_rule(b'{"a":0.1234567890123456}'),
    "control_char": raw_rule(b'{"a":"\\u007f"}'),
    "non_ascii": raw_rule(b'{"a":"\\u00e9"}'),
    "missing_file": lambda s: s.drop.add(rf.NESTING_FILE),
    "unexpected_file": lambda s: s.extra.update({"tables/extra.json": b"{}\n"}),
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
    "file_name_mismatch": lambda s: s.extra.update(
        {RULE_FILE: rf.storage_bytes(dict(base_rule(), rule_id="flavonoid.other.neg"))}),
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
    "logic_depth": set_rule(["logic", "all", 0], deep_logic(16)),
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
    "history_format": lambda s: setattr(s, "history", f"{VERSION}\tnot-a-sha\t{HISTORY_DATE}\n".encode("ascii")),
    "ruleset_unpinned": lambda s: setattr(s, "history", b""),
    "ruleset_sha_mismatch": lambda s: setattr(s, "history", f"{VERSION}\t{'0' * 64}\t{HISTORY_DATE}\n".encode()),
}


def patch_of(base, files):
    """{path: b64 or None}: what turns the base file map into this one (None removes a path)."""
    out = {p: b64(b) for p, b in sorted(files.items()) if base.get(p) != b}
    out.update({p: None for p in sorted(base) if p not in files})
    return dict(sorted(out.items()))


def loader_vectors():
    base = SynSet().files()
    res = rf.load_files(base, draft=False)
    assert res.errors == [], res.errors
    cases = [{"name": "clean", "code": None, "patch": {}, "errors": []}]
    missing = sorted(set(rf.ERROR_CODES) - set(DEFECTS))
    assert not missing, missing
    for code in sorted(DEFECTS):
        s = SynSet()
        DEFECTS[code](s)
        files = s.files()
        errors = rf.load_files(files, draft=False).errors
        listed = [[e.code, e.file, e.pointer] for e in errors]
        assert code in {e.code for e in errors}, (code, listed)
        cases.append({"name": code, "code": code, "patch": patch_of(base, files), "errors": listed})
    return {"call": "rule_format.load_files(files, draft=False): every *.json at a layout path is hashed, the "
                    "motif and nesting tables are required, staging folders are unexpected, and the version pin "
                    "is checked against ruleset-history.tsv",
            "base_files": {p: b64(b) for p, b in base.items()}, "cases": cases}


def rich_set():
    """A clean set with three rules: a reviewed o-glycoside with an ion series, a flavylium [M]+ rule, and an
    unmeasured caffeoyl hexose with a complementary pair and a [2M-H]- adduct."""
    s = SynSet()
    s.adducts.append(adduct("[2M-H]-", -1, remove="H", common=False, multiplier=2))
    r = s.rule
    r["terms"].append(term("ys", "ion_series", "upgrade", unit_classes=["sugar"], charge_carrier="proton",
                           form="even", plus=None, minus=None))
    flav = core_rule("anthocyanin.glucoside.pos", {"all": [{"term": "a"}, {"term": "b"}]},
                     [term("a", "product_ion", "required", formula="C15H11O6"),
                      term("b", "loss_from_precursor", "required", formula="C6H10O5")], ["a", "b"],
                     core="C21H21O11", mode="pos", adducts=("[M]+",))
    caf = core_rule("hydroxycinnamoyl.hexose_ester.neg", {"all": [{"term": "loss_hex"}, {"term": "caffeate"}]},
                    [term("loss_hex", "loss_from_precursor", "required", formula="C6H10O5"),
                     term("caffeate", "product_ion", "required", formula="C9H7O4")],
                    ["loss_hex", "caffeate"], core="C15H18O9", mode="neg", adducts=("[M-H]-", "[2M-H]-"),
                    measure=False, motif_id="hydroxycinnamoyl_hexose")
    caf["complementary"] = [["loss_hex", "caffeate"]]
    s.rules += [flav, caf]
    reviewed_with(s)
    return s


def ruleset_vectors():
    out = []
    for name, s in (("loader_base", SynSet()), ("three_rules", rich_set())):
        files = s.files()
        res = rf.load_files(files, draft=False)
        assert res.errors == [], (name, res.errors)
        out.append({
            "name": name, "files": {p: b64(b) for p, b in files.items()}, "errors": [],
            "rule_set_version": res.rule_set_version, "ruleset_sha256": res.ruleset_sha256,
            "file_sha256": dict(sorted(res.file_hashes.items())),
            "rules": {rid: {"rule_semantic_sha256": res.semantic_sha256(rid),
                            "resolved": rf.rule_resolved(res.rules[rid], res.tables)} for rid in sorted(res.rules)},
        })
    return out


def semantic_vectors():
    s = rich_set()
    res = rf.load_files(s.files(), draft=False)
    out = []
    for rid in sorted(res.rules):
        rule = res.rules[rid].doc
        resolved = rf.rule_resolved(res.rules[rid], res.tables)
        out.append({"name": rid, "rule": rule, "resolved": resolved, "engine_semantics": rf.ENGINE_SEMANTICS,
                    "sha256": rf.rule_semantic_sha256(rule, resolved)})
    rid = "flavonoid.o_glycoside.neg"
    rule = res.rules[rid].doc
    resolved = rf.rule_resolved(res.rules[rid], res.tables)
    unreviewed = dict(rule, author_reviewed=False, review=None)
    same = rf.rule_semantic_sha256(rule, resolved)
    assert rf.rule_semantic_sha256(unreviewed, resolved) == same
    out.append({"name": "review_fields_ignored", "rule": unreviewed, "resolved": resolved,
                "engine_semantics": rf.ENGINE_SEMANTICS, "sha256": same})
    out.append({"name": "engine_semantics_bumped", "rule": rule, "resolved": resolved,
                "engine_semantics": "msfinder-diagnostic-rules/2",
                "sha256": rf.rule_semantic_sha256(rule, resolved, "msfinder-diagnostic-rules/2")})
    no_motif = dict(resolved, motif=None, scaffold_list=None)
    out.append({"name": "null_motif_and_scaffold_list", "rule": rule, "resolved": no_motif,
                "engine_semantics": rf.ENGINE_SEMANTICS, "sha256": rf.rule_semantic_sha256(rule, no_motif)})
    for case in out:
        assert rf.rule_semantic_sha256(case["rule"], case["resolved"], case["engine_semantics"]) == case["sha256"]
    return out


def file_hash_vectors():
    """Per-file canonical sha256 and the rule-set sha256 over small file maps (every given file is hashed)."""
    ruleset = b'{"format":"msfinder-diagnostic-ruleset","format_version":1,"rule_set_version":"0.0.1"}'
    ruleset_crlf = (b'{\r\n  "rule_set_version" : "0.0.1",\r\n\t"format_version": 1,\r\n'
                    b'  "format": "msfinder-diagnostic-ruleset"\r\n}\r\n')
    elements = rf.storage_bytes(ELEMENTS)
    elements_compact = rf.canonical_bytes(ELEMENTS)
    floats = b'{"x":[0.050,5.00,10.0,1.0078250320],"y":{"b":1,"a":null}}'
    floats_canon = b'{"y":{"a":null,"b":1},"x":[0.05,5.0,10.0,1.007825032]}'
    sets = [
        ("one_file", {"ruleset.json": ruleset}),
        ("one_file_crlf_reordered", {"ruleset.json": ruleset_crlf}),
        ("two_files", {"ruleset.json": ruleset, "tables/elements.json": elements}),
        ("two_files_other_layout", {"ruleset.json": ruleset_crlf, "tables/elements.json": elements_compact}),
        ("float_forms", {"rules/x.json": floats}),
        ("float_forms_canonical", {"rules/x.json": floats_canon}),
        ("path_order", {"rules/b.json": b'{"b":1}', "rules/B.json": b'{"B":1}', "rules/a.json": b'{"a":1}',
                        "scaffolds/x.json": b'{"list_id":"x"}', "tables/elements.json": elements,
                        "ruleset.json": ruleset}),
    ]
    out = []
    for name, files in sets:
        out.append({"name": name, "files": {p: b64(b) for p, b in sorted(files.items())},
                    "file_sha256": {p: rf.file_sha256(b) for p, b in sorted(files.items())},
                    "ruleset_sha256": rf.ruleset_sha256(files)})
    by = {c["name"]: c for c in out}
    assert by["one_file"]["ruleset_sha256"] == by["one_file_crlf_reordered"]["ruleset_sha256"]
    assert by["two_files"]["ruleset_sha256"] == by["two_files_other_layout"]["ruleset_sha256"]
    assert by["float_forms"]["ruleset_sha256"] == by["float_forms_canonical"]["ruleset_sha256"]
    return out


# -------------------------------------------------------------------------------------------------- masses
LEXEME_PDA = ["12.0", "1.007825032", "14.003074", "15.99491462", "22.98976928", "30.97376163", "31.972071",
              "34.96885268", "38.96370668", "0.000549", "0.0001", "1.000000000001", "999.999999999999",
              "0.000100000000001", "1.0000000000001", "1.0078250321234", "0.0005485799"]
FORMULAS = ["C", "H", "N", "O", "P", "S", "Na", "K", "Cl", "H2O", "CO2", "CHO2", "C2H3O2", "H4N", "O3S", "ClH",
            "C6H10O5", "C6H10O4", "C5H8O4", "C6H8O6", "C9H6O3", "C9H7O4", "C7H5O4", "C7H3O4", "C15H10O5",
            "C15H11O6", "C21H20O10", "C21H21O11", "C15H18O9", "C7H12NO9S2", "C4H10N", "C7H15NO3", "C700", "C10H20N2O4P2S2", "C99H99Cl9K9N9Na9O99P9S9"]
FORMULA_ERRORS = ["C1H4", "C0", "c6", "C1000", "C06", "", "C6 H6", "C6H6+", "HCl", "NH4", "H2OC", "CH4C", "OH",
                  "C6H10O5Br", "Xx", "D2O", "C6H10O5-"]
ION_FORMULAS = ["H", "C7H5O4", "C7H3O4", "C9H7O4", "C15H11O5", "C15H9O5", "C21H21O11", "C4H10N", "HO4S",
                "C99H99Cl9K9N9Na9O99P9S9"]
PRECURSOR_ADDUCTS = [
    adduct("[M+H]+", 1, add="H"), adduct("[M+Na]+", 1, add="Na"), adduct("[M+NH4]+", 1, add="H4N"),
    adduct("[M+K]+", 1, add="K"), adduct("[M+H-H2O]+", 1, add="H", remove="H2O", common=False),
    adduct("[M+H-2H2O]+", 1, add="H", remove="H4O2", common=False), adduct("[M]+", 1, common=False),
    adduct("[M-H]-", -1, remove="H"), adduct("[M+HCOO]-", -1, add="CHO2"), adduct("[M+CH3COO]-", -1, add="C2H3O2"),
    adduct("[M+Cl]-", -1, add="Cl"), adduct("[2M-H]-", -1, remove="H", common=False, multiplier=2),
    adduct("[2M+Na]+", 1, add="Na", common=False, multiplier=2),
]
PRECURSOR_NEUTRALS = ["C15H10O5", "C21H20O10", "C21H21O11", "C15H18O9", "C7H15NO3", "C2H5NO2", "C27H30O16"]


def mass_vectors():
    errors = []
    el = rf.load_elements(rf.Ctx(errors, rf.ELEMENTS_FILE), ELEMENTS)
    assert el is not None and not errors, errors
    lex = []
    for t in LEXEME_PDA:
        try:
            lex.append({"lexeme": t, "pda": rf.lexeme_pda(rf.load_strict(b'{"a":' + t.encode() + b"}")["a"])})
        except rf.RuleFormatError as e:
            lex.append({"lexeme": t, "code": e.code})
    formulas = []
    for f in FORMULAS:
        comp = rf.parse_formula(f, el.pda)
        pda = el.mass_pda(comp)
        assert 0 < pda < rf.MASS_LIMIT_PDA, f
        formulas.append({"formula": f, "pda": pda, "mass_ieee_hex": ieee_hex(rf.mz_of(pda)),
                         "mass_text": fixed(rf.mz_of(pda), 6)})
    ferr = []
    for f in FORMULA_ERRORS:
        try:
            rf.parse_formula(f, el.pda)
        except rf.RuleFormatError as e:
            ferr.append({"formula": f, "code": e.code})
        else:
            raise AssertionError(f)
    ions = []
    for f in ION_FORMULAS:
        comp = rf.parse_formula(f, el.pda)
        for z in (1, -1):
            pda = el.ion_pda(comp, z)
            ions.append({"formula": f, "charge": z, "pda": pda, "mz_ieee_hex": ieee_hex(rf.mz_of(pda)),
                         "mz_text": fixed(rf.mz_of(pda), 6)})
    prec = []
    for a_row in PRECURSOR_ADDUCTS:
        c = rf.Ctx(errors, rf.ADDUCTS_FILE)
        a = rf.v_adduct_row(c, a_row, "", el)
        assert a is not None and not errors, (a_row, errors)
        for f in PRECURSOR_NEUTRALS:
            m = el.mass_pda(rf.parse_formula(f, el.pda))
            pda = a.multiplier * m + el.mass_pda(a.add) - el.mass_pda(a.remove) - a.charge * el.electron_pda
            prec.append({"neutral": f, "adduct": a.name, "multiplier": a.multiplier, "charge": a.charge,
                         "add": a_row["add"], "remove": a_row["remove"], "neutral_pda": m, "precursor_pda": pda,
                         "mz_ieee_hex": ieee_hex(rf.mz_of(pda)), "mz_text": fixed(rf.mz_of(pda), 6)})
    return {"elements": ELEMENTS, "electron_pda": el.electron_pda, "element_pda": dict(sorted(el.pda.items())),
            "lexeme_pda": lex, "formula_pda": formulas, "formula_errors": ferr, "ions": ions,
            "adduct_precursors": prec}


# ------------------------------------------------------------------------------------------------ the file
SECTIONS = {
    "strict_json": "input bytes -> ok with canonical text and sha256 of it, or the reader's error code",
    "float_lexeme": "accepted float lexeme -> IEEE bits; mantissa / 10^scale is the one division of design 3.1",
    "fixed_decimal": "double, decimals -> exact half-even text; all-zero digits carry no sign (design 4.11)",
    "result_double": "double, decimals -> text, or null when not finite (design 4.11)",
    "string_escape": "UTF-16 code units -> Python json.dumps(ensure_ascii=True) (design 3.2)",
    "ieee_bits": "IEEE pattern -> bits() with every NaN canonicalised (design 4.11)",
    "instrument_class": "UTF-16 code units or null -> the ASCII-only instrument class (design 4.9)",
    "file_hashes": "file map -> per-file canonical sha256 and ruleset_sha256; every given file is hashed",
    "ruleset_load": "clean synthetic sets -> ruleset, file and rule semantic sha256 with each rule's resolved set",
    "rule_semantic": "rule, resolved, engine_semantics -> rule_semantic_sha256",
    "masses": "lexeme and formula -> integer pDa; ions and adduct precursors with the electron; refused formulas",
    "loader_errors": "base_files updated by each patch (null removes a path) -> the exact sorted error list "
                     "[code, file, RFC 6901 pointer]; one case per error code",
}


def generate() -> dict:
    return {
        "format": FORMAT, "format_version": 1, "generator": GENERATOR,
        "oracle": "CPython json, struct and hashlib through tools/ws8_p3/rule_format.py",
        "engine_semantics": rf.ENGINE_SEMANTICS, "sections": SECTIONS,
        "error_codes": dict(sorted(rf.ERROR_CODES.items())),
        "strict_json": strict_vectors(), "float_lexeme": float_vectors(), "fixed_decimal": fixed_vectors(),
        "result_double": result_vectors(), "string_escape": string_vectors(), "ieee_bits": bits_vectors(),
        "instrument_class": instrument_vectors(), "file_hashes": file_hash_vectors(),
        "ruleset_load": ruleset_vectors(), "rule_semantic": semantic_vectors(), "masses": mass_vectors(),
        "loader_errors": loader_vectors(),
    }


def render() -> bytes:
    data = rf.storage_bytes(generate())
    rf.load_strict(data)  # the C# strict reader can read the vector file itself
    return data


def git_state():
    repo = Path(__file__).resolve().parents[2]

    def run(*args):
        try:
            out = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.SubprocessError):
            return None
        return out.stdout.strip() if out.returncode == 0 else None

    status = run("status", "--porcelain", "--", "tools/ws8_p3")
    return run("rev-parse", "HEAD") or "unknown", (bool(status) if status is not None else None)


def source_md(sha):
    commit, dirty = git_state()
    return (
        "# cjson-vectors.json\n\n"
        "Generated by `tools/ws8_p3/cjson_vectors.py` in `systemsomicslab/msdial_spectrum_catalog` "
        f"(commit `{commit}`{', tools/ws8_p3 modified' if dirty else ''}), with "
        f"{platform.python_implementation()} {platform.python_version()} as the oracle.\n\n"
        f"- sha256 of `cjson-vectors.json`: `{sha}`\n"
        f"- engine semantics: `{rf.ENGINE_SEMANTICS}`\n"
        "- regenerate from the catalog repository root: "
        "`python -B tools/ws8_p3/cjson_vectors.py write <out>/cjson-vectors.json --source-md`\n"
        "- check: `python -B tools/ws8_p3/cjson_vectors.py check <out>/cjson-vectors.json`\n\n"
        "All inputs are synthetic. The file is owned by the catalog (CPython is the oracle); copy it, do not edit it.\n"
    ).encode("ascii")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cjson_vectors.py", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("write", help="write cjson-vectors.json and print its sha256")
    w.add_argument("out")
    w.add_argument("--source-md", action="store_true", help="also write SOURCE.md beside it")
    c = sub.add_parser("check", help="regenerate and compare with an existing file")
    c.add_argument("path")
    args = ap.parse_args(argv)
    data = render()
    sha = hashlib.sha256(data).hexdigest()
    if args.cmd == "check":
        path = Path(args.path)
        same = path.is_file() and path.read_bytes() == data
        print(f"{'identical' if same else 'differs'}\t{sha}\t{path.name}")
        return 0 if same else 1
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    if args.source_md:
        (out.parent / "SOURCE.md").write_bytes(source_md(sha))
    print(f"{sha}\t{out.name}\t{len(data)} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
