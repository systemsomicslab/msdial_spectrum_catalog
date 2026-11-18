"""WS8 P3 rule-file format: the stdlib reference of the MS-FINDER diagnostic-rule set (stdlib only, Python 3.10).

This is the Python reference that the C# design (WS8-P3 C# design, section 2.5) names. It implements the format of
section 3 and the mass and hypothesis arithmetic of sections 4.1-4.7 that the load-time checks need:

- the strict JSON reader (3.1), canonical JSON and every sha256 of 3.2 (file, rule set, rule semantics);
- the storage form writer (3.3) and the version pin against ruleset-history.tsv;
- the closed schemas of the tables (3.4) and the rule record (3.5), with the refused keys;
- the formula grammar, integer pico-dalton masses (4.1), composition series (4.2), adduct hypotheses for both
  modes (4.3) and term member targets (4.5);
- every load-time check of 3.6, each reported as RuleLoadError(code, file, RFC 6901 pointer), sorted by
  (file, pointer, code).

It never evaluates a spectrum: the evaluator is the C# engine (a later Python port follows it).

Layout of a rule set directory (paths relative to RuleSet/, every *.json is hashed):

    ruleset.json                      {"format":"msfinder-diagnostic-ruleset","format_version":1,"rule_set_version":..}
    ruleset-history.tsv               rule_set_version<TAB>ruleset_sha256<TAB>yyyy-mm-dd (not hashed)
    tables/elements.json              exactly as design 3.4
    tables/adducts.json               {"format":"msfinder-diagnostic-adducts","format_version":1,"adducts":[rows]}
    tables/units.json                 {"format":"msfinder-diagnostic-units","format_version":1,"units":[rows]}
    tables/motifs.json                {"format":"msfinder-diagnostic-motifs","format_version":1,"motifs":[rows]}
    tables/motif_nesting.json         {"format":"msfinder-diagnostic-motif-nesting","format_version":1,"nesting":[rows]}
    scaffolds/<list_id>.json          {"list_id":..,"entries":[..]}
    rules/<rule_id>.json              the rule record

The design describes the four tables as arrays, but its reader refuses a top-level array (top_level_not_object)
and every *.json is hashed through that reader, so each table is wrapped in an object with a format header, as
elements.json is. A draft directory (--draft) may also hold the staging folders _motifs/*.json and _nesting/*.json
(top-level arrays of rows) and _requests/*.json ({"units":[rows],"adducts":[rows]}); they are merged in memory
after the table rows and are never hashed. _audit/ is ignored in a draft.

    python -B tools/ws8_p3/rule_format.py validate <ruleset_dir> [--draft] [--rules <glob>]
    python -B tools/ws8_p3/rule_format.py fmt <file>
    python -B tools/ws8_p3/rule_format.py sha <ruleset_dir> [--draft]
    python -B tools/ws8_p3/rule_format.py masses <ruleset_dir> --rule <rule_id> [--draft]
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import datetime
import decimal
import fnmatch
import hashlib
import itertools
import json
from pathlib import Path
import re
import sys
from typing import NamedTuple

ENGINE_SEMANTICS = "msfinder-diagnostic-rules/1"
MAX_JSON_DEPTH = 128
MAX_COMPOSITIONS = 100_000
MAX_LOGIC_DEPTH = 16
MAX_LOGIC_NODES = 256
PDA_PER_DA = 10 ** 12
MASS_LIMIT_PDA = 2 ** 53
# A member within 1.0 Da of the hypothesis precursor is never evidence (design 3.6 item 6, 4.4); the comparison
# carries the engine's 1e-10 Da slack, which is 100 pDa.
PRECURSOR_REGION_PDA = PDA_PER_DA + 100
NOTE_MAX = 300

HISTORY_FILE = "ruleset-history.tsv"
RULESET_FILE = "ruleset.json"
ELEMENTS_FILE = "tables/elements.json"
ADDUCTS_FILE = "tables/adducts.json"
UNITS_FILE = "tables/units.json"
MOTIFS_FILE = "tables/motifs.json"
NESTING_FILE = "tables/motif_nesting.json"
TABLE_WRAPPERS = {
    ADDUCTS_FILE: ("msfinder-diagnostic-adducts", "adducts"),
    UNITS_FILE: ("msfinder-diagnostic-units", "units"),
    MOTIFS_FILE: ("msfinder-diagnostic-motifs", "motifs"),
    NESTING_FILE: ("msfinder-diagnostic-motif-nesting", "nesting"),
}
STAGING_DIRS = ("_motifs", "_nesting", "_requests")
DRAFT_IGNORED_DIRS = ("_audit",)

ORIGINS = ("as_published", "tuned_vs20", "tuned_pooled")
ORIGIN_RANK = {o: i for i, o in enumerate(ORIGINS)}
MOTIF_ORIGINS = ("source_defined", "tuned_vs20", "tuned_pooled")
UNIT_CLASSES = ("sugar", "acyl", "conjugate")
REQUIRED_UNITS = ("Hex", "dHex", "Pen", "HexA")
VARIABLES = ("n_units", "n_sugars", "n_hex", "n_dhex", "n_pen", "n_hexa", "n_acyl")
OUT_OF_SCOPE_VARIABLES = ("n_chains",)
ROLES = ("required", "counted", "ordering", "veto_present", "veto_stronger_than", "upgrade")
SUPPORT_ROLES = ("required", "counted", "ordering")
VETO_ROLES = ("veto_present", "veto_stronger_than")
KINDS = ("product_ion", "loss_from_precursor", "aglycone_ion", "ion_series", "window")
OUT_OF_SCOPE_KINDS = ("chain", "chains", "chain_set")
CHOICE_FIELDS = {"product_ion": ("formula", "unit", "scaffold_ions"), "loss_from_precursor": ("formula", "unit")}
FIXED_FIELDS = {
    "aglycone_ion": ("form", "charge_carrier", "plus", "minus"),
    "ion_series": ("unit_classes", "charge_carrier", "form", "plus", "minus"),
    "window": ("formula", "as", "half_width_da"),
}
ALL_KIND_FIELDS = ("formula", "unit", "scaffold_ions", "form", "charge_carrier", "plus", "minus", "unit_classes",
                   "as", "half_width_da")
TERM_KEYS = ("term_id", "kind", "role", "term_source", "min_rel_intensity", "floor_source", "than", "ratio",
             "ratio_source", "note")
RULE_KEYS = ("format", "format_version", "rule_id", "version", "source", "transcribed_by", "source_checked",
             "author_reviewed", "measure", "review", "parameters_origin", "primary_precondition_mode",
             "calibrated_on", "convention", "product_tolerance", "loss_reference", "veto_excludes_motif", "claim",
             "preconditions", "terms", "complementary", "logic", "notes")
REVIEW_KEYS = ("by", "at", "basis", "run_id", "rules_precision_sha256", "review_proposal_sha256",
               "rule_semantic_sha256", "instrument_classes", "product_tolerance", "precursor_tolerance_da")
INSTRUMENT_CLASSES = ("QTOF-type", "Orbitrap/FT", "other", "unknown")
CALIBRATION_CLASSES = ("QTOF-type", "Orbitrap/FT", "other")
REL_SCALE = {"p2": 1.0, "msdial": 100.0}

FORMULA_DEPENDENT_KEYS = ("element_gate", "oracle_formula", "predicted_formula", "precursor_formula",
                          "assigned_formula", "library_formula", "inchikey")
TYPED_MASS_KEYS = ("mz", "mass", "exact_mass", "monoisotopic_mass", "printed_mz")
OUT_OF_SCOPE_KEYS = ("chain", "chains", "chain_set", "max_rel_intensity")
REFUSED_KEYS = {**{k: "formula_dependent_key" for k in FORMULA_DEPENDENT_KEYS},
                **{k: "typed_mass" for k in TYPED_MASS_KEYS},
                **{k: "out_of_scope_kind" for k in OUT_OF_SCOPE_KEYS}}

# Every code this module reports. The first block is the strict reader (design 3.1); the codes the design names
# are kept verbatim, and the ones it leaves open are named here so the C# loader can mirror them.
ERROR_CODES = {
    "bom": "UTF-8 byte-order mark",
    "invalid_utf8": "bytes are not UTF-8",
    "syntax": "not RFC 8259 JSON (comments, trailing commas, raw control characters ...)",
    "top_level_not_object": "the top level is not an object",
    "duplicate_key": "a key occurs twice in one object",
    "nonfinite": "NaN or Infinity",
    "nesting_depth": "arrays and objects nested deeper than 128",
    "number_range": "an integer above 2^53-1 or a float outside 0 and [1e-4, 1e15)",
    "negative_zero": "-0, -0.0",
    "number_form": "a float with an exponent",
    "float_precision": "a float with more than 15 significant digits",
    "control_char": "a string with a code point below 0x20 or 0x7F",
    "non_ascii": "a string with a code point above 0x7F",
    "missing_file": "a required file is absent",
    "unexpected_file": "a *.json outside the layout (or a staging file without --draft)",
    "unknown_key": "a key outside the closed key set",
    "missing_key": "a required key is absent",
    "bad_type": "a value of the wrong JSON type (floats must be written with a decimal point)",
    "bad_value": "a value outside its enum, pattern, length or range",
    "kind_field_choice": "not exactly one of the alternative fields (term kind source, unit class/name)",
    "formula_dependent_key": "a key that reads a formula (directive: MS2 rules never consult a formula)",
    "typed_mass": "a typed mass; rule files hold formulas only",
    "out_of_scope_kind": "a lipid chain construct",
    "out_of_scope_variable": "n_chains",
    "bad_formula": "a formula outside the grammar",
    "formula_not_hill": "a formula not in strict Hill order, or an element repeated",
    "unknown_element": "an element missing from elements.json",
    "mass_precision": "an element or electron mass with more than 12 fraction digits",
    "missing_source": "an as_published (or source_defined) row with no source",
    "duplicate_id": "a table id (adduct name or alias, unit, motif, scaffold, nesting pair) defined twice",
    "adduct_polarity": "an adduct whose charge contradicts its name or the rule's ion mode",
    "adduct_name_mismatch": "an adduct name that disagrees with multiplier, add and remove",
    "missing_required_unit": "units.json lacks Hex, dHex, Pen or HexA",
    "nesting_cycle": "motif nesting with a cycle or a self-loop",
    "file_name_mismatch": "rule_id or list_id differs from the file stem",
    "mode_mismatch": "the rule_id mode differs from preconditions.ion_mode",
    "unknown_term": "a term id that no term defines",
    "unknown_unit": "a unit name missing from units.json",
    "unknown_scaffold_list": "a scaffold list that does not exist, or scaffold_ions without a scaffold list",
    "unknown_ion_key": "a scaffold ion key missing (or empty for the mode) in some entry",
    "unknown_adduct": "an adduct name missing from adducts.json",
    "unknown_motif": "a motif id missing from motifs.json",
    "unknown_source_ref": "a provenance token src:<i> past the end of source",
    "unknown_variable": "a k variable outside the closed set",
    "duplicate_term": "a term id defined twice",
    "unused_term": "a non-upgrade term that the logic never uses",
    "claim_term_role": "a claim term that is repeated or is not a support term",
    "claim_terms_too_few": "a measured rule claiming on fewer than two terms (design 3.5)",
    "role_placement": "a role outside its place in the logic (design 3.6 item 3)",
    "logic_node": "a logic node that is not exactly one of all, any, not, k_of, term",
    "logic_depth": "logic deeper than 16 or larger than 256 nodes",
    "k_range": "k or its min below 1",
    "series_too_large": "more than 100 000 compositions",
    "series_empty": "no composition satisfies the unit ranges and totals",
    "negative_count": "a negative count (unit range, ion or precursor composition)",
    "mass_range": "a mass that is not positive or not below 2^53 pDa",
    "complementary_pair": "a complementary pair of one term, of one kind, of non-support terms, or repeated",
    "complementary_multi_member": "a complementary term with more than one member in some hypothesis",
    "single_ion_claim": "the logic holds with no evidence unit or with one (two-ion check, pass 1)",
    "claim_single_ion": "with every non-claim unit present, no or one claim unit satisfies the logic (pass 2)",
    "support_target_is_precursor": "a support or upgrade member equals a hypothesis precursor m/z",
    "origin_inconsistent": "parameters_origin differs from the derived maximum",
    "review_missing": "author_reviewed and review disagree",
    "review_not_allowed": "author_reviewed without source_checked and measure",
    "review_stale": "review.rule_semantic_sha256 differs from the recomputed one",
    "history_format": "a malformed or repeated ruleset-history.tsv row",
    "ruleset_unpinned": "rule_set_version has no ruleset-history.tsv row",
    "ruleset_sha_mismatch": "the history row for rule_set_version names another sha256",
}


# ------------------------------------------------------------------- strict reader and hashes (design 3.1, 3.2)
_FLT = re.compile(r"-?(?:0|[1-9][0-9]*)\.[0-9]+\Z")


class RuleFormatError(ValueError):
    def __init__(self, code):
        super().__init__(code)
        self.code = code


def _pi(lex):
    if lex == "-0":
        raise RuleFormatError("negative_zero")
    if len(lex.lstrip("-")) > 16:  # before int(): no 4300-digit ValueError
        raise RuleFormatError("number_range")
    v = int(lex)
    if abs(v) > 2 ** 53 - 1:
        raise RuleFormatError("number_range")
    return v


def _pf(lex):
    if not _FLT.match(lex):
        raise RuleFormatError("number_form")
    if len(lex.lstrip("-").replace(".", "").lstrip("0").rstrip("0")) > 15:
        raise RuleFormatError("float_precision")
    v = float(lex)
    if v == 0.0 and lex.startswith("-"):
        raise RuleFormatError("negative_zero")
    if v != 0.0 and not (1e-4 <= abs(v) < 1e15):
        raise RuleFormatError("number_range")
    return v


def _pairs(ps):
    d = {}
    for k, v in ps:
        if k in d:
            raise RuleFormatError("duplicate_key")
        d[k] = v
    return d


def _nonfinite(c):
    raise RuleFormatError("nonfinite")


def _check_printable_ascii(o):  # depth first, document order, each key before its value
    if isinstance(o, str):
        for ch in o:  # the first offending code point decides
            c = ord(ch)
            if c < 0x20 or c == 0x7F:
                raise RuleFormatError("control_char")
            if c > 0x7F:
                raise RuleFormatError("non_ascii")
    elif isinstance(o, dict):
        for k, v in o.items():
            _check_printable_ascii(k)
            _check_printable_ascii(v)
    elif isinstance(o, list):
        for v in o:
            _check_printable_ascii(v)


def _depth(o):  # iterative, so it cannot itself recurse too deep
    best, stack = 0, [(o, 1)]
    while stack:
        v, d = stack.pop()
        if isinstance(v, (dict, list)):
            best = max(best, d)
            stack.extend((c, d + 1) for c in (v.values() if isinstance(v, dict) else v))
    return best


def _load(b, allow_array=False):
    if b.startswith(b"\xef\xbb\xbf"):
        raise RuleFormatError("bom")
    try:
        s = b.decode("utf-8")
    except UnicodeDecodeError:
        raise RuleFormatError("invalid_utf8") from None
    try:
        o = json.loads(s, object_pairs_hook=_pairs, parse_int=_pi, parse_float=_pf, parse_constant=_nonfinite)
    except json.JSONDecodeError:
        raise RuleFormatError("syntax") from None
    except RecursionError:
        raise RuleFormatError("nesting_depth") from None
    if not (isinstance(o, dict) or (allow_array and isinstance(o, list))):
        raise RuleFormatError("top_level_not_object")
    if _depth(o) > MAX_JSON_DEPTH:
        raise RuleFormatError("nesting_depth")
    _check_printable_ascii(o)
    return o


def load_strict(b: bytes) -> dict:
    """Design 3.2 load_strict: the only reader of hashed files."""
    return _load(b)


def load_strict_staging(b: bytes):
    """The same reader, also accepting a top-level array (draft staging files only; never hashed)."""
    return _load(b, allow_array=True)


def canonical_bytes(o) -> bytes:
    return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("ascii")


def file_sha256(b: bytes) -> str:
    return hashlib.sha256(canonical_bytes(load_strict(b))).hexdigest()


def ruleset_sha256(files: dict) -> str:
    """{posix path relative to RuleSet/: bytes}, every hashed *.json."""
    m = {"files": {p: file_sha256(b) for p, b in files.items()},
         "format": "msfinder-diagnostic-ruleset-hash", "format_version": 1}
    return hashlib.sha256(canonical_bytes(m)).hexdigest()


def rule_semantic_sha256(rule: dict, resolved: dict, engine_semantics: str = ENGINE_SEMANTICS) -> str:
    r = {k: v for k, v in rule.items() if k not in ("author_reviewed", "review")}
    m = {"engine_semantics": engine_semantics, "format": "msfinder-diagnostic-rule-semantics",
         "format_version": 1, "resolved": resolved, "rule": r}
    return hashlib.sha256(canonical_bytes(m)).hexdigest()


def storage_bytes(o) -> bytes:
    """Storage form (design 3.3): sorted keys, indent 1, ASCII, LF, one trailing newline."""
    return (json.dumps(o, sort_keys=True, indent=1, ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def write_storage(path, o) -> bool:
    """Write o in storage form after checking that the strict reader accepts the result. True if changed."""
    data = storage_bytes(o)
    _load(data, allow_array=True)
    path = Path(path)
    if path.exists() and path.read_bytes() == data:
        return False
    path.write_bytes(data)
    return True


# --------------------------------------------------------------------------------------------- errors, pointers
class RuleLoadError(NamedTuple):
    code: str
    file: str
    pointer: str

    def __str__(self):
        return f"{self.code}\t{self.file}\t{self.pointer}"


def sort_errors(errors):
    return sorted(set(errors), key=lambda e: (e.file, e.pointer, e.code))


def jp(base, *tokens):
    """Extend an RFC 6901 pointer."""
    for t in tokens:
        base += "/" + str(t).replace("~", "~0").replace("/", "~1")
    return base


class Ctx:
    """Error sink for one file."""

    def __init__(self, sink, file):
        self.sink = sink
        self.file = file

    def err(self, code, pointer):
        assert code in ERROR_CODES, code
        self.sink.append(RuleLoadError(code, self.file, pointer))

    @property
    def n(self):
        return len(self.sink)


# -------------------------------------------------------------------------------------------- formulas, masses
_FORMULA_RE = re.compile(r"(?:[A-Z][a-z]?(?:[2-9]|[1-9][0-9]{1,2})?)+\Z")
_TOKEN_RE = re.compile(r"([A-Z][a-z]?)([0-9]*)")
_LOOSE_PART_RE = re.compile(r"(?:[A-Z][a-z]?(?:[1-9][0-9]*)?)+\Z")


def hill_order(symbols):
    s = set(symbols)
    if "C" in s:
        return ["C"] + (["H"] if "H" in s else []) + sorted(s - {"C", "H"})
    return sorted(s)


def parse_formula(text, elements) -> dict:
    """Strict Hill formula -> {element: count}. Raises RuleFormatError(bad_formula|unknown_element|formula_not_hill)."""
    if not isinstance(text, str) or not _FORMULA_RE.match(text):
        raise RuleFormatError("bad_formula")
    tokens = _TOKEN_RE.findall(text)
    symbols = [s for s, _ in tokens]
    for s in symbols:
        if s not in elements:
            raise RuleFormatError("unknown_element")
    if len(set(symbols)) != len(symbols) or symbols != hill_order(symbols):
        raise RuleFormatError("formula_not_hill")
    return {s: int(n) if n else 1 for s, n in tokens}


def hill_format(comp) -> str:
    items = {s: n for s, n in comp.items() if n}
    if any(n < 0 for n in items.values()):
        raise ValueError("negative count")
    return "".join(s + (str(items[s]) if items[s] != 1 else "") for s in hill_order(items))


def comp_add(a, b, k=1):
    out = dict(a)
    for s, n in b.items():
        out[s] = out.get(s, 0) + k * n
    return out


def comp_negative(comp):
    return any(n < 0 for n in comp.values())


def comp_empty(comp):
    return not any(comp.values())


def lexeme_pda(x) -> int:
    """A parsed float -> exact pico-daltons. repr() has the stripped lexeme's digits under the number subset."""
    q = decimal.Decimal(repr(x)).scaleb(12)
    if q != q.to_integral_value():
        raise RuleFormatError("mass_precision")
    return int(q)


def mz_of(pda) -> float:
    return pda / 1e12


@dataclass
class Elements:
    electron: float
    electron_pda: int
    masses: dict
    pda: dict

    def mass_pda(self, comp):
        return sum(n * self.pda[s] for s, n in comp.items())

    def ion_pda(self, comp, charge):
        return self.mass_pda(comp) - charge * self.electron_pda


# ------------------------------------------------------------------------------------------- generic validators
_MISSING = object()


def is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def check_keys(c, obj, ptr, required, optional=()):
    """Closed key set: report unknown and missing keys (refused keys are reported by scan_refused)."""
    for k in obj:
        if k in REFUSED_KEYS:
            continue
        if k not in required and k not in optional:
            c.err("unknown_key", jp(ptr, k))
    for k in required:
        if k not in obj:
            c.err("missing_key", jp(ptr, k))


def scan_refused(c, o, ptr=""):
    stack = [(o, ptr)]
    while stack:
        v, p = stack.pop()
        if isinstance(v, dict):
            for k, x in v.items():
                if k in REFUSED_KEYS:
                    c.err(REFUSED_KEYS[k], jp(p, k))
                stack.append((x, jp(p, k)))
        elif isinstance(v, list):
            stack.extend((x, jp(p, i)) for i, x in enumerate(v))


def v_obj(c, v, ptr):
    if not isinstance(v, dict):
        c.err("bad_type", ptr)
        return False
    return True


def v_list(c, v, ptr, min_len=0):
    if not isinstance(v, list):
        c.err("bad_type", ptr)
        return False
    if len(v) < min_len:
        c.err("bad_value", ptr)
        return False
    return True


def v_str(c, v, ptr, pattern=None, maxlen=None, nullable=False, nonempty=True):
    if v is _MISSING:
        return False
    if v is None and nullable:
        return True
    if not isinstance(v, str):
        c.err("bad_type", ptr)
        return False
    if (nonempty and v == "") or (maxlen is not None and len(v) > maxlen) or (pattern and not pattern.match(v)):
        c.err("bad_value", ptr)
        return False
    return True


def v_enum(c, v, ptr, allowed):
    if v is _MISSING:
        return False
    if not isinstance(v, str):
        c.err("bad_type", ptr)
        return False
    if v not in allowed:
        c.err("bad_value", ptr)
        return False
    return True


def v_bool(c, v, ptr):
    if v is _MISSING:
        return False
    if not isinstance(v, bool):
        c.err("bad_type", ptr)
        return False
    return True


def v_int(c, v, ptr, minimum=None, code="bad_value"):
    if v is _MISSING:
        return False
    if not is_int(v):
        c.err("bad_type", ptr)
        return False
    if minimum is not None and v < minimum:
        c.err(code, ptr)
        return False
    return True


def v_float(c, v, ptr, positive=False, maximum=None, nullable=False):
    if v is _MISSING:
        return False
    if v is None and nullable:
        return True
    if not isinstance(v, float):
        c.err("bad_type", ptr)
        return False
    if (positive and not v > 0.0) or v < 0.0 or (maximum is not None and v > maximum):
        c.err("bad_value", ptr)
        return False
    return True


def v_formula(c, v, ptr, elements, nullable=False):
    """-> composition dict, None for an allowed null, or _MISSING after an error."""
    if v is _MISSING:
        return _MISSING
    if v is None and nullable:
        return None
    if not isinstance(v, str):
        c.err("bad_type", ptr)
        return _MISSING
    try:
        return parse_formula(v, elements)
    except RuleFormatError as e:
        c.err(e.code, ptr)
        return _MISSING


def v_format_header(c, obj, fmt):
    if obj.get("format", fmt) != fmt:
        c.err("bad_value", "/format")
    v = obj.get("format_version", 1)
    if not is_int(v):
        c.err("bad_type", "/format_version")
    elif v != 1:
        c.err("bad_value", "/format_version")


_DOI_RE = re.compile(r"10\.[0-9]{4,9}/[!-~]+\Z")
_SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
_ID_RE = re.compile(r"[a-z0-9_]+\Z")
_TERM_ID_RE = re.compile(r"[A-Za-z0-9_]{1,64}\Z")
_UNIT_RE = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,31}\Z")
_RULE_ID_RE = re.compile(r"[a-z0-9_]+\.[a-z0-9_]+\.(pos|neg)\Z")
_VERSION_RE = re.compile(r"[a-z0-9][a-z0-9._-]*\Z")
_CHEMONT_RE = re.compile(r"CHEMONTID:[0-9]{7}\Z")
_SYMBOL_RE = re.compile(r"[A-Z][a-z]?\Z")
_TOKEN_SRC_RE = re.compile(r"src:(0|[1-9][0-9]{0,5})\Z")
_DATE_RE = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}\Z")
_ADDUCT_RE = re.compile(r"\[([1-9][0-9]*)?M((?:[+-](?:[1-9][0-9]*)?[A-Z][A-Za-z0-9]*)*)\]([1-9][0-9]*)?([+-])\Z")
_ADDUCT_PART_RE = re.compile(r"([+-])([1-9][0-9]*)?([A-Z][A-Za-z0-9]*)")


def v_citation(c, v, ptr):
    if not v_obj(c, v, ptr):
        return
    if "doi" in v:
        check_keys(c, v, ptr, ("doi", "locator"))
        v_str(c, v["doi"], jp(ptr, "doi"), pattern=_DOI_RE)
    elif "reference" in v:
        check_keys(c, v, ptr, ("reference", "locator"))
        v_str(c, v["reference"], jp(ptr, "reference"), maxlen=NOTE_MAX)
    else:
        c.err("missing_key", jp(ptr, "doi"))
        check_keys(c, v, ptr, ("locator",))
        return
    v_str(c, v.get("locator", _MISSING), jp(ptr, "locator"), maxlen=NOTE_MAX)


def v_citations(c, v, ptr, min_len=0):
    if not v_list(c, v, ptr, min_len):
        return
    for i, x in enumerate(v):
        v_citation(c, x, jp(ptr, i))


def v_token(c, v, ptr, n_source, src_only=False):
    """Provenance token -> origin, or None after an error."""
    if v is _MISSING:
        return None
    if not isinstance(v, str):
        c.err("bad_type", ptr)
        return None
    if v in ("tuned_vs20", "tuned_pooled") and not src_only:
        return v
    m = _TOKEN_SRC_RE.match(v)
    if not m:
        c.err("bad_value", ptr)
        return None
    if n_source is not None and int(m.group(1)) >= n_source:
        c.err("unknown_source_ref", ptr)
        return None
    return "as_published"


def v_range_pair(c, v, ptr):
    """[min, max] of non-negative ints -> tuple, or None."""
    if not isinstance(v, list) or len(v) != 2 or not all(is_int(x) for x in v):
        c.err("bad_type", ptr)
        return None
    if v[0] < 0 or v[1] < 0:
        c.err("negative_count", ptr)
        return None
    if v[0] > v[1]:
        c.err("bad_value", ptr)
        return None
    return v[0], v[1]


def v_row_origin(c, row, ptr, origins):
    ok = v_enum(c, row.get("origin", _MISSING), jp(ptr, "origin"), origins)
    src = row.get("source", _MISSING)
    if src is not _MISSING:
        v_citations(c, src, jp(ptr, "source"))
    if ok and row["origin"] == origins[0] and isinstance(src, list) and not src:
        c.err("missing_source", jp(ptr, "source"))


# ------------------------------------------------------------------------------------------------------ tables
@dataclass
class Adduct:
    row: dict
    name: str
    charge: int
    multiplier: int
    add: dict
    remove: dict
    common: bool
    origin: str


@dataclass
class Unit:
    row: dict
    name: str
    cls: str
    residue: dict
    origin: str


@dataclass
class ScaffoldList:
    doc: dict
    list_id: str
    entries: list  # (scaffold_id, formula comp, ions {key: {"pos": [comp], "neg": [comp]}}, origin, entry dict)


@dataclass
class Tables:
    elements: Elements
    adducts: list = field(default_factory=list)
    adduct_by_name: dict = field(default_factory=dict)
    units: list = field(default_factory=list)
    unit_by_name: dict = field(default_factory=dict)
    motifs: dict = field(default_factory=dict)
    nesting: list = field(default_factory=list)
    scaffolds: dict = field(default_factory=dict)


def load_elements(c, doc):
    check_keys(c, doc, "", ("format", "format_version", "electron", "elements"))
    n0 = c.n
    v_format_header(c, doc, "msfinder-diagnostic-elements")
    electron = doc.get("electron", _MISSING)
    e_pda = None
    if v_float(c, electron, "/electron", positive=True):
        try:
            e_pda = lexeme_pda(electron)
        except RuleFormatError as e:
            c.err(e.code, "/electron")
    masses, pda = {}, {}
    els = doc.get("elements", _MISSING)
    if els is not _MISSING and v_obj(c, els, "/elements"):
        if not els:
            c.err("bad_value", "/elements")
        for sym, m in els.items():
            p = jp("/elements", sym)
            if not _SYMBOL_RE.match(sym):
                c.err("bad_value", p)
                continue
            if v_float(c, m, p, positive=True):
                try:
                    pda[sym] = lexeme_pda(m)
                    masses[sym] = m
                except RuleFormatError as e:
                    c.err(e.code, p)
    if c.n != n0 or e_pda is None:
        return None
    return Elements(electron, e_pda, masses, pda)


def adduct_name_net(name, elements):
    """Canonical adduct name -> (multiplier, charge count, sign, net composition) or None if not formula-based."""
    m = _ADDUCT_RE.match(name)
    if not m:
        return None
    net = {}
    for sign, count, part in _ADDUCT_PART_RE.findall(m.group(2)):
        if not _LOOSE_PART_RE.match(part):
            return None
        for sym, n in _TOKEN_RE.findall(part):
            if sym not in elements:
                return None
            net[sym] = net.get(sym, 0) + (1 if sign == "+" else -1) * int(count or 1) * int(n or 1)
    return int(m.group(1) or 1), int(m.group(3) or 1), m.group(4), {s: n for s, n in net.items() if n}


def v_adduct_row(c, row, ptr, el):
    if not v_obj(c, row, ptr):
        return None
    n0 = c.n
    check_keys(c, row, ptr, ("adduct", "aliases", "charge", "multiplier", "add", "remove", "common", "origin",
                             "source"))
    name = row.get("adduct", _MISSING)
    name_ok = v_str(c, name, jp(ptr, "adduct"), pattern=_ADDUCT_RE)
    aliases = row.get("aliases", _MISSING)
    if aliases is not _MISSING and v_list(c, aliases, jp(ptr, "aliases")):
        for i, a in enumerate(aliases):
            v_str(c, a, jp(ptr, "aliases", i), maxlen=64)
    charge = row.get("charge", _MISSING)
    if v_int(c, charge, jp(ptr, "charge")) and charge not in (1, -1):
        c.err("bad_value", jp(ptr, "charge"))
    v_int(c, row.get("multiplier", _MISSING), jp(ptr, "multiplier"), minimum=1)
    add = v_formula(c, row.get("add", _MISSING), jp(ptr, "add"), el.pda, nullable=True)
    remove = v_formula(c, row.get("remove", _MISSING), jp(ptr, "remove"), el.pda, nullable=True)
    v_bool(c, row.get("common", _MISSING), jp(ptr, "common"))
    v_row_origin(c, row, ptr, ORIGINS)
    if c.n != n0:
        return None
    add, remove = add or {}, remove or {}
    parsed = adduct_name_net(name, el.pda) if name_ok else None
    if parsed is None:
        c.err("adduct_name_mismatch", jp(ptr, "adduct"))
        return None
    mult, zcount, sign, net = parsed
    if (sign == "+") != (row["charge"] == 1) or zcount != 1:
        c.err("adduct_polarity", jp(ptr, "charge"))
        return None
    own_net = {s: n for s, n in comp_add(add, remove, -1).items() if n}
    if mult != row["multiplier"] or net != own_net:
        c.err("adduct_name_mismatch", jp(ptr, "adduct"))
        return None
    return Adduct(row, name, row["charge"], row["multiplier"], add, remove, row["common"], row["origin"])


def v_unit_row(c, row, ptr, el):
    if not v_obj(c, row, ptr):
        return None
    n0 = c.n
    check_keys(c, row, ptr, ("name", "class", "residue", "origin", "source"))
    v_str(c, row.get("name", _MISSING), jp(ptr, "name"), pattern=_UNIT_RE)
    v_enum(c, row.get("class", _MISSING), jp(ptr, "class"), UNIT_CLASSES)
    residue = v_formula(c, row.get("residue", _MISSING), jp(ptr, "residue"), el.pda)
    v_row_origin(c, row, ptr, ORIGINS)
    if c.n != n0:
        return None
    return Unit(row, row["name"], row["class"], residue, row["origin"])


def v_motif_row(c, row, ptr):
    if not v_obj(c, row, ptr):
        return None
    n0 = c.n
    check_keys(c, row, ptr, ("motif_id", "smarts", "origin", "source", "chemont"))
    v_str(c, row.get("motif_id", _MISSING), jp(ptr, "motif_id"), pattern=_ID_RE)
    v_str(c, row.get("smarts", _MISSING), jp(ptr, "smarts"), maxlen=2000)
    v_row_origin(c, row, ptr, MOTIF_ORIGINS)
    chem = row.get("chemont", _MISSING)
    if chem is not _MISSING and v_obj(c, chem, jp(ptr, "chemont")):
        check_keys(c, chem, jp(ptr, "chemont"), ("compatible", "broader"))
        for key in ("compatible", "broader"):
            nodes = chem.get(key, _MISSING)
            if nodes is _MISSING or not v_list(c, nodes, jp(ptr, "chemont", key)):
                continue
            for i, node in enumerate(nodes):
                p = jp(ptr, "chemont", key, i)
                if v_obj(c, node, p):
                    check_keys(c, node, p, ("id", "name"))
                    v_str(c, node.get("id", _MISSING), jp(p, "id"), pattern=_CHEMONT_RE)
                    v_str(c, node.get("name", _MISSING), jp(p, "name"), maxlen=NOTE_MAX)
    return row if c.n == n0 else None


def v_nesting_row(c, row, ptr):
    if not v_obj(c, row, ptr):
        return None
    n0 = c.n
    check_keys(c, row, ptr, ("child", "parent", "source"))
    v_str(c, row.get("child", _MISSING), jp(ptr, "child"), pattern=_ID_RE)
    v_str(c, row.get("parent", _MISSING), jp(ptr, "parent"), pattern=_ID_RE)
    if "source" in row:
        v_citations(c, row["source"], jp(ptr, "source"), min_len=1)
    return row if c.n == n0 else None


def v_scaffold_list(c, doc, stem, el):
    n0 = c.n
    check_keys(c, doc, "", ("list_id", "entries"))
    list_id = doc.get("list_id", _MISSING)
    if v_str(c, list_id, "/list_id", pattern=_ID_RE) and list_id != stem:
        c.err("file_name_mismatch", "/list_id")
    entries = doc.get("entries", _MISSING)
    out, seen = [], set()
    if entries is not _MISSING and v_list(c, entries, "/entries", min_len=1):
        for i, e in enumerate(entries):
            p = jp("/entries", i)
            if not v_obj(c, e, p):
                continue
            check_keys(c, e, p, ("scaffold_id", "name", "formula", "class", "ions", "origin", "source"))
            sid = e.get("scaffold_id", _MISSING)
            if v_str(c, sid, jp(p, "scaffold_id"), pattern=_ID_RE):
                if sid in seen:
                    c.err("duplicate_id", jp(p, "scaffold_id"))
                seen.add(sid)
            v_str(c, e.get("name", _MISSING), jp(p, "name"), maxlen=NOTE_MAX)
            comp = v_formula(c, e.get("formula", _MISSING), jp(p, "formula"), el.pda)
            v_str(c, e.get("class", _MISSING), jp(p, "class"), pattern=_ID_RE)
            ions_out = {}
            ions = e.get("ions", _MISSING)
            if ions is not _MISSING and v_obj(c, ions, jp(p, "ions")):
                for key, modes in ions.items():
                    pk = jp(p, "ions", key)
                    if not _ID_RE.match(key):
                        c.err("bad_value", pk)
                        continue
                    if not v_obj(c, modes, pk):
                        continue
                    check_keys(c, modes, pk, ("pos", "neg"))
                    ions_out[key] = {}
                    for mode in ("pos", "neg"):
                        lst = modes.get(mode, _MISSING)
                        if lst is _MISSING or not v_list(c, lst, jp(pk, mode)):
                            continue
                        ions_out[key][mode] = [v_formula(c, f, jp(pk, mode, j), el.pda) for j, f in enumerate(lst)]
            v_row_origin(c, e, p, ORIGINS)
            out.append((sid, comp, ions_out, e.get("origin"), e))
    if c.n != n0:
        return None
    return ScaffoldList(doc, list_id, out)


# ---------------------------------------------------------------------------------------------- rule checking
@dataclass
class Term:
    index: int
    tid: str
    kind: str
    role: str
    obj: dict
    formula: dict | None = None
    unit_sel: tuple | None = None  # ("class"|"name", value)
    plus: dict | None = None
    minus: dict | None = None
    scaffold_key: str | None = None


@dataclass
class Hypothesis:
    entry_id: str
    entry_index: int
    counts: tuple
    adduct: Adduct
    adduct_source: str
    agly: dict
    comp: dict
    neutral_pda: int
    precursor_pda: int
    composition_id: str
    hypothesis_id: str


@dataclass
class RuleInfo:
    """A rule whose schema and references passed: enough to enumerate hypotheses and targets."""
    rule_id: str
    doc: dict
    file: str
    mode: str
    charge: int
    terms: list
    term_by_id: dict
    scaffold: ScaffoldList | None
    core: dict
    unit_ranges: list  # (lo, hi) per units-table row
    totals: dict
    tried: list  # [(Adduct, "rule"|"common")]
    logic: tuple
    k_nodes: list  # [(pointer, kspec)]
    claim_terms: tuple
    groups: list  # list of sets of term indices
    pair_index: dict  # term index -> [pair indices]


def walk_logic_schema(c, node, ptr, n_source, stats, depth):
    stats["nodes"] += 1
    stats["depth"] = max(stats["depth"], depth)
    if stats["nodes"] > MAX_LOGIC_NODES or depth > MAX_LOGIC_DEPTH:
        return
    if not isinstance(node, dict) or len(node) != 1:
        c.err("logic_node", ptr)
        return
    (key, val), = node.items()
    p = jp(ptr, key)
    if key == "term":
        if not isinstance(val, str):
            c.err("bad_type", p)
    elif key in ("all", "any"):
        if not isinstance(val, list) or not val:
            c.err("logic_node", p)
            return
        for i, child in enumerate(val):
            walk_logic_schema(c, child, jp(p, i), n_source, stats, depth + 1)
    elif key == "not":
        walk_logic_schema(c, val, p, n_source, stats, depth + 1)
    elif key == "k_of":
        if not isinstance(val, dict):
            c.err("logic_node", p)
            return
        check_keys(c, val, p, ("k", "k_source", "of"))
        k = val.get("k", _MISSING)
        if is_int(k):
            if k < 1:
                c.err("k_range", jp(p, "k"))
        elif isinstance(k, dict):
            pk = jp(p, "k")
            check_keys(c, k, pk, ("var", "offset", "min"))
            var = k.get("var", _MISSING)
            if var is not _MISSING:
                if not isinstance(var, str):
                    c.err("bad_type", jp(pk, "var"))
                elif var in OUT_OF_SCOPE_VARIABLES:
                    c.err("out_of_scope_variable", jp(pk, "var"))
                elif var not in VARIABLES:
                    c.err("unknown_variable", jp(pk, "var"))
            v_int(c, k.get("offset", _MISSING), jp(pk, "offset"))
            v_int(c, k.get("min", _MISSING), jp(pk, "min"), minimum=1, code="k_range")
        elif k is not _MISSING:
            c.err("bad_type", jp(p, "k"))
        v_token(c, val.get("k_source", _MISSING), jp(p, "k_source"), n_source)
        of = val.get("of", _MISSING)
        if of is _MISSING:
            return
        if not isinstance(of, list) or not of:
            c.err("logic_node", jp(p, "of"))
            return
        for i, child in enumerate(of):
            walk_logic_schema(c, child, jp(p, "of", i), n_source, stats, depth + 1)
    else:
        c.err("logic_node", ptr)


def v_term(c, t, ptr, n_source, convention, el):
    if not v_obj(c, t, ptr):
        return
    kind = t.get("kind", _MISSING)
    kind_ok = False
    if kind is not _MISSING:
        if not isinstance(kind, str):
            c.err("bad_type", jp(ptr, "kind"))
        elif kind in OUT_OF_SCOPE_KINDS:
            c.err("out_of_scope_kind", jp(ptr, "kind"))
        elif kind not in KINDS:
            c.err("bad_value", jp(ptr, "kind"))
        else:
            kind_ok = True
    if not kind_ok:
        check_keys(c, t, ptr, TERM_KEYS, ALL_KIND_FIELDS)
        return
    if kind in CHOICE_FIELDS:
        present = [k for k in CHOICE_FIELDS[kind] if k in t]
        if len(present) != 1:
            c.err("kind_field_choice", ptr)
        fields = tuple(present[:1])
        check_keys(c, t, ptr, TERM_KEYS + fields, tuple(present[1:]))
    else:
        fields = FIXED_FIELDS[kind]
        check_keys(c, t, ptr, TERM_KEYS + fields)
    v_str(c, t.get("term_id", _MISSING), jp(ptr, "term_id"), pattern=_TERM_ID_RE)
    role = t.get("role", _MISSING)
    if v_enum(c, role, jp(ptr, "role"), ROLES) and kind == "window" and role not in VETO_ROLES:
        c.err("role_placement", jp(ptr, "role"))
    v_token(c, t.get("term_source", _MISSING), jp(ptr, "term_source"), n_source)
    v_token(c, t.get("floor_source", _MISSING), jp(ptr, "floor_source"), n_source)
    v_float(c, t.get("min_rel_intensity", _MISSING), jp(ptr, "min_rel_intensity"),
            maximum=REL_SCALE.get(convention))
    than = t.get("than", _MISSING)
    if than is not _MISSING and than is not None and v_list(c, than, jp(ptr, "than"), min_len=1):
        for i, x in enumerate(than):
            v_str(c, x, jp(ptr, "than", i), pattern=_TERM_ID_RE)
    if than is not _MISSING and role in ROLES:
        if (role in ("ordering", "veto_stronger_than")) != (than is not None):
            c.err("bad_value", jp(ptr, "than"))
    ratio = t.get("ratio", _MISSING)
    if v_float(c, ratio, jp(ptr, "ratio"), positive=True, nullable=True) and ratio is not None and than is None:
        c.err("bad_value", jp(ptr, "ratio"))
    rs = t.get("ratio_source", _MISSING)
    if rs is not _MISSING and ratio is not _MISSING:
        if (rs is None) != (ratio is None):
            c.err("bad_value", jp(ptr, "ratio_source"))
        elif rs is not None:
            v_token(c, rs, jp(ptr, "ratio_source"), n_source)
    v_str(c, t.get("note", _MISSING), jp(ptr, "note"), maxlen=NOTE_MAX, nullable=True)
    if "formula" in fields:
        v_formula(c, t.get("formula", _MISSING), jp(ptr, "formula"), el.pda)
    if "unit" in fields:
        u = t.get("unit", _MISSING)
        pu = jp(ptr, "unit")
        if u is not _MISSING and v_obj(c, u, pu):
            sel = [k for k in ("class", "name") if k in u]
            if len(sel) != 1:
                c.err("kind_field_choice", pu)
            check_keys(c, u, pu, tuple(sel[:1]) + ("plus", "minus"), tuple(sel[1:]))
            if "class" in sel[:1]:
                v_enum(c, u["class"], jp(pu, "class"), UNIT_CLASSES)
            if "name" in sel[:1]:
                v_str(c, u["name"], jp(pu, "name"), pattern=_UNIT_RE)
            v_formula(c, u.get("plus", _MISSING), jp(pu, "plus"), el.pda, nullable=True)
            v_formula(c, u.get("minus", _MISSING), jp(pu, "minus"), el.pda, nullable=True)
    if "scaffold_ions" in fields:
        v_str(c, t.get("scaffold_ions", _MISSING), jp(ptr, "scaffold_ions"), pattern=_ID_RE)
    if "form" in fields:
        form_ok = v_enum(c, t.get("form", _MISSING), jp(ptr, "form"), ("even", "radical"))
        carrier_ok = v_enum(c, t.get("charge_carrier", _MISSING), jp(ptr, "charge_carrier"), ("proton", "adduct"))
        if form_ok and carrier_ok and t["form"] == "radical" and t["charge_carrier"] == "adduct":
            c.err("bad_value", jp(ptr, "charge_carrier"))
        v_formula(c, t.get("plus", _MISSING), jp(ptr, "plus"), el.pda, nullable=True)
        v_formula(c, t.get("minus", _MISSING), jp(ptr, "minus"), el.pda, nullable=True)
    if "unit_classes" in fields:
        uc = t.get("unit_classes", _MISSING)
        if uc is not _MISSING and v_list(c, uc, jp(ptr, "unit_classes"), min_len=1):
            for i, x in enumerate(uc):
                v_enum(c, x, jp(ptr, "unit_classes", i), UNIT_CLASSES)
            if len(set(map(str, uc))) != len(uc):
                c.err("bad_value", jp(ptr, "unit_classes"))
    if "as" in fields:
        v_enum(c, t.get("as", _MISSING), jp(ptr, "as"), ("ion", "loss"))
        v_float(c, t.get("half_width_da", _MISSING), jp(ptr, "half_width_da"), positive=True)


def v_review(c, r, ptr):
    if r is None or not v_obj(c, r, ptr):
        return
    check_keys(c, r, ptr, REVIEW_KEYS)
    v_str(c, r.get("by", _MISSING), jp(ptr, "by"), maxlen=NOTE_MAX)
    at = r.get("at", _MISSING)
    if v_str(c, at, jp(ptr, "at"), pattern=_DATE_RE):
        try:
            datetime.date.fromisoformat(at)
        except ValueError:
            c.err("bad_value", jp(ptr, "at"))
    v_enum(c, r.get("basis", _MISSING), jp(ptr, "basis"), ("benchmark",))
    v_str(c, r.get("run_id", _MISSING), jp(ptr, "run_id"), maxlen=NOTE_MAX)
    for k in ("rules_precision_sha256", "review_proposal_sha256", "rule_semantic_sha256"):
        v_str(c, r.get(k, _MISSING), jp(ptr, k), pattern=_SHA_RE)
    ic = r.get("instrument_classes", _MISSING)
    if ic is not _MISSING and v_list(c, ic, jp(ptr, "instrument_classes"), min_len=1):
        for i, x in enumerate(ic):
            v_enum(c, x, jp(ptr, "instrument_classes", i), INSTRUMENT_CLASSES)
        if len(set(map(str, ic))) != len(ic):
            c.err("bad_value", jp(ptr, "instrument_classes"))
    pt = r.get("product_tolerance", _MISSING)
    pp = jp(ptr, "product_tolerance")
    if pt is not _MISSING and v_obj(c, pt, pp):
        check_keys(c, pt, pp, ("abs_da", "ppm", "source"))
        v_float(c, pt.get("abs_da", _MISSING), jp(pp, "abs_da"), positive=True)
        v_float(c, pt.get("ppm", _MISSING), jp(pp, "ppm"))
        v_enum(c, pt.get("source", _MISSING), jp(pp, "source"), ("p2", "msfinder"))
    v_float(c, r.get("precursor_tolerance_da", _MISSING), jp(ptr, "precursor_tolerance_da"), positive=True)


def v_examples(c, v, ptr, min_len):
    if not v_list(c, v, ptr, min_len):
        return
    for i, x in enumerate(v):
        p = jp(ptr, i)
        if v_obj(c, x, p):
            check_keys(c, x, p, ("pubchem_cid", "smiles"))
            v_int(c, x.get("pubchem_cid", _MISSING), jp(p, "pubchem_cid"), minimum=1)
            v_str(c, x.get("smiles", _MISSING), jp(p, "smiles"), maxlen=2000)


def v_rule_schema(c, doc, stem, el):
    """Phase A: schema, types, enums, patterns, formula grammar (design 3.6 item 1)."""
    check_keys(c, doc, "", RULE_KEYS)
    v_format_header(c, doc, "msfinder-diagnostic-rule")
    rid = doc.get("rule_id", _MISSING)
    if v_str(c, rid, "/rule_id", pattern=_RULE_ID_RE) and rid != stem:
        c.err("file_name_mismatch", "/rule_id")
    v_int(c, doc.get("version", _MISSING), "/version", minimum=1)
    src = doc.get("source", _MISSING)
    n_source = None
    if src is not _MISSING:
        v_citations(c, src, "/source", min_len=1)
        n_source = len(src) if isinstance(src, list) else None
    v_str(c, doc.get("transcribed_by", _MISSING), "/transcribed_by", maxlen=NOTE_MAX)
    for k in ("source_checked", "author_reviewed", "measure", "veto_excludes_motif"):
        v_bool(c, doc.get(k, _MISSING), "/" + k)
    if "review" in doc:
        v_review(c, doc["review"], "/review")
    v_enum(c, doc.get("parameters_origin", _MISSING), "/parameters_origin", ORIGINS)
    v_enum(c, doc.get("primary_precondition_mode", _MISSING), "/primary_precondition_mode",
           ("adduct_ambiguous", "mz_only"))
    cal = doc.get("calibrated_on", _MISSING)
    if cal is not _MISSING and v_list(c, cal, "/calibrated_on"):
        for i, x in enumerate(cal):
            v_enum(c, x, jp("/calibrated_on", i), CALIBRATION_CLASSES)
        if len(set(map(str, cal))) != len(cal):
            c.err("bad_value", "/calibrated_on")
    convention = doc.get("convention", _MISSING)
    if not v_enum(c, convention, "/convention", ("p2", "msdial")):
        convention = None
    v_enum(c, doc.get("product_tolerance", _MISSING), "/product_tolerance", ("p2", "msfinder"))
    v_enum(c, doc.get("loss_reference", _MISSING), "/loss_reference", ("hypothesis", "observed"))
    v_str(c, doc.get("notes", _MISSING), "/notes", maxlen=NOTE_MAX, nonempty=False)
    claim = doc.get("claim", _MISSING)
    if claim is not _MISSING and v_obj(c, claim, "/claim"):
        check_keys(c, claim, "/claim", ("motif_id", "claim_term_ids", "claim_level", "non_claims",
                                        "positive_examples", "negative_examples"))
        v_str(c, claim.get("motif_id", _MISSING), "/claim/motif_id", pattern=_ID_RE)
        ids = claim.get("claim_term_ids", _MISSING)
        if ids is not _MISSING and v_list(c, ids, "/claim/claim_term_ids"):
            for i, x in enumerate(ids):
                v_str(c, x, jp("/claim/claim_term_ids", i), pattern=_TERM_ID_RE)
        v_enum(c, claim.get("claim_level", _MISSING), "/claim/claim_level", ("class", "scaffold", "substituent"))
        nc = claim.get("non_claims", _MISSING)
        if nc is not _MISSING and v_list(c, nc, "/claim/non_claims"):
            for i, x in enumerate(nc):
                v_str(c, x, jp("/claim/non_claims", i), maxlen=NOTE_MAX)
        if "positive_examples" in claim:
            v_examples(c, claim["positive_examples"], "/claim/positive_examples", 2)
        if "negative_examples" in claim:
            v_examples(c, claim["negative_examples"], "/claim/negative_examples", 1)
    pre = doc.get("preconditions", _MISSING)
    if pre is not _MISSING and v_obj(c, pre, "/preconditions"):
        check_keys(c, pre, "/preconditions", ("ion_mode", "adducts", "precursor_tolerance_da", "composition",
                                              "priors", "competitors"))
        mode = pre.get("ion_mode", _MISSING)
        if v_enum(c, mode, "/preconditions/ion_mode", ("pos", "neg")) and isinstance(rid, str) \
                and _RULE_ID_RE.match(rid) and rid.rsplit(".", 1)[1] != mode:
            c.err("mode_mismatch", "/preconditions/ion_mode")
        ad = pre.get("adducts", _MISSING)
        if ad is not _MISSING and v_list(c, ad, "/preconditions/adducts", min_len=1):
            for i, x in enumerate(ad):
                v_str(c, x, jp("/preconditions/adducts", i), maxlen=64)
            if len(set(map(str, ad))) != len(ad):
                c.err("bad_value", "/preconditions/adducts")
        v_float(c, pre.get("precursor_tolerance_da", _MISSING), "/preconditions/precursor_tolerance_da",
                positive=True)
        comp = pre.get("composition", _MISSING)
        pc = "/preconditions/composition"
        if comp is not _MISSING and v_obj(c, comp, pc):
            check_keys(c, comp, pc, ("core", "scaffold_list", "units", "unit_totals"))
            core = v_formula(c, comp.get("core", _MISSING), jp(pc, "core"), el.pda, nullable=True)
            sl = comp.get("scaffold_list", _MISSING)
            v_str(c, sl, jp(pc, "scaffold_list"), pattern=_ID_RE, nullable=True)
            if core is None and sl is None:
                c.err("bad_value", pc)
            units = comp.get("units", _MISSING)
            if units is not _MISSING and v_obj(c, units, jp(pc, "units")):
                for name, rng in units.items():
                    if not _UNIT_RE.match(name):
                        c.err("bad_value", jp(pc, "units", name))
                    else:
                        v_range_pair(c, rng, jp(pc, "units", name))
            totals = comp.get("unit_totals", _MISSING)
            if totals is not _MISSING and v_obj(c, totals, jp(pc, "unit_totals")):
                for cls, rng in totals.items():
                    if cls not in UNIT_CLASSES:
                        c.err("bad_value", jp(pc, "unit_totals", cls))
                    else:
                        v_range_pair(c, rng, jp(pc, "unit_totals", cls))
        priors = pre.get("priors", _MISSING)
        if priors is not _MISSING and v_list(c, priors, "/preconditions/priors"):
            for i, x in enumerate(priors):
                p = jp("/preconditions/priors", i)
                if v_obj(c, x, p):
                    check_keys(c, x, p, ("source", "note"))
                    v_token(c, x.get("source", _MISSING), jp(p, "source"), n_source, src_only=True)
                    v_str(c, x.get("note", _MISSING), jp(p, "note"), maxlen=NOTE_MAX)
        comps = pre.get("competitors", _MISSING)
        if comps is not _MISSING and v_list(c, comps, "/preconditions/competitors"):
            for i, x in enumerate(comps):
                p = jp("/preconditions/competitors", i)
                if v_obj(c, x, p):
                    check_keys(c, x, p, ("motif_id", "note"))
                    v_str(c, x.get("motif_id", _MISSING), jp(p, "motif_id"), pattern=_ID_RE)
                    v_str(c, x.get("note", _MISSING), jp(p, "note"), maxlen=NOTE_MAX)
    terms = doc.get("terms", _MISSING)
    if terms is not _MISSING and v_list(c, terms, "/terms", min_len=1):
        for i, t in enumerate(terms):
            v_term(c, t, jp("/terms", i), n_source, convention, el)
    comp_pairs = doc.get("complementary", _MISSING)
    if comp_pairs is not _MISSING and v_list(c, comp_pairs, "/complementary"):
        for i, pair in enumerate(comp_pairs):
            p = jp("/complementary", i)
            if not isinstance(pair, list) or len(pair) != 2:
                c.err("bad_type", p)
                continue
            for j, x in enumerate(pair):
                v_str(c, x, jp(p, j), pattern=_TERM_ID_RE)
    if "logic" in doc:
        stats = {"nodes": 0, "depth": 0}
        walk_logic_schema(c, doc["logic"], "/logic", n_source, stats, 1)
        if stats["nodes"] > MAX_LOGIC_NODES or stats["depth"] > MAX_LOGIC_DEPTH:
            c.err("logic_depth", "/logic")


def place_logic(c, node, ptr, terms, inside_not, parent_kof, used, out_k):
    """Phase B: role placement (design 3.6 item 3). Returns a compiled node."""
    (key, val), = node.items()
    p = jp(ptr, key)
    if key == "term":
        t = terms.get(val)
        if t is None:
            c.err("unknown_term", p)
            return ("term", None)
        used.add(t.index)
        role = t.role
        bad = (role == "upgrade"
               or (role == "counted") != parent_kof
               or (role in VETO_ROLES) != inside_not
               or (role in ("required", "ordering") and inside_not))
        if bad:
            c.err("role_placement", p)
        return ("term", t.index)
    if key in ("all", "any"):
        return (key, [place_logic(c, ch, jp(p, i), terms, inside_not, False, used, out_k) for i, ch in enumerate(val)])
    if key == "not":
        if inside_not:
            c.err("role_placement", p)
        return ("not", place_logic(c, val, p, terms, True, False, used, out_k))
    # k_of
    if inside_not:
        c.err("role_placement", p)
    slot = len(out_k)
    out_k.append((p, val["k"]))
    kids = [place_logic(c, ch, jp(p, "of", i), terms, inside_not, True, used, out_k) for i, ch in enumerate(val["of"])]
    return ("k_of", slot, kids)


def logic_kof_parents(node, parent_slot, out):
    """term index -> list of k_of slots that hold it as a direct child."""
    if node[0] == "term":
        if parent_slot is not None and node[1] is not None:
            out.setdefault(node[1], []).append(parent_slot)
    elif node[0] in ("all", "any"):
        for ch in node[1]:
            logic_kof_parents(ch, None, out)
    elif node[0] == "not":
        logic_kof_parents(node[1], None, out)
    else:
        for ch in node[2]:
            logic_kof_parents(ch, node[1], out)


def multi_member_capable(t):
    if t.kind == "ion_series" or t.scaffold_key is not None:
        return True
    return t.unit_sel is not None and t.unit_sel[0] == "class"


def check_rule_refs(c, doc, file, tables):
    """Phase B: references, roles, claim terms, complementary pairs. Returns RuleInfo or None."""
    n0 = c.n
    el = tables.elements
    pre = doc["preconditions"]
    comp = pre["composition"]
    mode = pre["ion_mode"]
    charge = 1 if mode == "pos" else -1
    terms, by_id = [], {}
    for i, t in enumerate(doc["terms"]):
        ptr = jp("/terms", i)
        term = Term(i, t["term_id"], t["kind"], t["role"], t)
        if term.tid in by_id:
            c.err("duplicate_term", jp(ptr, "term_id"))
        else:
            by_id[term.tid] = term
        terms.append(term)
        if "formula" in t:
            term.formula = parse_formula(t["formula"], el.pda)
        if "plus" in t and t["plus"] is not None:
            term.plus = parse_formula(t["plus"], el.pda)
        if "minus" in t and t["minus"] is not None:
            term.minus = parse_formula(t["minus"], el.pda)
        if "unit" in t:
            u = t["unit"]
            sel = "class" if "class" in u else "name"
            term.unit_sel = (sel, u[sel])
            term.plus = parse_formula(u["plus"], el.pda) if u["plus"] is not None else None
            term.minus = parse_formula(u["minus"], el.pda) if u["minus"] is not None else None
            if sel == "name" and u["name"] not in tables.unit_by_name:
                c.err("unknown_unit", jp(ptr, "unit", "name"))
        if "scaffold_ions" in t:
            term.scaffold_key = t["scaffold_ions"]
        if term.kind in ("aglycone_ion", "ion_series") and "H" not in el.pda:
            c.err("unknown_element", jp(ptr, "form"))  # Y0 forms add or remove hydrogen
        for j, b in enumerate(t["than"] or ()):
            if b == term.tid:
                c.err("bad_value", jp(ptr, "than", j))
    for i, t in enumerate(terms):
        for j, b in enumerate(t.obj["than"] or ()):
            if b not in by_id:
                c.err("unknown_term", jp("/terms", i, "than", j))
    # composition
    pc = "/preconditions/composition"
    core = parse_formula(comp["core"], el.pda) if comp["core"] is not None else {}
    scaffold = None
    if comp["scaffold_list"] is not None:
        scaffold = tables.scaffolds.get(comp["scaffold_list"])
        if scaffold is None:
            c.err("unknown_scaffold_list", jp(pc, "scaffold_list"))
    for name in comp["units"]:
        if name not in tables.unit_by_name:
            c.err("unknown_unit", jp(pc, "units", name))
    for t in terms:
        if t.scaffold_key is None:
            continue
        p = jp("/terms", t.index, "scaffold_ions")
        if comp["scaffold_list"] is None:
            c.err("unknown_scaffold_list", p)
        elif scaffold is not None and any(not e[2].get(t.scaffold_key, {}).get(mode) for e in scaffold.entries):
            c.err("unknown_ion_key", p)
    # adducts tried in adduct_ambiguous order: the rule's list, then common adducts of the polarity (design 4.3)
    tried = []
    for i, name in enumerate(pre["adducts"]):
        a = tables.adduct_by_name.get(name)
        if a is None:
            c.err("unknown_adduct", jp("/preconditions/adducts", i))
        elif a.charge != charge:
            c.err("adduct_polarity", jp("/preconditions/adducts", i))
        else:
            tried.append((a, "rule"))
    listed = set(pre["adducts"])
    tried += [(a, "common") for a in tables.adducts if a.common and a.charge == charge and a.name not in listed]
    # motifs
    claim = doc["claim"]
    if claim["motif_id"] not in tables.motifs:
        c.err("unknown_motif", "/claim/motif_id")
    for i, comp_m in enumerate(pre["competitors"]):
        if comp_m["motif_id"] not in tables.motifs:
            c.err("unknown_motif", jp("/preconditions/competitors", i, "motif_id"))
    # logic placement
    used, k_nodes = set(), []
    logic = place_logic(c, doc["logic"], "/logic", by_id, False, False, used, k_nodes)
    for t in terms:
        if t.role != "upgrade" and t.index not in used and by_id.get(t.tid) is t:
            c.err("unused_term", jp("/terms", t.index))
    # claim terms
    claim_idx, seen = [], set()
    for i, tid in enumerate(claim["claim_term_ids"]):
        p = jp("/claim/claim_term_ids", i)
        t = by_id.get(tid)
        if t is None:
            c.err("unknown_term", p)
        elif tid in seen or t.role not in SUPPORT_ROLES:
            c.err("claim_term_role", p)
        else:
            claim_idx.append(t.index)
        seen.add(tid)
    # complementary pairs, groups by transitive closure
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x

    pair_index, seen_pairs = {}, set()
    for i, (a, b) in enumerate(doc["complementary"]):
        p = jp("/complementary", i)
        ta, tb = by_id.get(a), by_id.get(b)
        for j, (tid, t) in enumerate(((a, ta), (b, tb))):
            if t is None:
                c.err("unknown_term", jp(p, j))
        if ta is None or tb is None:
            continue
        key = frozenset((a, b))
        if a == b or ta.kind == tb.kind or ta.role not in SUPPORT_ROLES or tb.role not in SUPPORT_ROLES \
                or key in seen_pairs:
            c.err("complementary_pair", p)
            continue
        seen_pairs.add(key)
        pair_index.setdefault(ta.index, []).append(i)
        pair_index.setdefault(tb.index, []).append(i)
        parent[find(ta.index)] = find(tb.index)
    groups = {}
    for x in pair_index:
        groups.setdefault(find(x), set()).add(x)
    if c.n != n0:
        return None
    # claim_terms_too_few (design 3.5): >= 2 terms, or one multi-member term counted under k >= 2 everywhere.
    # It does not stop the two-ion passes: the design's example expects it together with claim_single_ion.
    if doc["measure"]:
        ok = len(claim_idx) >= 2
        if len(claim_idx) == 1:
            t = terms[claim_idx[0]]
            parents = {}
            logic_kof_parents(logic, None, parents)
            slots = parents.get(t.index, [])

            def k_at_least_two(k):
                return (is_int(k) and k >= 2) or (isinstance(k, dict) and k["min"] >= 2)

            ok = bool(t.role == "counted" and multi_member_capable(t) and slots
                      and all(k_at_least_two(k_nodes[s][1]) for s in slots))
        if not ok:
            c.err("claim_terms_too_few", "/claim/claim_term_ids")
    unit_ranges =[tuple(comp["units"].get(u.name, (0, 0))) for u in tables.units]
    return RuleInfo(doc["rule_id"], doc, file, mode, charge, terms, by_id, scaffold, core, unit_ranges,
                    {k: tuple(v) for k, v in comp["unit_totals"].items()}, tried, logic, k_nodes,
                    tuple(sorted(set(claim_idx))), list(groups.values()), pair_index)


# ------------------------------------------------------------------------------------- series and hypotheses
def count_series(rule, tables):
    """Number of compositions, by convolution per unit class (never enumerated)."""
    by_class = {}
    for (lo, hi), u in zip(rule.unit_ranges, tables.units):
        by_class.setdefault(u.cls, []).append((lo, hi))
    total = 1
    for cls, ranges in by_class.items():
        ways = {0: 1}
        for lo, hi in ranges:
            nxt = {}
            for s, w in ways.items():
                for n in range(lo, hi + 1):
                    nxt[s + n] = nxt.get(s + n, 0) + w
            ways = nxt
        lo, hi = rule.totals.get(cls, (0, max(ways)))
        total *= sum(w for s, w in ways.items() if lo <= s <= hi)
    for cls, (lo, hi) in rule.totals.items():
        if cls not in by_class and lo > 0:
            return 0
    n_entries = len(rule.scaffold.entries) if rule.scaffold else 1
    return total * n_entries


def iter_compositions(rule, tables):
    """(entry index, entry id, agly composition, unit counts) in design 4.2 order."""
    entries = ([(i, e[0], comp_add(rule.core, e[1])) for i, e in enumerate(rule.scaffold.entries)]
               if rule.scaffold else [(0, "core", dict(rule.core))])
    classes = [u.cls for u in tables.units]
    for idx, eid, agly in entries:
        for counts in itertools.product(*(range(lo, hi + 1) for lo, hi in rule.unit_ranges)):
            ok = True
            for cls, (lo, hi) in rule.totals.items():
                s = sum(n for n, k in zip(counts, classes) if k == cls)
                if not lo <= s <= hi:
                    ok = False
                    break
            if ok:
                yield idx, eid, agly, counts


def units_key(counts, tables):
    return ",".join(f"{u.name}{n}" for u, n in zip(tables.units, counts) if n)


def iter_hypotheses(rule, tables):
    el = tables.elements
    for idx, eid, agly, counts in iter_compositions(rule, tables):
        comp = dict(agly)
        for u, n in zip(tables.units, counts):
            if n:
                comp = comp_add(comp, u.residue, n)
        neutral = el.mass_pda(comp)
        cid = f"{eid}|{units_key(counts, tables)}"
        for a, src in rule.tried:
            pre_comp = comp_add(comp_add({s: n * a.multiplier for s, n in comp.items()}, a.add), a.remove, -1)
            precursor = a.multiplier * neutral + el.mass_pda(a.add) - el.mass_pda(a.remove) \
                - a.charge * el.electron_pda
            h = Hypothesis(eid, idx, counts, a, src, agly, comp, neutral, precursor, cid, f"{cid}|{a.name}")
            yield h, pre_comp


def hypothesis_vars(h, tables):
    v = {"n_units": sum(h.counts), "n_sugars": 0, "n_acyl": 0, "n_hex": 0, "n_dhex": 0, "n_pen": 0, "n_hexa": 0}
    for u, n in zip(tables.units, h.counts):
        if u.cls == "sugar":
            v["n_sugars"] += n
        elif u.cls == "acyl":
            v["n_acyl"] += n
        name_var = {"Hex": "n_hex", "dHex": "n_dhex", "Pen": "n_pen", "HexA": "n_hexa"}.get(u.name)
        if name_var:
            v[name_var] += n
    return v


def k_value(kspec, vars_):
    if is_int(kspec):
        return kspec
    return max(kspec["min"], vars_[kspec["var"]] + kspec["offset"])


def ionise(agly, t, h, el):
    """Aglycone-based ion composition (design 4.5): even/radical, proton or adduct carrier, then plus/minus."""
    form, carrier = t.obj["form"], t.obj["charge_carrier"]
    z = h.adduct.charge
    if carrier == "adduct":
        ion = comp_add(comp_add(agly, h.adduct.add), h.adduct.remove, -1)
    elif form == "even":
        ion = comp_add(agly, {"H": 1}, z)
    else:
        ion = dict(agly) if z == 1 else comp_add(agly, {"H": 2}, -1)
    if t.plus:
        ion = comp_add(ion, t.plus)
    if t.minus:
        ion = comp_add(ion, t.minus, -1)
    return ion


def term_members(t, h, rule, tables):
    """-> [(target pDa, member key, composition, is_loss)] sorted by target then key, equal targets collapsed.
    A window (veto only) yields its centre."""
    el = tables.elements
    z = h.adduct.charge
    raw = []  # (comp, key, is_loss)
    if t.kind in ("product_ion", "loss_from_precursor", "window"):
        is_loss = t.kind == "loss_from_precursor" or (t.kind == "window" and t.obj["as"] == "loss")
        if t.formula is not None:
            raw.append((t.formula, t.obj["formula"], is_loss))
        elif t.unit_sel is not None:
            sel, val = t.unit_sel
            for u, n in zip(tables.units, h.counts):
                if n and (u.cls if sel == "class" else u.name) == val:
                    comp = u.residue
                    if t.plus:
                        comp = comp_add(comp, t.plus)
                    if t.minus:
                        comp = comp_add(comp, t.minus, -1)
                    raw.append((comp, u.name, is_loss))
        elif t.scaffold_key is not None:
            entry = rule.scaffold.entries[h.entry_index]
            for comp in entry[2][t.scaffold_key][rule.mode]:
                raw.append((comp, hill_format(comp), False))
    elif t.kind == "aglycone_ion":
        raw.append((ionise(h.agly, t, h, el), "A", False))
    elif t.kind == "ion_series":
        classes = set(t.obj["unit_classes"])
        sel = [(i, n) for i, (u, n) in enumerate(zip(tables.units, h.counts)) if n and u.cls in classes]
        for sub in itertools.product(*(range(n + 1) for _, n in sel)):
            if not any(sub):
                continue
            s_counts = [0] * len(h.counts)
            for (i, _), k in zip(sel, sub):
                s_counts[i] = k
            if tuple(s_counts) == tuple(h.counts):
                continue
            base = dict(h.agly)
            for u, k in zip(tables.units, s_counts):
                if k:
                    base = comp_add(base, u.residue, k)
            raw.append((ionise(base, t, h, el), units_key(s_counts, tables), False))
    out = []
    for comp, key, is_loss in raw:
        mass = el.mass_pda(comp)
        target = h.precursor_pda - mass if is_loss else mass - z * el.electron_pda
        out.append((target, key, comp, is_loss))
    out.sort(key=lambda m: (m[0], m[1]))
    collapsed = []
    for m in out:
        if collapsed and collapsed[-1][0] == m[0]:
            prev = collapsed[-1]
            collapsed[-1] = (prev[0], prev[1] + "+" + m[1], prev[2], prev[3])
        else:
            collapsed.append(m)
    return collapsed


def eval_logic(node, present, kvals, terms):
    """Two-ion setting (design 3.6 item 6): vetoes absent (every `not` is true), ordering satisfied.
    present: term index -> set of present evidence-unit ids."""
    kind = node[0]
    if kind == "term":
        return bool(present.get(node[1]))
    if kind == "all":
        return all([eval_logic(ch, present, kvals, terms) for ch in node[1]])
    if kind == "any":
        return any([eval_logic(ch, present, kvals, terms) for ch in node[1]])
    if kind == "not":
        return True
    slot, kids = node[1], node[2]
    units, extra = set(), 0
    for ch in kids:
        if ch[0] == "term":
            units |= present.get(ch[1], set())
        elif eval_logic(ch, present, kvals, terms):
            extra += 1
    return len(units) + extra >= max(1, kvals[slot])


def check_series(c, rule, tables):
    """Phase C: series size, counts and masses, complementary members, the two-ion claim, precursor targets
    (design 3.6 items 4-7)."""
    n = count_series(rule, tables)
    pc = "/preconditions/composition"
    if n > MAX_COMPOSITIONS:
        c.err("series_too_large", pc)
        return
    if n == 0:
        c.err("series_empty", pc)
        return
    el = tables.elements
    support = [t for t in rule.terms if t.role in SUPPORT_ROLES]
    targets_checked = [t for t in rule.terms if t.role in SUPPORT_ROLES or t.role == "upgrade"]
    group_of = {}
    for gi, g in enumerate(rule.groups):
        for ti in g:
            group_of[ti] = gi
    claim_set = set(rule.claim_terms)
    measure = rule.doc["measure"]
    flagged = set()
    signatures = {}

    def flag(code, ptr):
        if (code, ptr) not in flagged:
            flagged.add((code, ptr))
            c.err(code, ptr)

    for h, pre_comp in iter_hypotheses(rule, tables):
        if comp_negative(pre_comp):
            flag("negative_count", pc)
            continue
        if not 0 < h.neutral_pda < MASS_LIMIT_PDA or not 0 < h.precursor_pda < MASS_LIMIT_PDA:
            flag("mass_range", pc)
            continue
        members = {}
        bad = False
        for t in rule.terms:
            ms = term_members(t, h, rule, tables)
            for target, _, comp, is_loss in ms:
                mass = el.mass_pda(comp)
                if comp_negative(comp):
                    flag("negative_count", jp("/terms", t.index))
                    bad = True
                elif not 0 < mass < MASS_LIMIT_PDA or (not is_loss and not 0 < target < MASS_LIMIT_PDA):
                    flag("mass_range", jp("/terms", t.index))
                    bad = True
            members[t.index] = ms
        if bad:
            continue
        for t in targets_checked:
            if any(m[0] == h.precursor_pda for m in members[t.index]):
                flag("support_target_is_precursor", jp("/terms", t.index))
        multi = False
        for ti, pairs in rule.pair_index.items():
            if len(members[ti]) > 1:
                multi = True
                for i in pairs:
                    flag("complementary_multi_member", jp("/complementary", i))
        if multi or not measure:
            continue
        zero = tuple(tuple(abs(m[0] - h.precursor_pda) <= PRECURSOR_REGION_PDA for m in members[t.index])
                     for t in support)
        vars_ = hypothesis_vars(h, tables)
        kvals = tuple(k_value(k, vars_) for _, k in rule.k_nodes)
        signatures.setdefault((zero, kvals), h.hypothesis_id)
    if not measure:
        return
    for (zero, kvals) in signatures:
        res = two_ion_passes(rule, support, zero, kvals, group_of, claim_set)
        for code, ptr in res:
            flag(code, ptr)


def two_ion_passes(rule, support, zero, kvals, group_of, claim_set):
    """Both passes for one hypothesis signature -> [(code, pointer)]."""
    unit_of = {}  # (term index, member index) -> unit id
    claim_units, units = set(), []
    group_unit = {}
    for t, zflags in zip(support, zero):
        for mi, _ in enumerate(zflags):
            gi = group_of.get(t.index)
            if gi is not None:
                if gi not in group_unit:
                    group_unit[gi] = len(units)
                    units.append(gi)
                uid = group_unit[gi]
            else:
                uid = len(units)
                units.append((t.index, mi))
            unit_of[(t.index, mi)] = uid
            if t.index in claim_set:
                claim_units.add(uid)

    def present_map(present_units):
        out = {}
        for t, zflags in zip(support, zero):
            s = {unit_of[(t.index, mi)] for mi, z in enumerate(zflags) if not z and unit_of[(t.index, mi)] in present_units}
            out[t.index] = s
        return out

    def holds(present_units):
        return eval_logic(rule.logic, present_map(present_units), kvals, rule.terms)

    errors = []
    all_units = range(len(units))
    if holds(set()) or any(holds({u}) for u in all_units):
        errors.append(("single_ion_claim", "/logic"))
    base = {u for u in all_units if u not in claim_units}
    if holds(base) or any(holds(base | {u}) for u in sorted(claim_units)):
        errors.append(("claim_single_ion", "/claim/claim_term_ids"))
    return errors


def rule_resolved(rule, tables):
    """The `resolved` argument of rule_semantic_sha256 (design 3.2)."""
    doc = rule.doc
    comp = doc["preconditions"]["composition"]
    adducts = [a.row for a, _ in rule.tried]
    units = [u.row for u in tables.units if u.name in comp["units"]]
    motif = tables.motifs.get(doc["claim"]["motif_id"])
    symbols = set()

    def add(f):
        if isinstance(f, str):
            symbols.update(s for s, _ in _TOKEN_RE.findall(f))

    add(comp["core"])
    if rule.scaffold is not None:
        for e in rule.scaffold.doc["entries"]:
            add(e["formula"])
            for modes in e["ions"].values():
                for lst in modes.values():
                    for f in lst:
                        add(f)
    for u in units:
        add(u["residue"])
    for t in doc["terms"]:
        for k in ("formula", "plus", "minus"):
            add(t.get(k))
        if "unit" in t:
            add(t["unit"]["plus"])
            add(t["unit"]["minus"])
    for a in adducts:
        add(a["add"])
        add(a["remove"])
    el = tables.elements
    return {
        "adducts": adducts,
        "units": units,
        "scaffold_list": rule.scaffold.doc if rule.scaffold is not None else None,
        "motif": {k: v for k, v in motif.items() if k != "chemont"} if motif is not None else None,
        "elements": {"electron": el.electron, "masses": {s: el.masses[s] for s in sorted(symbols)}},
    }


def check_provenance_and_review(c, rule, tables):
    """Phase D: parameters_origin (design 3.5) and the review block (design 3.6 item 9)."""
    doc = rule.doc
    origins = []

    def tok(v):
        if v is not None:
            origins.append("as_published" if v.startswith("src:") else v)

    for t in doc["terms"]:
        tok(t["term_source"])
        tok(t["floor_source"])
        tok(t["ratio_source"])
    stack = [doc["logic"]]
    while stack:
        node = stack.pop()
        (key, val), = node.items()
        if key in ("all", "any"):
            stack.extend(val)
        elif key == "not":
            stack.append(val)
        elif key == "k_of":
            tok(val["k_source"])
            stack.extend(val["of"])
    if rule.scaffold is not None:
        origins += [e[3] for e in rule.scaffold.entries]
    comp = doc["preconditions"]["composition"]
    origins += [u.origin for u in tables.units if u.name in comp["units"]]
    origins += [a.origin for a, _ in rule.tried]
    motif = tables.motifs.get(doc["claim"]["motif_id"])
    if motif is not None:
        origins.append("as_published" if motif["origin"] == "source_defined" else motif["origin"])
    derived = max(origins, key=ORIGIN_RANK.__getitem__) if origins else "as_published"
    if doc["parameters_origin"] != derived:
        c.err("origin_inconsistent", "/parameters_origin")
    review = doc["review"]
    if doc["author_reviewed"] != (review is not None):
        c.err("review_missing", "/review")
    if doc["author_reviewed"] and not (doc["source_checked"] and doc["measure"]):
        c.err("review_not_allowed", "/author_reviewed")
    if review is not None:
        if review["rule_semantic_sha256"] != rule_semantic_sha256(doc, rule_resolved(rule, tables)):
            c.err("review_stale", "/review/rule_semantic_sha256")


# ------------------------------------------------------------------------------------------------ the rule set
@dataclass
class LoadResult:
    errors: list
    draft: bool
    file_hashes: dict = field(default_factory=dict)
    ruleset_sha256: str | None = None
    rule_set_version: str | None = None
    n_rules: int = 0
    tables: Tables | None = None
    rules: dict = field(default_factory=dict)  # rule_id -> RuleInfo (schema and references passed)
    staging_files: list = field(default_factory=list)

    @property
    def ok(self):
        return not self.errors

    def semantic_sha256(self, rule_id):
        rule = self.rules[rule_id]
        return rule_semantic_sha256(rule.doc, rule_resolved(rule, self.tables))


def classify(rel, draft):
    parts = rel.split("/")
    top = parts[0]
    is_json = rel.endswith(".json")
    if top.startswith("_") and len(parts) > 1:
        if draft and top in DRAFT_IGNORED_DIRS:
            return "ignore"
        if draft and top in STAGING_DIRS and len(parts) == 2 and is_json:
            return "staging"
        return "unexpected" if is_json else "ignore"
    if not is_json:
        return "ignore"
    if rel in (RULESET_FILE, ELEMENTS_FILE) or rel in TABLE_WRAPPERS:
        return "hashed"
    if len(parts) == 2 and top in ("rules", "scaffolds"):
        return "hashed"
    return "unexpected"


def read_directory(root) -> dict:
    """{posix path: bytes} for every *.json and the history file under root (other files map to None)."""
    root = Path(root)
    files = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            rel = p.relative_to(root).as_posix()
            keep = rel.endswith(".json") or rel == HISTORY_FILE
            files[rel] = p.read_bytes() if keep and not rel.startswith(DRAFT_IGNORED_DIRS[0] + "/") else None
    return files


def unwrap_table(c, doc, path):
    fmt, key = TABLE_WRAPPERS[path]
    check_keys(c, doc, "", ("format", "format_version", key))
    v_format_header(c, doc, fmt)
    rows = doc.get(key, _MISSING)
    if rows is _MISSING:
        return []
    if not isinstance(rows, list):
        c.err("bad_type", "/" + key)
        return []
    return [(i, r) for i, r in enumerate(rows)]


def load_files(files: dict, draft=False, rules_glob=None) -> LoadResult:
    """Validate a rule set given as {posix path relative to RuleSet/: bytes}."""
    errors = []
    res = LoadResult(errors, draft)
    docs, hashed = {}, {}
    staging = {}
    for rel in sorted(files):
        kind = classify(rel, draft)
        if kind == "ignore":
            continue
        if kind == "unexpected":
            errors.append(RuleLoadError("unexpected_file", rel, ""))
            continue
        b = files[rel]
        try:
            doc = load_strict(b) if kind == "hashed" else load_strict_staging(b)
        except RuleFormatError as e:
            errors.append(RuleLoadError(e.code, rel, ""))
            if kind == "hashed":
                hashed[rel] = None
            continue
        scan_refused(Ctx(errors, rel), doc)
        if kind == "hashed":
            hashed[rel] = b
            docs[rel] = doc
        else:
            staging[rel] = doc
    res.staging_files = sorted(staging)
    if all(b is not None for b in hashed.values()):
        res.file_hashes = {p: file_sha256(b) for p, b in hashed.items()}
        res.ruleset_sha256 = ruleset_sha256(hashed)
    # ruleset.json
    required = [RULESET_FILE, ELEMENTS_FILE, ADDUCTS_FILE, UNITS_FILE]
    if not draft:
        required += [MOTIFS_FILE, NESTING_FILE]
    for rel in required:
        if rel not in files:
            errors.append(RuleLoadError("missing_file", rel, ""))
    if RULESET_FILE in docs:
        c = Ctx(errors, RULESET_FILE)
        doc = docs[RULESET_FILE]
        check_keys(c, doc, "", ("format", "format_version", "rule_set_version"))
        v_format_header(c, doc, "msfinder-diagnostic-ruleset")
        ver = doc.get("rule_set_version", _MISSING)
        if v_str(c, ver, "/rule_set_version", pattern=_VERSION_RE):
            res.rule_set_version = ver
    if ELEMENTS_FILE not in docs:
        res.errors = sort_errors(errors)
        return res
    el = load_elements(Ctx(errors, ELEMENTS_FILE), docs[ELEMENTS_FILE])
    if el is None:
        res.errors = sort_errors(errors)
        return res
    tables = Tables(el)
    res.tables = tables
    # row sources: table rows first, then staging files in path order (draft only)
    rows = {"adducts": [], "units": [], "motifs": [], "nesting": []}
    for path, key in ((ADDUCTS_FILE, "adducts"), (UNITS_FILE, "units"), (MOTIFS_FILE, "motifs"),
                      (NESTING_FILE, "nesting")):
        if path in docs:
            for i, r in unwrap_table(Ctx(errors, path), docs[path], path):
                rows[key].append((path, jp("/" + key, i), r))
    for rel, doc in sorted(staging.items()):
        c = Ctx(errors, rel)
        top = rel.split("/")[0]
        if top in ("_motifs", "_nesting"):
            if not isinstance(doc, list):
                c.err("bad_type", "")
                continue
            key = "motifs" if top == "_motifs" else "nesting"
            rows[key] += [(rel, jp("", i), r) for i, r in enumerate(doc)]
        else:
            if not isinstance(doc, dict):
                c.err("bad_type", "")
                continue
            check_keys(c, doc, "", ("units", "adducts"))
            for key in ("units", "adducts"):
                lst = doc.get(key, _MISSING)
                if lst is _MISSING:
                    continue
                if not isinstance(lst, list):
                    c.err("bad_type", "/" + key)
                    continue
                rows[key] += [(rel, jp("/" + key, i), r) for i, r in enumerate(lst)]
    alias_seen = {}
    for rel, ptr, r in rows["adducts"]:
        c = Ctx(errors, rel)
        a = v_adduct_row(c, r, ptr, el)
        if a is None:
            continue
        spellings = [(a.name, jp(ptr, "adduct"))] + [(x, jp(ptr, "aliases", i)) for i, x in enumerate(r["aliases"])]
        dup = False
        for s, p in spellings:
            k = s.replace(" ", "")
            if k in alias_seen:
                c.err("duplicate_id", p)
                dup = True
            alias_seen.setdefault(k, a.name)
        if not dup:
            tables.adducts.append(a)
            tables.adduct_by_name[a.name] = a
    for rel, ptr, r in rows["units"]:
        c = Ctx(errors, rel)
        u = v_unit_row(c, r, ptr, el)
        if u is None:
            continue
        if u.name in tables.unit_by_name:
            c.err("duplicate_id", jp(ptr, "name"))
            continue
        tables.units.append(u)
        tables.unit_by_name[u.name] = u
    if UNITS_FILE in docs:
        for name in REQUIRED_UNITS:
            if name not in tables.unit_by_name:
                errors.append(RuleLoadError("missing_required_unit", UNITS_FILE, "/units"))
                break
    for rel, ptr, r in rows["motifs"]:
        c = Ctx(errors, rel)
        m = v_motif_row(c, r, ptr)
        if m is None:
            continue
        if m["motif_id"] in tables.motifs:
            c.err("duplicate_id", jp(ptr, "motif_id"))
            continue
        tables.motifs[m["motif_id"]] = m
    pairs, edges = set(), {}
    for rel, ptr, r in rows["nesting"]:
        c = Ctx(errors, rel)
        n = v_nesting_row(c, r, ptr)
        if n is None:
            continue
        bad = False
        for k in ("child", "parent"):
            if n[k] not in tables.motifs:
                c.err("unknown_motif", jp(ptr, k))
                bad = True
        if bad:
            continue
        if n["child"] == n["parent"]:
            c.err("nesting_cycle", ptr)
            continue
        if (n["child"], n["parent"]) in pairs:
            c.err("duplicate_id", ptr)
            continue
        pairs.add((n["child"], n["parent"]))
        edges.setdefault(n["child"], set()).add(n["parent"])
        if _reaches(edges, n["parent"], n["child"]):
            c.err("nesting_cycle", ptr)
            edges[n["child"]].discard(n["parent"])
            continue
        tables.nesting.append(n)
    for rel in sorted(docs):
        if rel.startswith("scaffolds/"):
            sl = v_scaffold_list(Ctx(errors, rel), docs[rel], rel[len("scaffolds/"):-len(".json")], el)
            if sl is not None:
                tables.scaffolds[sl.list_id] = sl
    # rules
    for rel in sorted(docs):
        if not rel.startswith("rules/"):
            continue
        stem = rel[len("rules/"):-len(".json")]
        if rules_glob and not (fnmatch.fnmatchcase(stem, rules_glob) or fnmatch.fnmatchcase(stem + ".json", rules_glob)):
            continue
        res.n_rules += 1
        c = Ctx(errors, rel)
        doc = docs[rel]
        n0 = c.n
        v_rule_schema(c, doc, stem, el)
        if c.n != n0 or any(e.file == rel for e in errors):
            continue
        rule = check_rule_refs(c, doc, rel, tables)
        if rule is None:
            continue
        res.rules[rule.rule_id] = rule
        check_series(c, rule, tables)
        check_provenance_and_review(c, rule, tables)
    # version pin (embedded set only; a draft is a directory set, unpinned)
    if not draft and res.rule_set_version is not None and res.ruleset_sha256 is not None:
        check_pin(errors, files.get(HISTORY_FILE), res.rule_set_version, res.ruleset_sha256)
    res.errors = sort_errors(errors)
    return res


def _reaches(edges, start, goal):
    stack, seen = [start], set()
    while stack:
        x = stack.pop()
        if x == goal:
            return True
        if x in seen:
            continue
        seen.add(x)
        stack.extend(edges.get(x, ()))
    return False


_HISTORY_ROW = re.compile(r"([a-z0-9][a-z0-9._-]*)\t([0-9a-f]{64})\t([0-9]{4}-[0-9]{2}-[0-9]{2})\Z")


def check_pin(errors, history, version, sha):
    rows = {}
    if history is not None:
        for line in history.decode("ascii", errors="replace").split("\n"):
            line = line.rstrip("\r")
            if not line:
                continue
            m = _HISTORY_ROW.match(line)
            if not m:
                errors.append(RuleLoadError("history_format", HISTORY_FILE, ""))
                continue
            try:
                datetime.date.fromisoformat(m.group(3))
            except ValueError:
                errors.append(RuleLoadError("history_format", HISTORY_FILE, ""))
                continue
            if m.group(1) in rows:
                errors.append(RuleLoadError("history_format", HISTORY_FILE, ""))
                continue
            rows[m.group(1)] = m.group(2)
    if version not in rows:
        errors.append(RuleLoadError("ruleset_unpinned", RULESET_FILE, "/rule_set_version"))
    elif rows[version] != sha:
        errors.append(RuleLoadError("ruleset_sha_mismatch", RULESET_FILE, "/rule_set_version"))


def validate(root, draft=False, rules_glob=None) -> LoadResult:
    return load_files(read_directory(root), draft=draft, rules_glob=rules_glob)


# ------------------------------------------------------------------------------------------------------- CLI
def fmt6(x):
    s = f"{x:.6f}"
    return s[1:] if s.startswith("-") and set(s[1:]) <= set("0.") else s


def describe_masses(res, rule_id):
    """Every hypothesis precursor m/z and term member target of one rule, as text lines."""
    rule = res.rules[rule_id]
    tables = res.tables
    el = tables.elements
    lines = [f"rule\t{rule_id}\tcompositions={count_series(rule, tables)}\tadducts_tried="
             + ",".join(f"{a.name}({src})" for a, src in rule.tried)]
    if count_series(rule, tables) > MAX_COMPOSITIONS:
        lines.append("series_too_large")
        return lines
    for h, pre_comp in iter_hypotheses(rule, tables):
        lines.append(f"H\t{h.hypothesis_id}\tadduct_source={h.adduct_source}\tneutral={hill_format(h.comp)}"
                     f"\tM={fmt6(mz_of(h.neutral_pda))}\tprecursor_mz={fmt6(mz_of(h.precursor_pda))}")
        for t in rule.terms:
            for target, key, comp, is_loss in term_members(t, h, rule, tables):
                region = "\tprecursor_region" if abs(target - h.precursor_pda) <= PRECURSOR_REGION_PDA else ""
                if t.kind == "window":
                    lines.append(f"  W\t{t.tid}\t{t.role}\tcentre_mz={fmt6(mz_of(target))}"
                                 f"\thalf_width_da={t.obj['half_width_da']}")
                elif is_loss:
                    lines.append(f"  T\t{t.tid}\t{t.role}\t{t.kind}\t{key}\tloss={fmt6(mz_of(el.mass_pda(comp)))}"
                                 f"\ttarget_mz={fmt6(mz_of(target))}{region}")
                else:
                    lines.append(f"  T\t{t.tid}\t{t.role}\t{t.kind}\t{key}\tion={hill_format(comp)}"
                                 f"\ttarget_mz={fmt6(mz_of(target))}{region}")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rule_format.py", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="validate a rule set directory")
    v.add_argument("ruleset_dir")
    v.add_argument("--draft", action="store_true", help="merge _motifs, _nesting, _requests; skip the version pin")
    v.add_argument("--rules", help="validate only rule files whose stem or name matches this glob")
    f = sub.add_parser("fmt", help="rewrite a file in storage form")
    f.add_argument("file")
    s = sub.add_parser("sha", help="print the rule-set sha256 and per-file hashes")
    s.add_argument("ruleset_dir")
    s.add_argument("--draft", action="store_true")
    m = sub.add_parser("masses", help="print every hypothesis precursor m/z and term target of a rule")
    m.add_argument("ruleset_dir")
    m.add_argument("--rule", required=True)
    m.add_argument("--draft", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "fmt":
        path = Path(args.file)
        try:
            obj = load_strict_staging(path.read_bytes())
            changed = write_storage(path, obj)
        except RuleFormatError as e:
            print(f"{e.code}\t{path.as_posix()}\t")
            return 1
        print(("formatted " if changed else "unchanged ") + path.as_posix())
        return 0
    root = Path(args.ruleset_dir)
    if not root.is_dir():
        print(f"not a directory: {root}", file=sys.stderr)
        return 2
    if args.cmd == "validate":
        res = validate(root, draft=args.draft, rules_glob=args.rules)
        for e in res.errors:
            print(e)
        print(f"ERRORS n={len(res.errors)}" if res.errors else f"OK n_rules={res.n_rules}")
        return 1 if res.errors else 0
    res = validate(root, draft=args.draft)
    if args.cmd == "sha":
        for e in res.errors:
            print(e)
        if res.ruleset_sha256 is None:
            print(f"ERRORS n={len(res.errors)}")
            return 1
        print(f"ruleset_sha256\t{res.ruleset_sha256}")
        for p, h in sorted(res.file_hashes.items()):
            print(f"{h}\t{p}")
        for rid in sorted(res.rules):
            print(f"rule_semantic_sha256\t{res.semantic_sha256(rid)}\t{rid}")
        for p in res.staging_files:
            print(f"not_hashed\t{p}")
        return 1 if res.errors else 0
    # masses
    if args.rule not in res.rules:
        for e in res.errors:
            print(e, file=sys.stderr)
        print(f"rule not loaded: {args.rule}", file=sys.stderr)
        return 1
    for e in res.errors:
        print(e, file=sys.stderr)
    for line in describe_masses(res, args.rule):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
