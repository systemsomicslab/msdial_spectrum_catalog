"""WS8 P3 candidate normaliser: literature rule candidates into one table (stdlib only, Python 3.10).

Reads the rule-candidate JSON files given with --inputs (label=path), maps kinds, families, ion modes,
adducts and instruments onto one vocabulary through the TSV tables beside this file, splits the formula
field into term rows, recomputes every mass from its formula, checks printed values against it, and groups
term rows across papers by (family, kind, ion mode, adduct, formula). Nothing is dropped: values that could
not be mapped are written to unmapped.tsv and every input record appears exactly once in records.tsv.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
KIND_MAP = HERE / "kind-map-p3.tsv"
FAMILY_MAP = HERE / "family-map-p3.tsv"
INSTRUMENT_MAP = HERE / "instrument-map-p3.tsv"
PAPER_MAP = HERE / "paper-map-p3.tsv"

KINDS = ("product_ion", "neutral_loss", "ion_ratio", "absence", "veto", "sequence", "precondition",
         "adduct_behaviour", "manual")
SCOPES = ("in_scope", "out_of_scope", "generic")
INSTRUMENT_CLASSES = ("qtof", "orbitrap_ft", "ion_trap", "qqq_or_qtrap", "mixed", "compilation", "unknown")
ION_MODES = ("positive", "negative", "both", "unknown")

# Element masses. tools/ws8_p2/p2_diagnostics.py:76 sums RDKit GetMostCommonIsotopeMass; these are those values
# (RDKit 2026.09.2), copied rather than imported so this tool stays stdlib only.
ELEMENT_MASS = {
    "C": 12.0,
    "H": 1.007825032,
    "N": 14.003074,
    "O": 15.99491462,
    "P": 30.97376163,
    "S": 31.972071,
    "Na": 22.98976928,
    "K": 38.96370668,
    "Cl": 34.96885268,
}
# Electron mass as in tools/ws8_p2/p2_diagnostics.py:77.
ELECTRON_MASS = 0.000549
PRODUCT_TOL_DA = 0.005
PRODUCT_TOL_PPM = 10.0
NEAR_WINDOW_DA = 0.5
LOWRES_DECIMALS = 3

PERIODIC = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
    "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg "
    "Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()
ATOMIC_NUMBER = {e: i + 1 for i, e in enumerate(PERIODIC)}
# Uncharged tokens without a digit are words more often than formulas (PC, SI, IS ...); only these are kept.
SHORT_NEUTRALS = {"CO", "NO", "OH", "HO", "CN", "HCN", "HCl", "CHO", "HCO", "NaCl", "KCl", "CS", "COS", "HNO"}
SINGLE_ELEMENT_NEUTRALS = {"H2", "O2", "N2"}
# Ion notations that parse as formulas but are not (B = nucleobase in [BH2]+).
NOTATIONS = {"BH", "BH2", "BH3"}
STANDARD_ADDUCTS = {
    "[M+H]+", "[M-H]-", "[M]+", "[M]-", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M+Li]+", "[M+HCOO]-",
    "[M+CH3COO]-", "[M+Cl]-", "[M-2H]2-", "[M+2H]2+", "[M+H-H2O]+", "[2M+H]+", "[2M-H]-", "[2M+Na]+",
    "[M+2Na-H]+", "[Agly+H]+", "[Agly-H]-",
}
ADDUCT_ALIASES = {
    "M+FA-H": "M+HCOO", "M+HCOOH-H": "M+HCOO", "M-H+HCOOH": "M+HCOO", "M+HCO2": "M+HCOO",
    "M+CH3COOH-H": "M+CH3COO", "M-H+CH3COOH": "M+CH3COO", "M+CH3CO2": "M+CH3COO",
    "M-H2O+H": "M+H-H2O", "107Ag+M": "M+107Ag", "109Ag+M": "M+109Ag",
}

DASHES = dict.fromkeys(map(ord, "−‐‑‒–—―﹣－"), "-")
DOTS = dict.fromkeys(map(ord, "·∙⋅・"), "•")
SUPERSCRIPT_DIGITS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")
RADICAL = "•"


class NormaliseError(Exception):
    """A refused run: bad arguments, an unusable output folder or an unreadable input."""


class ReconciliationError(NormaliseError):
    """records.tsv does not hold every input record exactly once."""


# ----------------------------------------------------------------------------------------------- utilities
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_table(path):
    """Read a TSV whose lines starting with '#' are comments."""
    with open(path, encoding="utf-8", newline="") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: cell(row.get(k)) for k in fields})


def cell(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6f}"
    if isinstance(value, (list, tuple)):
        return ";".join(cell(v) for v in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, ensure_ascii=False)
    return str(value)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
                          encoding="utf-8")


def normalise_text(text):
    """Unicode minus and dashes to '-', middle dots to the radical dot, superscript charges to plain text."""
    if text is None:
        return ""
    text = str(text).translate(DASHES).translate(DOTS)
    text = re.sub(r"([⁰¹²³⁴-⁹]+)([⁺⁻])",
                  lambda m: " " + m.group(1).translate(SUPERSCRIPT_DIGITS) + m.group(2), text)
    return text.replace("⁺", "+").replace("⁻", "-")


# --------------------------------------------------------------------------------------------- mapping tables
def load_kind_map(path=KIND_MAP):
    table = {}
    for row in read_table(path):
        if row["kind"] not in KINDS:
            raise NormaliseError(f"kind map: unknown kind {row['kind']!r}")
        key = (row["source"], row["kind_original"])
        if key in table:
            raise NormaliseError(f"kind map: duplicate row {key}")
        table[key] = {"kind": row["kind"], "evidence_sense": row["evidence_sense"]}
    return table


def load_family_map(path=FAMILY_MAP):
    table = {}
    for row in read_table(path):
        if row["scope"] not in SCOPES:
            raise NormaliseError(f"family map: unknown scope {row['scope']!r}")
        if not re.fullmatch(r"[a-z0-9_]+", row["family"]) or not re.fullmatch(r"[a-z0-9_]+", row["family_group"]):
            raise NormaliseError(f"family map: label must be a lower-case slug: {row['family']!r}")
        if row["family_original"] in table:
            raise NormaliseError(f"family map: duplicate row {row['family_original']!r}")
        table[row["family_original"]] = {k: row[k] for k in ("family", "family_group", "scope", "reason")}
    scopes = defaultdict(set)
    for v in table.values():
        scopes[v["family"]].add((v["scope"], v["family_group"]))
    clash = sorted(f for f, s in scopes.items() if len(s) > 1)
    if clash:
        raise NormaliseError(f"family map: a unified family has two scopes or groups: {clash}")
    return table


def load_instrument_map(path=INSTRUMENT_MAP):
    defines, rules = {}, []

    def expand(pattern):
        return re.sub(r"%([A-Z0-9_]+)%", lambda m: "(?:" + defines[m.group(1)] + ")", pattern)

    for row in read_table(path):
        if row["kind"] == "define":
            defines[row["name"]] = expand(row["pattern"])
        elif row["kind"] == "rule":
            if row["instrument_class"] not in INSTRUMENT_CLASSES or row["instrument_strat"] not in INSTRUMENT_CLASSES:
                raise NormaliseError(f"instrument map: unknown class in rule {row['name']!r}")
            rules.append((row["name"], re.compile(expand(row["pattern"]), re.IGNORECASE),
                          row["instrument_class"], row["instrument_strat"]))
        else:
            raise NormaliseError(f"instrument map: unknown row kind {row['kind']!r}")
    return rules


def load_paper_map(path=PAPER_MAP):
    return [(row["source"], re.compile(row["id_pattern"]), normalise_doi(row["doi"])) for row in read_table(path)]


def load_maps():
    return {"kind": load_kind_map(), "family": load_family_map(), "instrument": load_instrument_map(),
            "paper": load_paper_map()}


# ------------------------------------------------------------------------------------------ field normalisers
def normalise_doi(text):
    """First DOI in text as 'doi:10.xxx' in lower case, or None."""
    m = re.search(r"\b10\.\d{4,9}/\S+", normalise_text(text))
    if not m:
        return None
    doi = m.group(0)
    while doi and doi[-1] in ".,;:'\"":
        doi = doi[:-1]
    while doi.endswith(")") and doi.count(")") > doi.count("("):
        doi = doi[:-1]
    while doi.endswith("]") and doi.count("]") > doi.count("["):
        doi = doi[:-1]
    return "doi:" + doi.lower()


def normalise_ion_mode(raw):
    """Return (mode, inferred) or None when the value cannot be mapped. A leading polarity wins; a leading
    'not stated' stays unknown even when the extractor added an inference in parentheses."""
    if raw is None or not str(raw).strip():
        return "unknown", False
    text = normalise_text(raw).strip().lower()
    inferred = bool(re.search(r"infer|implied|likely|not printed|not restated", text))
    if re.match(r"(not stated|not printed|not specified|unspecified|unknown)\b", text):
        return "unknown", False
    if re.match(r"(both|dual polarity)\b", text):
        return "both", inferred
    if re.match(r"(positive and negative|negative and positive|positive\s*/\s*negative|negative\s*/\s*positive)\b",
                text):
        return "both", inferred
    if re.match(r"positive only\b", text):
        return "positive", inferred
    if re.match(r"negative only\b", text):
        return "negative", inferred
    if re.match(r"(positive\b[^;]*;\s*negative|negative\b[^;]*;\s*positive)\b", _strip_parentheticals(text)):
        return "both", inferred
    if re.match(r"positive\b", text):
        return "positive", inferred
    if re.match(r"negative\b", text):
        return "negative", inferred
    return None


def _strip_parentheticals(text):
    """Drop (...) outside square brackets; adduct notations such as [(M+Na)-aglycone]+ keep theirs."""
    out, square, paren = [], 0, 0
    for ch in text:
        if ch == "[" and paren == 0:
            square += 1
        elif ch == "]" and paren == 0 and square:
            square -= 1
        elif ch == "(" and square == 0:
            paren += 1
            continue
        elif ch == ")" and square == 0 and paren:
            paren -= 1
            continue
        if paren == 0:
            out.append(ch)
    return "".join(out)


def _canonical_adduct(inner, digits, sign, radical):
    inner = re.sub(r"\s+", "", inner)
    inner = ADDUCT_ALIASES.get(inner, inner)
    inner = re.sub(r"^(?:Agly|Agl|A|aglycone)(?=[+-])", "Agly", inner)
    return f"[{inner}]{digits if digits not in ('', '1') else ''}{sign}{RADICAL if radical else ''}"


def parse_adducts(raw):
    """Return (adducts, flags). adducts is None when nothing could be mapped.

    Only the main text counts: parenthetical remarks are dropped, text after an arrow is a product, and in
    'X from Y' the precursor is Y. 'not stated' (or no value) gives ['unstated']."""
    if raw is None or not str(raw).strip():
        return ["unstated"], []
    text = normalise_text(raw).strip()
    if re.match(r"(not stated|not printed|not specified|unspecified)\b", text, re.IGNORECASE):
        return ["unstated"], []
    flags = []
    if re.search(r"in-source", text, re.IGNORECASE):
        flags.append("in_source")
    text = re.sub(r"(?<!\[)\((M[+-][^()]*)\)\s*(\d*[+-])(?![\w\]])", r"[\1]\2", text)
    text = _strip_parentheticals(text)
    text = re.split(r"->|→", text)[0]
    parts = re.split(r"\bfrom\b", text)
    if len(parts) > 1:
        text = parts[-1]
    found = []
    for m in re.finditer(r"\[([^\[\]]+)\]\s*(\d*)([+-])(" + RADICAL + r")?", text):
        found.append((m.start(), _canonical_adduct(m.group(1), m.group(2), m.group(3), bool(m.group(4)))))
    masked = re.sub(r"\[[^\[\]]*\]\s*\d*[+-]?", lambda m: " " * len(m.group(0)), text)
    for m in re.finditer(r"(?<![\w\]])M([+-])(\.|" + RADICAL + r")?(?![\w\[(])", masked):
        found.append((m.start(), "[M]" + m.group(1) + (RADICAL if m.group(2) else "")))
    adducts = []
    for _, a in sorted(found):
        if a not in adducts:
            adducts.append(a)
    if not adducts:
        if "in_source" in flags:
            return ["in_source"], flags
        return None, flags
    if any(a not in STANDARD_ADDUCTS for a in adducts):
        flags.append("nonstandard")
    return adducts, flags


def adduct_sign(adduct):
    m = re.search(r"\]\d*([+-])" + RADICAL + r"?$", adduct)
    return m.group(1) if m else ""


def classify_instrument(raw, rules):
    """Return (instrument_class, instrument_strat, rule) or None when no rule matches."""
    if raw is None or not str(raw).strip():
        return "unknown", "unknown", "missing"
    text = normalise_text(raw)
    for name, pattern, klass, strat in rules:
        if pattern.search(text):
            return klass, strat, name
    return None


# ---------------------------------------------------------------------------------------------- formulas
BODY_RE = re.compile(r"(?<![A-Za-z0-9])((?:[A-Z][a-z]?\d*)+)")
PART_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_body(body):
    """Element counts of a formula body, or None if a symbol is not an element or a count is zero."""
    counts = Counter()
    for el, n in PART_RE.findall(body):
        if el not in ATOMIC_NUMBER:
            return None
        n = int(n) if n else 1
        if n == 0:
            return None
        counts[el] += n
    return counts


def hill(counts):
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else []) + sorted(e for e in counts if e not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)


def charge_text(z):
    if z == 0:
        return ""
    return (str(abs(z)) if abs(z) > 1 else "") + ("+" if z > 0 else "-")


def canonical_formula(counts, z):
    """Hill order, then the charge. A multiply charged ion is bracketed so that the charge digit cannot be
    read as a count: C9H7O3+, C2H4O2, [C21H18O12]2-."""
    if abs(z) > 1:
        return f"[{hill(counts)}]{charge_text(z)}"
    return hill(counts) + charge_text(z)


def formula_mass(counts, z):
    """m/z of an ion, (sum of masses - z * electron) / |z|; the monoisotopic mass of a neutral; None when an
    element is outside the table."""
    if any(e not in ELEMENT_MASS for e in counts):
        return None
    total = sum(ELEMENT_MASS[e] * n for e, n in counts.items())
    if z == 0:
        return total
    return (total - z * ELECTRON_MASS) / abs(z)


def odd_electron(counts, z):
    return (sum(ATOMIC_NUMBER[e] * n for e, n in counts.items()) - z) % 2 == 1


def _read_charge(text, pos):
    """Charge written at text[pos]: returns (z, radical_dot, dot_after_charge, end)."""
    radical = False
    if pos < len(text) and text[pos] == RADICAL:
        radical, pos = True, pos + 1
    z = 0
    m = re.match(r"([+-])(?![A-Za-z0-9(\[])", text[pos:])
    if m:
        z, pos = (1 if m.group(1) == "+" else -1), pos + 1
    else:
        m = re.match(r" (\d)([+-])(?![A-Za-z0-9(\[])", text[pos:])
        if m:
            z = int(m.group(1)) * (1 if m.group(2) == "+" else -1)
            pos += len(m.group(0))
    dot = False
    if pos < len(text) and text[pos] == RADICAL:
        radical, pos = True, pos + 1
    elif z and pos < len(text) and text[pos] == "." and (pos + 1 == len(text) or text[pos + 1] in " ,;)/"):
        dot = True
    return z, radical, dot, pos


def _label_after(text, pos):
    while pos < len(text) and text[pos] == " ":
        pos += 1
    if pos >= len(text) or text[pos] != "(":
        return ""
    depth = 0
    for i in range(pos, len(text)):
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        if depth == 0:
            return text[pos + 1:i].strip()[:120]
    return text[pos + 1:].strip()[:120]


def _make_token(body, counts, z, radical, dot, label, original):
    odd = odd_electron(counts, z)
    unsupported = sorted(e for e in counts if e not in ELEMENT_MASS)
    return {
        "token_original": original,
        "label": label,
        "hill": hill(counts),
        "charge": z,
        "formula": canonical_formula(counts, z),
        "odd_electron": odd,
        "radical_printed": radical or (dot and odd),
        "radical_mark_mismatch": radical and not odd,
        "unsupported_elements": unsupported,
        "mass": formula_mass(counts, z),
        "body": body,
    }


def tokenise_formulas(text, from_value=False):
    """Formula tokens in a text. Returns (tokens, rejected): rejected lists formula-like strings (starting with
    an element, with a digit) that use a symbol that is not an element. Tokens inside ion notations such as
    [M-H-CO2]- are not terms; a bracket holding only a formula, followed by a charge ([C4H5O2]+), is one.
    Uncharged single-element tokens and short digit-free words (PC, SI ...) are not formulas, except as the
    leading token of a formula field ('O', 'I (calc ...)'). With from_value, a token must start with carbon
    and contain hydrogen."""
    text = normalise_text(text)
    tokens, rejected = [], []

    def accept(body, counts, z, radical, dot, label, original, leading=False):
        if body in NOTATIONS:
            return
        if counts is None:
            first = PART_RE.match(body)
            if re.search(r"\d", body) and first and first.group(1) in ATOMIC_NUMBER:
                rejected.append(body)
            return
        if z == 0 and not leading:
            if len(counts) == 1 and body not in SINGLE_ELEMENT_NEUTRALS:
                return
            if not re.search(r"\d", body) and len(body) <= 3 and body not in SHORT_NEUTRALS:
                return
        if from_value and (PART_RE.match(body).group(1) != "C" or "H" not in counts):
            return
        tokens.append(_make_token(body, counts, z, radical, dot, label, original))

    masked = list(text)
    for m in re.finditer(r"\[([^\[\]]*)\]", text):
        end = m.end()
        inner = m.group(1).strip()
        if re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", inner):
            z, radical, dot, end2 = _read_charge(text, end)
            if z:
                accept(inner, parse_body(inner), z, radical, dot, _label_after(text, end2), text[m.start():end2])
                end = end2
        for i in range(m.start(), end):
            masked[i] = " "
    masked = "".join(masked)

    for m in BODY_RE.finditer(masked):
        body, start, end = m.group(1), m.start(), m.end()
        before = masked[start - 1] if start else ""
        if before == "=":
            continue
        if before == "-" and start >= 2 and masked[start - 2] not in " ([/;,":
            continue
        after = masked[end] if end < len(masked) else ""
        if after == "=" or (after == "(" and end + 1 < len(masked) and masked[end + 1].isupper()):
            continue
        z, radical, dot, end2 = _read_charge(masked, end)
        nxt = masked[end2] if end2 < len(masked) else ""
        if nxt.isalnum():
            continue
        rest = masked[end2:].lstrip()
        leading = not from_value and not masked[:start].strip() and (not rest or rest[0] in "(;,")
        accept(body, parse_body(body), z, radical, dot, _label_after(masked, end2), masked[start:end2].strip(),
               leading)
    return tokens, rejected


# ------------------------------------------------------------------------------------------- mass check
NUMBER_RE = re.compile(r"(?<![\d.])(\d+)\.(\d+)(?![\d.])(?!\s*(?:%|ppm|eV|kV|V\b|min\b|mDa|mM|uM|µ|"
                       r"°|ms\b|s\b|h\b))")


def printed_numbers(*texts):
    """Decimal numbers in the texts as (value, decimals, printed text); intensities, ppm, eV and times are
    skipped."""
    out = []
    for text in texts:
        for m in NUMBER_RE.finditer(normalise_text(text)):
            out.append((float(m.group(0)), len(m.group(2)), m.group(0)))
    return out


def product_tolerance(mass):
    return max(PRODUCT_TOL_DA, PRODUCT_TOL_PPM * 1e-6 * mass)


def mass_check(mass, printed, fallback=()):
    """Compare a computed m/z (or neutral mass) with printed decimals (value, decimals, text).

    Numbers within NEAR_WINDOW_DA of the mass are compared; the fallback numbers (the note) are used only
    when the primary ones (the value) hold none, so a misprint in the value is not hidden by a recomputed
    mass in the note. A number printed with d decimals matches when |printed - mass| <= tol + 0.5 * 10^-d.
    When none matches, the flag is printed_mismatch if any near number could resolve the tolerance, and
    lowres_printed if all have fewer than LOWRES_DECIMALS decimals and half their last digit is >= tol.
    Status: matched | printed_mismatch | lowres_printed | no_printed_value | no_mass."""
    if mass is None:
        return {"status": "no_mass"}
    tol = product_tolerance(mass)
    near, origin = [], ""
    for origin, numbers in (("value", printed), ("note", fallback)):
        near = [x for x in numbers if abs(x[0] - mass) <= NEAR_WINDOW_DA]
        if near:
            break
    if not near:
        return {"status": "no_printed_value", "tolerance": tol}
    nearest = min(near, key=lambda x: (abs(x[0] - mass), -x[1]))
    result = {"tolerance": tol, "printed_from": origin, "printed_near": sorted({x[2] for x in near}),
              "nearest": nearest[2], "nearest_decimals": nearest[1], "nearest_diff": nearest[0] - mass}
    if any(abs(p - mass) <= tol + 0.5 * 10 ** -d for p, d, _ in near):
        result["status"] = "matched"
        return result
    resolving = [x for x in near if x[1] >= LOWRES_DECIMALS or 0.5 * 10 ** -x[1] < tol]
    result["status"] = "printed_mismatch" if resolving else "lowres_printed"
    return result


# --------------------------------------------------------------------------------------------- the build
RECORD_FIELDS = [
    "source", "id", "paper_id", "paper_key", "doi", "family_original", "family", "family_group", "scope",
    "kind_original", "kind", "evidence_sense", "manual_queue", "ion_mode_original", "ion_mode", "ion_mode_inferred",
    "adduct_original", "adduct", "adduct_flags", "instrument_original", "instrument_class", "instrument_strat",
    "instrument_rule", "standards_based", "formula_original", "formula_from", "n_formula_tokens", "n_term_rows",
    "term_status", "class_scope", "value", "condition", "collision", "standards", "location", "note", "citation",
    "source_text", "primary_refs", "author_reviewed",
]
TERM_FIELDS = [
    "term_id", "record_id", "source", "paper_key", "family", "family_group", "scope", "kind", "kind_original",
    "ion_mode", "expanded_from_both", "adduct", "token_original", "token_label", "formula_from", "formula", "hill",
    "charge", "term_type", "odd_electron", "radical_printed", "mass", "unsupported_elements", "flags",
    "mass_check", "tolerance", "printed_from", "nearest_printed", "nearest_diff", "instrument_class",
    "instrument_strat",
    "standards_based", "group_id",
]
MASS_FLAG_FIELDS = ["term_id", "record_id", "source", "family", "kind", "formula", "term_type", "mass", "tolerance",
                    "flag", "printed_from", "printed_near", "nearest_printed", "nearest_decimals", "nearest_diff"]
MANUAL_FIELDS = ["source", "id", "paper_key", "family", "family_original", "kind_original", "reason", "ion_mode",
                 "adduct", "instrument_class", "standards_based", "value", "formula_original", "note", "citation"]
UNMAPPED_FIELDS = ["field", "value", "n_records", "record_ids"]


def load_inputs(pairs):
    """pairs: list of 'label=path'. Returns ordered {label: (path, records)}."""
    inputs = {}
    for pair in pairs:
        if "=" not in pair:
            raise NormaliseError(f"--inputs takes label=path, got {pair!r}")
        label, path = pair.split("=", 1)
        if not re.fullmatch(r"[a-z0-9_]+", label) or label in inputs:
            raise NormaliseError(f"bad or repeated input label {label!r}")
        path = Path(path)
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise NormaliseError(f"cannot read input {label}: {exc}") from exc
        if not isinstance(data, list) or not all(isinstance(r, dict) and r.get("id") for r in data):
            raise NormaliseError(f"input {label} must be a JSON list of records with an id")
        inputs[label] = (path, data)
    if not inputs:
        raise NormaliseError("no inputs")
    return inputs


def check_out_dir(out):
    out = Path(out)
    if out.parent.name != "ws8-p3":
        raise NormaliseError(f"output parent folder must be named ws8-p3, not {out.parent.name!r}")
    if not re.fullmatch(r"candidates-\d{8}", out.name):
        raise NormaliseError(f"output folder must be named candidates-<yyyymmdd>, not {out.name!r}")
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise NormaliseError(f"refusing to write into an existing non-empty {out.name}")
    if not out.parent.is_dir():
        raise NormaliseError("output parent folder does not exist")
    return out


def reconcile(inputs, record_rows):
    """Every input record appears exactly once in record_rows; raises ReconciliationError otherwise."""
    expected = Counter((label, r["id"]) for label, (_, recs) in inputs.items() for r in recs)
    got = Counter((row["source"], row["id"]) for row in record_rows)
    problems = []
    dup_in = sorted(k for k, n in expected.items() if n > 1)
    if dup_in:
        problems.append(f"{len(dup_in)} input ids repeated, e.g. {dup_in[:3]}")
    missing = sorted(set(expected) - set(got))
    if missing:
        problems.append(f"{len(missing)} input records missing from records.tsv, e.g. {missing[:3]}")
    extra = sorted(set(got) - set(expected))
    if extra:
        problems.append(f"{len(extra)} rows not in the inputs, e.g. {extra[:3]}")
    twice = sorted(k for k, n in got.items() if n > 1 and expected.get(k) == 1)
    if twice:
        problems.append(f"{len(twice)} records written more than once, e.g. {twice[:3]}")
    per_source = {label: {"input": len(recs), "records_tsv": sum(1 for r in record_rows if r["source"] == label)}
                  for label, (_, recs) in inputs.items()}
    if problems:
        raise ReconciliationError("reconciliation failed: " + "; ".join(problems))
    return {"status": "ok", "per_source": per_source, "n_input": sum(expected.values()), "n_records": len(record_rows)}


def _std(value):
    return value if isinstance(value, bool) else None


def build(inputs, maps):
    unmapped = defaultdict(list)
    records, terms, manual = [], [], []
    for label, (_, recs) in inputs.items():
        for rec in recs:
            rid = rec["id"]
            row = {"source": label, "id": rid, "paper_id": rec.get("paper_id")}
            # kind
            kind_original = rec.get("kind")
            km = maps["kind"].get((label, kind_original))
            if km is None:
                unmapped[("kind", f"{label}:{kind_original}")].append(rid)
                km = {"kind": "manual", "evidence_sense": "manual"}
            row.update(kind_original=kind_original, kind=km["kind"], evidence_sense=km["evidence_sense"])
            # family
            family_original = rec.get("family")
            fm = maps["family"].get(family_original)
            if fm is None:
                unmapped[("family", str(family_original))].append(rid)
                fm = {"family": str(family_original or ""), "family_group": "", "scope": "unmapped"}
            row.update(family_original=family_original, family=fm["family"], family_group=fm["family_group"],
                       scope=fm["scope"])
            # ion mode
            mode = normalise_ion_mode(rec.get("ion_mode"))
            if mode is None:
                unmapped[("ion_mode", str(rec.get("ion_mode")))].append(rid)
                mode = ("unknown", False)
            row.update(ion_mode_original=rec.get("ion_mode"), ion_mode=mode[0], ion_mode_inferred=mode[1])
            # adducts
            adducts, aflags = parse_adducts(rec.get("adduct"))
            if adducts is None:
                unmapped[("adduct", str(rec.get("adduct")))].append(rid)
                adducts = ["unmapped"]
            signs = {adduct_sign(a) for a in adducts} - {""}
            if (mode[0] == "positive" and "-" in signs) or (mode[0] == "negative" and "+" in signs):
                aflags.append("adduct_mode_conflict")
            row.update(adduct_original=rec.get("adduct"), adduct=adducts, adduct_flags=aflags)
            # instrument
            inst = classify_instrument(rec.get("instrument"), maps["instrument"])
            if inst is None:
                unmapped[("instrument", str(rec.get("instrument")))].append(rid)
                inst = ("unknown", "unknown", "unmapped")
            row.update(instrument_original=rec.get("instrument"), instrument_class=inst[0],
                       instrument_strat=inst[1], instrument_rule=inst[2])
            # paper and citations
            doi = None
            if label in {s for s, _, _ in maps["paper"]}:
                for source, pattern, mapped_doi in maps["paper"]:
                    if source == label and pattern.search(rid):
                        doi = mapped_doi
                        break
            if doi is None:
                doi = normalise_doi(rec.get("citation") or rec.get("source") or "")
            paper_key = rec.get("paper_id") or doi
            if not paper_key:
                unmapped[("paper", rid)].append(rid)
            if not doi:
                unmapped[("doi", str(rec.get("citation") or rec.get("source")))].append(rid)
            row.update(paper_key=paper_key or "", doi=doi or "", standards_based=_std(rec.get("standards_based")))
            # originals
            row.update(class_scope=rec.get("scope"), value=rec.get("value"), condition=rec.get("condition"),
                       collision=rec.get("collision"), standards=rec.get("standards"), location=rec.get("location"),
                       note=rec.get("note"), citation=rec.get("citation"), source_text=rec.get("source"),
                       primary_refs=json.dumps(rec.get("primary_refs"), ensure_ascii=False)
                       if rec.get("primary_refs") is not None else "",
                       author_reviewed=rec.get("author_reviewed"), formula_original=rec.get("formula"))
            # terms
            formula_field = rec.get("formula")
            if isinstance(formula_field, str) and formula_field.strip():
                tokens, rejected = tokenise_formulas(formula_field)
                formula_from = "formula"
                for body in rejected:
                    unmapped[("formula_token", body)].append(rid)
                if tokens:
                    status = "has_terms"
                else:
                    status = "formula_unparsed"
                    unmapped[("formula", formula_field)].append(rid)
            else:
                tokens, _ = tokenise_formulas(rec.get("value") or "", from_value=True)
                formula_from = "value" if tokens else ""
                status = "has_terms" if tokens else "no_formula"
            unique = {}
            for tok in tokens:
                if tok["formula"] in unique:
                    prev = unique[tok["formula"]]
                    if tok["label"] and tok["label"] not in prev["label"]:
                        prev["label"] = (prev["label"] + " | " + tok["label"]).strip(" |")
                    continue
                unique[tok["formula"]] = dict(tok)
            tokens = list(unique.values())
            printed_value = printed_numbers(rec.get("value"))
            printed_note = printed_numbers(rec.get("note"))
            n_rows = 0
            for k, tok in enumerate(tokens, 1):
                check = mass_check(tok["mass"], printed_value, printed_note)
                if mode[0] == "both":
                    pols = ["positive"] if tok["charge"] > 0 else ["negative"] if tok["charge"] < 0 else \
                        ["positive", "negative"]
                else:
                    pols = [mode[0]]
                for pol in pols:
                    if mode[0] == "both":
                        sign = "+" if pol == "positive" else "-"
                        pol_adducts = [a for a in adducts if adduct_sign(a) in (sign, "")] or ["unstated"]
                        term_id = f"{rid}#t{k}.{pol[:3]}"
                    else:
                        pol_adducts = adducts
                        term_id = f"{rid}#t{k}"
                    flags = []
                    if km["kind"] == "neutral_loss" and tok["charge"]:
                        flags.append("kind_conflict")
                    if km["kind"] == "product_ion" and not tok["charge"]:
                        flags.append("kind_conflict")
                    if tok["unsupported_elements"]:
                        flags.append("unsupported_element")
                    if (pol == "positive" and tok["charge"] < 0) or (pol == "negative" and tok["charge"] > 0):
                        flags.append("charge_mode_conflict")
                    if tok["radical_mark_mismatch"]:
                        flags.append("radical_mark_mismatch")
                    if check["status"] in ("printed_mismatch", "lowres_printed"):
                        flags.append(check["status"])
                    terms.append({
                        "term_id": term_id, "record_id": rid, "source": label, "paper_key": row["paper_key"],
                        "doi": row["doi"], "family": row["family"], "family_group": row["family_group"],
                        "scope": row["scope"], "kind": row["kind"], "kind_original": kind_original, "ion_mode": pol,
                        "expanded_from_both": mode[0] == "both", "adduct": sorted(set(pol_adducts)),
                        "token_original": tok["token_original"], "token_label": tok["label"],
                        "formula_from": formula_from, "formula": tok["formula"], "hill": tok["hill"],
                        "charge": tok["charge"], "term_type": "ion" if tok["charge"] else "neutral",
                        "odd_electron": tok["odd_electron"], "radical_printed": tok["radical_printed"],
                        "mass": tok["mass"], "unsupported_elements": tok["unsupported_elements"], "flags": flags,
                        "mass_check": check["status"], "check": check, "tolerance": check.get("tolerance"),
                        "printed_from": check.get("printed_from"),
                        "nearest_printed": check.get("nearest"), "nearest_diff": check.get("nearest_diff"),
                        "instrument_class": row["instrument_class"], "instrument_strat": row["instrument_strat"],
                        "standards_based": row["standards_based"],
                        "secondary": json.loads(row["primary_refs"]) if row["primary_refs"] else [],
                        "citation_text": rec.get("citation") or rec.get("source") or "",
                    })
                    n_rows += 1
            row.update(formula_from=formula_from, n_formula_tokens=len(tokens), n_term_rows=n_rows,
                       term_status=status, manual_queue=row["kind"] == "manual")
            records.append(row)
            if row["kind"] == "manual":
                reason = "kind_other" if kind_original == "other" else "kind_unmapped"
                manual.append({**row, "reason": reason, "adduct": adducts})
    groups = group_terms(terms)
    return {"records": records, "terms": terms, "groups": groups, "manual": manual, "unmapped": unmapped}


def group_key(term):
    return (term["family"], term["kind"], term["ion_mode"], ";".join(term["adduct"]), term["formula"])


def group_terms(terms):
    by_key = defaultdict(list)
    for t in terms:
        by_key[group_key(t)].append(t)
    groups = []
    for key in sorted(by_key):
        members = by_key[key]
        gid = "g-" + hashlib.sha256(json.dumps(key).encode("utf-8")).hexdigest()[:12]
        papers = sorted({t["paper_key"] for t in members if t["paper_key"]})
        std_papers = sorted({t["paper_key"] for t in members if t["paper_key"] and t["standards_based"] is True})
        flag_counts = Counter(f for t in members for f in t["flags"])
        first = members[0]
        for t in members:
            t["group_id"] = gid
        groups.append({
            "group_id": gid, "family": key[0], "kind": key[1], "ion_mode": key[2],
            "adduct": key[3].split(";") if key[3] else [], "formula": key[4],
            "family_group": first["family_group"], "scope": first["scope"], "term_type": first["term_type"],
            "mass": round(first["mass"], 6) if first["mass"] is not None else None,
            "n_term_rows": len(members), "record_ids": sorted({t["record_id"] for t in members}),
            "term_ids": sorted(t["term_id"] for t in members),
            "sources": sorted({t["source"] for t in members}),
            "n_papers": len(papers), "paper_keys": papers,
            "n_standards_papers": len(std_papers), "standards_paper_keys": std_papers,
            "citations": sorted({t["doi"] for t in members if t["doi"]}),
            "citations_without_doi": sorted({t["citation_text"] for t in members if not t["doi"]}),
            "secondary_citations": sorted({s for t in members for s in t["secondary"]}),
            "instrument_classes": dict(sorted(Counter(t["instrument_class"] for t in members).items())),
            "instrument_strata": dict(sorted(Counter(t["instrument_strat"] for t in members).items())),
            "kind_original": sorted({str(t["kind_original"]) for t in members}),
            "expanded_from_both": sum(1 for t in members if t["expanded_from_both"]),
            "flags": dict(sorted(flag_counts.items())),
        })
    return groups


def counts_family_kind(records, sources):
    table = defaultdict(Counter)
    meta = {}
    for r in records:
        table[(r["family"], r["kind"])][r["source"]] += 1
        meta[r["family"]] = (r["family_group"], r["scope"])
    rows = []
    for (family, kind) in sorted(table):
        c = table[(family, kind)]
        rows.append({"family": family, "family_group": meta[family][0], "scope": meta[family][1], "kind": kind,
                     **{s: c[s] for s in sources}, "total": sum(c.values())})
    total = Counter()
    for (family, kind), c in table.items():
        total.update(c)
    rows.append({"family": "TOTAL", "family_group": "", "scope": "", "kind": "",
                 **{s: total[s] for s in sources}, "total": sum(total.values())})
    return ["family", "family_group", "scope", "kind", *sources, "total"], rows


def counts_family_group(records, terms, groups):
    def summarise(recs, trms, grps):
        papers = {r["paper_key"] for r in recs if r["paper_key"]}
        std = {r["paper_key"] for r in recs if r["paper_key"] and r["standards_based"] is True}
        return {"n_records": len(recs), "n_term_rows": len(trms), "n_groups": len(grps), "n_papers": len(papers),
                "n_standards_papers": len(std),
                "n_records_standards_based": sum(1 for r in recs if r["standards_based"] is True),
                "n_records_no_term": sum(1 for r in recs if not r["n_term_rows"])}

    rec_by, term_by, group_by = defaultdict(list), defaultdict(list), defaultdict(list)
    scope = {}
    for r in records:
        rec_by[(r["family_group"], r["family"])].append(r)
        scope[(r["family_group"], r["family"])] = r["scope"]
    for t in terms:
        term_by[(t["family_group"], t["family"])].append(t)
    for g in groups:
        group_by[(g["family_group"], g["family"])].append(g)
    rows = []
    for fg in sorted({k[0] for k in rec_by}):
        keys = sorted(k for k in rec_by if k[0] == fg)
        scopes = sorted({scope[k] for k in keys})
        rows.append({"family_group": fg, "family": "*", "scope": "/".join(scopes),
                     **summarise([r for k in keys for r in rec_by[k]], [t for k in keys for t in term_by[k]],
                                 [g for k in keys for g in group_by[k]])})
        for k in keys:
            rows.append({"family_group": fg, "family": k[1], "scope": scope[k],
                         **summarise(rec_by[k], term_by[k], group_by[k])})
    rows.append({"family_group": "TOTAL", "family": "", "scope": "", **summarise(records, terms, groups)})
    fields = ["family_group", "family", "scope", "n_records", "n_term_rows", "n_groups", "n_papers",
              "n_standards_papers", "n_records_standards_based", "n_records_no_term"]
    return fields, rows


def git_state():
    def run(*args):
        try:
            out = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.SubprocessError):
            return None
        return out.stdout.strip() if out.returncode == 0 else None

    commit = run("rev-parse", "HEAD")
    status = run("status", "--porcelain", "--", "tools/ws8_p3")
    return {"commit": commit or "unknown", "tools_ws8_p3_dirty": bool(status) if status is not None else None}


def write_outputs(out, inputs, result, argv):
    out = Path(out)
    records, terms, groups = result["records"], result["terms"], result["groups"]
    reconciliation = reconcile(inputs, records)
    out.mkdir(exist_ok=True)
    write_tsv(out / "records.tsv", RECORD_FIELDS, records)
    with open(out / "records.tsv", encoding="utf-8", newline="") as f:
        written = list(csv.DictReader(f, delimiter="\t"))
    reconciliation = reconcile(inputs, written)
    reconciliation["checked"] = "in memory and after re-reading records.tsv"
    write_tsv(out / "terms.tsv", TERM_FIELDS, terms)
    write_json(out / "groups.json", groups)
    sources = list(inputs)
    fields, rows = counts_family_kind(records, sources)
    write_tsv(out / "counts_family_kind.tsv", fields, rows)
    fields, rows = counts_family_group(records, terms, groups)
    write_tsv(out / "counts_family_group.tsv", fields, rows)
    write_tsv(out / "manual_queue.tsv", MANUAL_FIELDS, result["manual"])
    flag_rows = []
    for t in terms:
        if t["mass_check"] in ("printed_mismatch", "lowres_printed"):
            c = t["check"]
            flag_rows.append({**t, "flag": t["mass_check"], "printed_from": c["printed_from"],
                              "printed_near": c["printed_near"], "nearest_decimals": c["nearest_decimals"]})
    write_tsv(out / "mass_flags.tsv", MASS_FLAG_FIELDS, flag_rows)
    unmapped_rows = [{"field": f, "value": v, "n_records": len(ids), "record_ids": sorted(set(ids))}
                     for (f, v), ids in sorted(result["unmapped"].items())]
    write_tsv(out / "unmapped.tsv", UNMAPPED_FIELDS, unmapped_rows)
    counts = {
        "n_records_in": reconciliation["n_input"], "n_records_out": len(records), "n_terms": len(terms),
        "n_groups": len(groups), "n_manual_queue": len(result["manual"]), "n_mass_flags": len(flag_rows),
        "mass_flags": dict(Counter(r["flag"] for r in flag_rows)),
        "n_unmapped_values": len(unmapped_rows),
        "unmapped_by_field": dict(Counter(r["field"] for r in unmapped_rows)),
        "unmapped_records_by_field": {f: len({i for r in unmapped_rows if r["field"] == f for i in r["record_ids"]})
                                      for f in sorted({r["field"] for r in unmapped_rows})},
        "term_status": dict(Counter(r["term_status"] for r in records)),
        "formula_from": dict(Counter(r["formula_from"] or "none" for r in records)),
        "kind": dict(Counter(r["kind"] for r in records)),
        "scope": dict(Counter(r["scope"] for r in records)),
        "ion_mode": dict(Counter(r["ion_mode"] for r in records)),
        "instrument_class": dict(Counter(r["instrument_class"] for r in records)),
        "instrument_strat": dict(Counter(r["instrument_strat"] for r in records)),
        "term_flags": dict(Counter(f for t in terms for f in t["flags"])),
        "mass_check": dict(Counter(t["mass_check"] for t in terms)),
        "terms_expanded_from_both": sum(1 for t in terms if t["expanded_from_both"]),
        "adduct_tokens": dict(Counter(a for r in records for a in r["adduct"]).most_common()),
        "groups_by_scope": dict(Counter(g["scope"] for g in groups)),
    }
    doi_papers = defaultdict(set)
    for r in records:
        if r["doi"] and r["paper_key"]:
            doi_papers[r["doi"]].add(r["paper_key"])
    run = {
        "run_id": out.name,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tool": "tools/ws8_p3/normalise_candidates.py",
        "argv": argv,
        "inputs": [{"label": label, "path": str(path), "sha256": sha256_file(path), "n_records": len(recs)}
                   for label, (path, recs) in inputs.items()],
        "code": {**git_state(), "script_sha256": sha256_file(Path(__file__)),
                 "mapping_files": {p.name: sha256_file(p) for p in (KIND_MAP, FAMILY_MAP, INSTRUMENT_MAP, PAPER_MAP)}},
        "environment": {"python": platform.python_version(), "implementation": platform.python_implementation(),
                        "platform": platform.system()},
        "parameters": {
            "element_mass": ELEMENT_MASS, "element_mass_source": "tools/ws8_p2/p2_diagnostics.py:76 (RDKit "
            "GetMostCommonIsotopeMass, 2026.09.2 values)", "electron_mass": ELECTRON_MASS,
            "electron_mass_source": "tools/ws8_p2/p2_diagnostics.py:77",
            "mz": "(sum of element masses - z * electron_mass) / |z|; neutrals: sum of element masses",
            "product_tolerance": {"da": PRODUCT_TOL_DA, "ppm": PRODUCT_TOL_PPM, "rule": "max(da, ppm * mass)"},
            "near_window_da": NEAR_WINDOW_DA, "printed_precision_allowance": "0.5 * 10^-decimals",
            "lowres_decimals": LOWRES_DECIMALS,
            "lowres_rule": "a non-matching value with fewer decimals whose half last digit is >= the tolerance",
            "group_key": ["family", "kind", "ion_mode", "adduct", "formula"],
        },
        "counts": counts,
        "doi_shared_by_paper_ids": {d: sorted(p) for d, p in sorted(doi_papers.items()) if len(p) > 1},
        "reconciliation": reconciliation,
    }
    write_json(out / "run.json", run)
    return run


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--inputs", nargs="+", required=True, metavar="LABEL=PATH",
                        help="rule-candidate JSON files; labels must appear in kind-map-p3.tsv (npr, stage2, stage3)")
    parser.add_argument("--out", required=True, help=".../ws8-p3/candidates-<yyyymmdd>, new or empty")
    args = parser.parse_args(argv)
    try:
        out = check_out_dir(args.out)
        maps = load_maps()
        inputs = load_inputs(args.inputs)
        known = {s for s, _ in maps["kind"]}
        unknown = sorted(set(inputs) - known)
        if unknown:
            raise NormaliseError(f"input labels not in kind-map-p3.tsv: {unknown}")
        result = build(inputs, maps)
        run = write_outputs(out, inputs, result, argv)
    except NormaliseError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    c = run["counts"]
    print(f"{run['run_id']}: {c['n_records_in']} records in, {c['n_records_out']} out, {c['n_terms']} term rows, "
          f"{c['n_groups']} groups, {c['n_manual_queue']} manual, {c['n_mass_flags']} mass flags, "
          f"{c['n_unmapped_values']} unmapped values; reconciliation {run['reconciliation']['status']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
