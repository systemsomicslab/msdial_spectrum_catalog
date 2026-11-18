"""WS8 P3 candidate normaliser: literature rule candidates into one table (stdlib only, Python 3.10).

Reads the rule-candidate JSON files given with --inputs (label=path), maps kinds, families, ion modes,
adducts and instruments onto one vocabulary through the TSV tables beside this file, splits the formula
field into term rows (one per formula token, polarity and adduct), recomputes every mass from its formula,
checks printed values against it, and groups term rows across papers by (family, kind, ion mode, adduct,
formula), with a roll-up over adducts. Nothing is dropped: values that could not be mapped are written to
unmapped.tsv and every input record appears exactly once in records.tsv.
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
FAMILY_OVERRIDE = HERE / "family-override-p3.tsv"
TERM_OVERRIDE = HERE / "term-override-p3.tsv"
STRATUM_OVERRIDE = HERE / "stratum-override-p3.tsv"
MAP_FILES = (KIND_MAP, FAMILY_MAP, FAMILY_OVERRIDE, INSTRUMENT_MAP, PAPER_MAP, TERM_OVERRIDE, STRATUM_OVERRIDE)

KINDS = ("product_ion", "neutral_loss", "ion_ratio", "absence", "veto", "sequence", "precondition",
         "adduct_behaviour", "manual")
SCOPES = ("in_scope", "out_of_scope", "generic")
INSTRUMENT_CLASSES = ("qtof", "orbitrap_ft", "ion_trap", "qqq_or_qtrap", "mixed", "compilation", "unknown")
ACTIVATION_NAMES = ("strip", "trap", "cid", "beam", "beam_gas", "ft_cell", "energy", "instrument_trap",
                    "instrument_beam", "in_source", "in_source_strip", "tandem", "source_voltage")
# Collision regime of a record: in-source activation only, resonance CID in an ion trap, beam-type CID in a
# collision cell (QTOF, triple quadrupole, HCD), activation in the FT-ICR cell, more than one, or not known.
REGIMES = ("in_source", "trap", "beam", "ft_cell", "mixed", "unknown")
ION_MODES = ("positive", "negative", "both", "unknown")
# Kinds whose formula tokens may name the absent ion, the observed one, an isobar or a competitor: their term
# rows are flagged role_unverified and kept out of groups until a reader says which token is absent.
NEGATIVE_KINDS = ("absence", "veto")
# Token roles. A precursor is the precursor or intact molecule; not_observed is a formula the text names as not
# observed ('the analogues would be ...'); printed_as is a notation quoted as printed ('printed as C7H12O6'),
# which the extractor does not vouch for. All three are kept out of groups.
ROLES = ("precursor", "not_observed", "printed_as")
UNOBSERVED_ROLES = ("not_observed", "printed_as")

# Element masses. tools/ws8_p2/p2_diagnostics.py:73-77 holds no table: formula_mass() sums RDKit's
# GetMostCommonIsotopeMass at run time. These are the values RDKit 2026.09.2 returns (checked equal, difference
# 0.0, with the interpreter P2 records in its run.json, 2026-11-18), copied so that this tool stays stdlib only.
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
# Electron mass as in tools/ws8_p2/p2_diagnostics.py:77 (a literal there).
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


def load_family_override(path=FAMILY_OVERRIDE, family_map=None):
    """(source, id) -> {family_original, family, family_group, scope, reason}; the target family must be a
    unified family of the family map, which gives its group and scope."""
    family_map = load_family_map() if family_map is None else family_map
    unified = {v["family"]: v for v in family_map.values()}
    table = {}
    for row in read_table(path):
        key = (row["source"], row["id"])
        if key in table:
            raise NormaliseError(f"family override: duplicate row {key}")
        if row["family"] not in unified:
            raise NormaliseError(f"family override: {row['family']!r} is not a family of the family map")
        target = unified[row["family"]]
        table[key] = {"family_original": row["family_original"], "family": row["family"],
                      "family_group": target["family_group"], "scope": target["scope"], "reason": row["reason"]}
    return table


def load_instrument_map(path=INSTRUMENT_MAP):
    """{'rules': [(name, regex, class, strat, energy_only_strat)], 'activation': {name: regex}}.
    A rule whose strat is 'activation' takes its stratum from instrument_activation()."""
    defines, rules, activation = {}, [], {}

    def expand(pattern):
        return re.sub(r"%([A-Z0-9_]+)%", lambda m: "(?:" + defines[m.group(1)] + ")", pattern)

    for row in read_table(path):
        if row["kind"] == "define":
            defines[row["name"]] = expand(row["pattern"])
        elif row["kind"] == "activation":
            if row["name"] not in ACTIVATION_NAMES or row["name"] in activation:
                raise NormaliseError(f"instrument map: unknown or repeated activation row {row['name']!r}")
            activation[row["name"]] = re.compile(row["pattern"], re.IGNORECASE)
        elif row["kind"] == "rule":
            strat, energy_only = row["instrument_strat"], row.get("energy_only_strat") or ""
            if row["instrument_class"] not in INSTRUMENT_CLASSES or (
                    strat not in INSTRUMENT_CLASSES and strat != "activation"):
                raise NormaliseError(f"instrument map: unknown class in rule {row['name']!r}")
            if (strat == "activation") != bool(energy_only) or (energy_only and energy_only not in INSTRUMENT_CLASSES):
                raise NormaliseError(f"instrument map: rule {row['name']!r} needs energy_only_strat exactly when "
                                     "its strat is activation")
            rules.append((row["name"], re.compile(expand(row["pattern"]), re.IGNORECASE),
                          row["instrument_class"], strat, energy_only))
        else:
            raise NormaliseError(f"instrument map: unknown row kind {row['kind']!r}")
    missing = sorted(set(ACTIVATION_NAMES) - set(activation))
    if missing:
        raise NormaliseError(f"instrument map: activation rows missing: {missing}")
    return {"rules": rules, "activation": activation}


def load_paper_map(path=PAPER_MAP):
    return [(row["source"], re.compile(row["id_pattern"]), normalise_doi(row["doi"])) for row in read_table(path)]


def load_term_override(path=TERM_OVERRIDE):
    """(source, id, canonical formula) -> {role, reason}: the role of one formula token of one record, where the
    text says it but the tokeniser cannot read it. A row whose record is present but has no such token stops
    the run (build)."""
    table = {}
    for row in read_table(path):
        key = (row["source"], row["id"], row["formula"])
        if key in table:
            raise NormaliseError(f"term override: duplicate row {key}")
        if row["role"] not in ROLES:
            raise NormaliseError(f"term override: unknown role {row['role']!r}")
        table[key] = {"role": row["role"], "reason": row["reason"]}
    return table


def load_stratum_override(path=STRATUM_OVERRIDE):
    """[(source, id regex, instrument_strat, collision_regime, reason)]: the stratum and regime of the records of
    one paper, where one record names the activation and the others share its method."""
    rows = []
    for row in read_table(path):
        if row["instrument_strat"] not in INSTRUMENT_CLASSES or row["collision_regime"] not in REGIMES:
            raise NormaliseError(f"stratum override: unknown class or regime in {row['id_pattern']!r}")
        rows.append((row["source"], re.compile(row["id_pattern"]), row["instrument_strat"], row["collision_regime"],
                     row["reason"]))
    return rows


def load_maps():
    family = load_family_map()
    return {"kind": load_kind_map(), "family": family, "family_override": load_family_override(family_map=family),
            "instrument": load_instrument_map(), "paper": load_paper_map(), "term_override": load_term_override(),
            "stratum_override": load_stratum_override()}


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
    # a radical mark after the charge: '•', or a '.' directly after it that ends the field or a phrase ([M+H]+.)
    for m in re.finditer(r"\[([^\[\]]+)\]\s*(\d*)([+-])(" + RADICAL + r"|\.(?=$|[\s,;)/]))?", text):
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


def instrument_activation(instrument, collision, condition, activation):
    """The activation a record names, as (marks, origin). marks is a subset of {trap, beam, ft_cell} read from
    the collision field and from explicit phrases in the instrument text ('CID in the linear ion trap'); the
    condition field is read only when those give nothing. In-source and pre-trap activation and negated
    mentions ('not CID') are removed first. When no activation is named but the collision field gives a
    normalised or percentage energy, marks is {energy}.

    A bare 'CID' is trap CID on these hybrids, unless the text names a beam-cell collision gas (N2): the linear
    trap uses helium, so CID with N2 is beam-type activation in the HCD cell (origin 'collision gas')."""
    def clean(text):
        return activation["strip"].sub(" ", normalise_text(text))

    gas = []

    def read(text):
        marks = {n for n in ("trap", "beam", "ft_cell") if activation[n].search(text)}
        if activation["cid"].search(text):
            if activation["beam_gas"].search(text):
                marks.add("beam")
                gas.append(True)
            else:
                marks.add("trap")
        return marks

    marks, origin = read(clean(collision)), []
    if marks:
        origin.append("collision gas" if gas else "collision")
    inst = normalise_text(instrument)
    inst_marks = ({"trap"} if activation["instrument_trap"].search(inst) else set()) | (
        {"beam"} if activation["instrument_beam"].search(inst) else set())
    if inst_marks:
        marks |= inst_marks
        origin.append("instrument")
    if not marks:
        marks = read(clean(condition))
        if marks:
            origin.append("condition")
    if not marks and activation["energy"].search(clean(collision)):
        marks, origin = {"energy"}, ["collision"]
    return marks, "+".join(origin)


def activation_stratum(marks, energy_only):
    """trap CID -> ion_trap; HCD, beam-type or FT-cell activation -> orbitrap_ft; both -> mixed; a bare energy
    -> the rule's energy_only_strat; nothing -> mixed."""
    trap, other = "trap" in marks, bool(marks & {"beam", "ft_cell"})
    if trap and other:
        return "mixed"
    if trap:
        return "ion_trap"
    if other:
        return "orbitrap_ft"
    if "energy" in marks:
        return energy_only
    return "mixed"


def classify_instrument(raw, imap, collision=None, condition=None):
    """Return (instrument_class, instrument_strat, rule, activation) or None when no rule matches. activation is
    '' unless the rule derives its stratum from the activation, then 'trap', 'beam+trap' ... with its origin."""
    if raw is None or not str(raw).strip():
        return "unknown", "unknown", "missing", ""
    text = normalise_text(raw)
    for name, pattern, klass, strat, energy_only in imap["rules"]:
        if pattern.search(text):
            if strat != "activation":
                return klass, strat, name, ""
            marks, origin = instrument_activation(text, collision, condition, imap["activation"])
            label = "+".join(sorted(marks)) or "none"
            return klass, activation_stratum(marks, energy_only), name, label + (f" ({origin})" if origin else "")
    return None


# A list of numbers with optional units: '10-80 V', '100 V and 250 V', '20, 40, 60 V', '1.0-4.0 kV'.
NUMBER_LIST = (r"[-+]?\d+(?:\.\d+)?(?:\s*(?:%|eV|kV|V)\b)?"
               r"(?:\s*(?:-|to|,|and|or|/)\s*[-+]?\d+(?:\.\d+)?(?:\s*(?:%|eV|kV|V)\b)?)*")
VOLTAGE_RE = re.compile(r"\d\s*(?:eV|kV|V)\b")


def _source_voltage_only(text, act):
    """The collision field names a source voltage and no other voltage or energy ('cone voltage 30 V', not
    'SE 13 V, SG 10 V; fragmentor 200 V')."""
    if not act["source_voltage"].search(text):
        return False
    phrase = re.compile(r"(?:" + act["source_voltage"].pattern + r")(?:\s+(?:voltage|potential))?"
                        r"(?:\s+\d+(?=\s*=))?\s*(?:=|of|at)?\s*(?:" + NUMBER_LIST + r")?", re.IGNORECASE)
    return not VOLTAGE_RE.search(phrase.sub(" ", text))


def _marks_regime(marks, strat):
    kinds = marks & {"trap", "beam", "ft_cell"}
    if len(kinds) > 1:
        return "mixed"
    if kinds:
        return next(iter(kinds))
    if "energy" in marks and strat == "ion_trap":
        return "trap"
    return "unknown"


def collision_regime(inst, instrument, collision, condition, imap):
    """(regime, basis, in_source_mentioned) of a record; inst is its classify_instrument() result.

    in_source when the collision field names in-source activation (in-source CID, source-induced dissociation,
    front-end CID ...) and, once those phrases are removed, no collision-cell or trap activation, or when it names
    only a source voltage (cone, fragmentor, orifice ...) and no other voltage or energy. Otherwise the regime of
    the fragmenting instrument: the
    activation the record names on Orbitrap/FT hybrids with an ion trap; beam for QTOF, triple quadrupole and
    Orbitrap strata (ft_cell for IRMPD/SORI in an FT-ICR cell); trap for the ion-trap stratum; mixed for the
    mixed stratum; for compilations and unknown instruments only an explicit trap, HCD or FT-cell activation
    gives a regime. in_source_mentioned marks a record that also names in-source activation (pseudo-MS3)."""
    act = imap["activation"]
    text = normalise_text(collision)
    in_src = bool(act["in_source"].search(text))
    tandem = bool(act["tandem"].search(act["in_source_strip"].sub(" ", text)))
    if in_src and not tandem:
        return "in_source", "collision: in-source only", True
    if not in_src and not tandem and _source_voltage_only(text, act):
        return "in_source", "collision: source voltage only", False
    klass, strat, rule, label = inst
    if label:
        marks, _ = instrument_activation(instrument, collision, condition, act)
        regime = _marks_regime(marks, strat)
        return regime, "activation" if regime != "unknown" else "activation not named", in_src
    if strat in ("qtof", "qqq_or_qtrap"):
        return "beam", "instrument", in_src
    cleaned = act["strip"].sub(" ", text)
    if strat == "orbitrap_ft":
        if act["ft_cell"].search(cleaned) and not act["beam"].search(cleaned):
            return "ft_cell", "collision", in_src
        return "beam", "instrument", in_src
    if strat == "ion_trap":
        return "trap", "instrument", in_src
    if strat == "mixed":
        return "mixed", "instrument", in_src
    explicit = {n for n in ("trap", "beam", "ft_cell") if act[n].search(cleaned)}
    if len(explicit) == 1:
        return next(iter(explicit)), "collision", in_src
    return "unknown", "instrument", in_src


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


CHARGE_RE = re.compile(r"([+-])(?![A-Za-z0-9(\[])| (\d)([+-])(?![A-Za-z0-9(\[])")


def _read_charge(text, pos):
    """Charge written at text[pos]: returns (z, radical_mark, dot_after_charge, end).

    A radical mark is '•' anywhere around the charge, '*' or '.' directly before the charge ('C7H6O.-',
    'C19H15NO4*+'), '*' directly after it ('C8H15N-*'), or '*' after an uncharged body ('CH3*'). A '.' after
    the charge may end a sentence, so it is returned as dot_after_charge and counts only for an odd-electron
    ion."""
    radical = False
    if pos < len(text) and (text[pos] == RADICAL or (text[pos] in "*." and CHARGE_RE.match(text, pos + 1))):
        radical, pos = True, pos + 1
    z = 0
    m = CHARGE_RE.match(text, pos)
    if m:
        if m.group(1):
            z = 1 if m.group(1) == "+" else -1
        else:
            z = int(m.group(2)) * (1 if m.group(3) == "+" else -1)
        pos = m.end()
    dot = False
    if pos < len(text) and text[pos] == RADICAL:
        radical, pos = True, pos + 1
    elif pos < len(text) and text[pos] == "*" and (z or not radical) and (
            pos + 1 == len(text) or not text[pos + 1].isalnum()):
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


# Words and notations that say a formula is the precursor (or the intact molecule), not a product or a loss.
PRECURSOR_WORDS = (r"precursors?(?:\s+ions?)?|parent(?:\s+ions?)?|molecular\s+ions?|(?:de)?protonated\s+molecules?|"
                   r"pseudomolecular\s+ions?")
PRECURSOR_ADDUCT = (r"\[M(?:[+-](?:H|Na|K|NH4|Li|HCOO|CH3COO|Cl|2H|FA-H|2Na-H))?\]\d*[+-]" + RADICAL + r"?"
                    r"|(?<![\w\]])M[+-]" + RADICAL + r"?(?![\w\[(])|(?<![\w\]])MH[+-](?!\w)"
                    r"|(?<![\w\]])\(M[+-](?:H|Na|K|NH4)\)[+-]")
PRECURSOR_BEFORE_RE = re.compile(r"(?:(?i:" + PRECURSOR_WORDS + r")|" + PRECURSOR_ADDUCT + r")$")
PRECURSOR_LABEL_RE = re.compile(r"\s*(?:(?i:" + PRECURSOR_WORDS + r")\b|(?:" + PRECURSOR_ADDUCT + r")"
                                r"(?=\s*$|\s*[,;]|\s+(?:theor|calc|m/z|\d)))")
# A heading that introduces a list of precursors: 'Precursor formulas from Table 1: A, B, C'.
PRECURSOR_HEADER_RE = re.compile(r"\b(?:precursors|(?:precursor|parent)\s+(?:formulas?|ions|compositions?))\b"
                                 r"[^:;]{0,30}:", re.IGNORECASE)
# Words between a precursor notation and the token it introduces: '[M+H]+. at m/z 369.1686 (C18H27NO7)'.
CONNECTOR_TAIL_RE = re.compile(r"(?:\s+(?:ion|ions|observed|seen|detected|found|appearing|appears))*\s+at$",
                               re.IGNORECASE)
CLAUSE_SEPARATOR_RE = re.compile(r";|->|→")
# A full stop that ends a sentence: after a word, number, charge or ')', before a capitalised word or the end, not
# after an abbreviation.
SENTENCE_END_RE = re.compile(r"(?<=[A-Za-z0-9)+-])(?<!calc)(?<!theor)(?<!approx)(?<!\bca)(?<!\bcf)(?<!\bvs)"
                             r"(?<!e\.g)(?<!i\.e)(?<!\bal)(?<!\bFig)(?<!\bno)\.(?:\s+(?=[A-Z][a-z])|\s*$)")
LIST_GAP_RE = re.compile(r"(?:[\s,()]|m/z|\d+(?:\.\d+)?)*")
# Formulas the text names as not observed: 'the QA-path analogues would be C7H7O4- and C6H7O2-', and
# 'C21H27O+ would be 295.2056'.
HYPOTHETICAL_RE = re.compile(r"\bwould\s+(?:be|give|yield|appear|form)\b|\bhypothetical\b", re.IGNORECASE)
WOULD_AFTER_RE = re.compile(r"\s*would\s+(?:be|give|yield|appear|form)\b", re.IGNORECASE)
PRINTED_AS_RE = re.compile(r"\b(?:mis)?printed\s+as\s*['\"‘’“”]?\s*$", re.IGNORECASE)
# Isotope-labelled notations, which are not formulas of the unlabelled ion: 'C10H10(15N)+'.
ISOTOPE_PART = r"\(\d{1,3}[A-Z][a-z]?\d*\)"
ISOTOPE_RE = re.compile(r"(?<![A-Za-z0-9])(?:(?:[A-Z][a-z]?\d*)+" + ISOTOPE_PART + r"|" + ISOTOPE_PART
                        + r"(?=[A-Z]))(?:[A-Z][a-z]?\d*|" + ISOTOPE_PART + r")*")
# A token directly before an arrow (with only numbers and m/z between): 'C24H39O5- 407.2803 -> C23H35O2-'.
ARROW_AFTER_RE = re.compile(r"(?:[\s()]|m/z|\d+(?:\.\d+)?)*(?:->|→)")
# Polarity labels next to a token: 'C14H9O2 and C13H9O2 (negative)', 'in both modes', 'positive: ...'.
POLARITY_AFTER_RE = re.compile(r"\(\s*(?:in\s+)?(pos|neg)(?:itive|ative|\.)?(?:\s+(?:ion\s+)?mode|\s+ESI)?\s*[,;)]",
                               re.IGNORECASE)
POLARITY_GAP_RE = re.compile(r"(?:[\s,/()]|\band\b|\bor\b|m/z|\d+(?:\.\d+)?(?![\d.])|[A-Z][A-Za-z0-9]*(?![A-Za-z0-9])"
                             r"[+-]?)*")
BOTH_MODES_RE = re.compile(r"[^;]*?\bin\s+both\s+(?:ion\s+)?(?:modes|polarities)\b", re.IGNORECASE)
POLARITY_HEADER_RE = re.compile(r"\b(positive|negative)(?:\s+(?:ion\s+)?mode|\s+ESI)?\s*:", re.IGNORECASE)
POLARITY_WORD_RE = re.compile(r"\b(?:pos(?:itive)?|neg(?:ative)?)\b|ESI\s*[+-]", re.IGNORECASE)
POLARITY_WINDOW = 80
# 'A+B', 'A + 2B', 'RCOOH + NH3': formulas joined into one loss. A '+' attached to a formula and followed by a
# space or the end is a charge instead.
OPERAND = r"\d*(?:R\d*'?|[A-Z][a-z]?\d*)+"
COMPOSITE_RE = re.compile(r"(?<![A-Za-z0-9])" + OPERAND + r"(?:(?:\+|\s+\+\s+)" + OPERAND + r")+(?![a-z])")
MULTIPLIED_RE = re.compile(r"(?<![A-Za-z0-9.,])(\d+)((?:[A-Z][a-z]?\d*){2,})(?![a-z])")
TEMPLATE_RE = re.compile(r"(?:R\d*'?|[A-Z][a-z]?\d*)+")


def _clause_start(text, pos):
    """Start of the clause holding text[pos]: after the last ';', arrow or sentence end before it."""
    start = 0
    for regex in (CLAUSE_SEPARATOR_RE, SENTENCE_END_RE):
        for m in regex.finditer(text):
            if m.end() > pos:
                break
            start = max(start, m.end())
    return start


def _depth(text, pos):
    """Parenthesis depth at text[pos]."""
    depth = 0
    for ch in text[:pos]:
        if ch == "(":
            depth += 1
        elif ch == ")" and depth:
            depth -= 1
    return depth


def _top_level_clause_start(text, pos):
    """Start of the clause holding text[pos], counting only ';' outside parentheses."""
    start, depth = 0, 0
    for i, ch in enumerate(text[:pos]):
        if ch == "(":
            depth += 1
        elif ch == ")" and depth:
            depth -= 1
        elif ch == ";" and depth == 0:
            start = i + 1
    return start


def _scope_end(text, pos):
    """End of the parenthetical or clause holding text[pos]: the ')' that closes it, or a ';' or sentence end
    outside nested parentheses."""
    depth = 0
    for i in range(pos, len(text)):
        ch = text[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            if depth == 0:
                return i
            depth -= 1
        elif depth == 0 and (ch == ";" or (ch == "." and SENTENCE_END_RE.match(text, i))):
            return i
    return len(text)


def precursor_role(text, start, label):
    """('precursor', evidence) when the text says the token at text[start] is the precursor or the intact
    molecule: its clause starts with 'precursors' or holds a heading such as 'Precursor formulas ...:', the words
    just before it are 'precursor', 'parent', 'molecular ion', '[M-H]-', 'M+', 'MH+' ... (numbers, 'm/z' and '('
    in between are skipped, and 'at' with 'ion', 'observed' ... after a notation: '[M+H]+. at m/z 369.1686 ('), or
    its parenthetical label starts with one of them. A number followed by a colon in between ends the reach:
    in 'precursor m/z 353: C7H11O6-' the number is the precursor. Clauses end at ';', an arrow or a sentence end.
    Otherwise ('', '')."""
    clause = text[_clause_start(text, start):start]
    if re.match(r"\s*precursors\b", clause, re.IGNORECASE) or PRECURSOR_HEADER_RE.search(clause):
        return "precursor", "clause"
    tail = re.sub(r"(?:[\s(:=]|m/z|\d+(?:\.\d+)?)+$", "", clause)
    # 'precursor m/z 353: C7H11O6-': the number names the precursor and the colon introduces its products
    numbered_heading = re.search(r"\d\s*:", clause[len(tail):])
    if PRECURSOR_BEFORE_RE.search(tail) and not numbered_heading:
        return "precursor", "before"
    connector = CONNECTOR_TAIL_RE.search(tail)
    if (connector and not numbered_heading
            and PRECURSOR_BEFORE_RE.search(re.sub(r"[.•]$", "", tail[:connector.start()]))):
        return "precursor", "before"
    if label and PRECURSOR_LABEL_RE.match(label):
        return "precursor", "label"
    return "", ""


def polarity_label(text, start, end):
    """The polarity the text gives a token: a label after it, past a list of other formulas and numbers
    ('C14H9O2 and C13H9O2 (negative)'), 'in both modes' later in its clause, or a heading of its clause
    ('positive: ...'). Returns positive | negative | both | ''."""
    for m in POLARITY_AFTER_RE.finditer(text, end, min(len(text), end + POLARITY_WINDOW)):
        if POLARITY_GAP_RE.fullmatch(text, end, m.start()):
            return "positive" if m.group(1).lower() == "pos" else "negative"
        break
    if BOTH_MODES_RE.match(text, end):
        return "both"
    headers = list(POLARITY_HEADER_RE.finditer(text, _top_level_clause_start(text, start), start))
    if headers:
        return headers[-1].group(1).lower()
    return ""


def merge_polarity(a, b):
    if not a or a == b:
        return b
    if not b:
        return a
    return "both"


def tokenise_formulas(text, from_value=False, leading=True):
    """Formula tokens in a text. Returns (tokens, rejected, composites, isotopes).

    Tokens inside ion notations such as [M-H-CO2]- are not terms; a bracket holding only a formula, followed
    by a charge ([C4H5O2]+), is one. Uncharged single-element tokens and short digit-free words (PC, SI ...)
    are not formulas, except as the leading token of a formula field ('O', 'I (calc ...)') when leading is true.
    With from_value, a token must start with carbon and contain hydrogen. Each token carries its position
    (start, end) in the normalised text, a role (precursor, not_observed, printed_as or '', see precursor_role,
    HYPOTHETICAL_RE and PRINTED_AS_RE; a role carries over a list that continues it at the same parenthesis
    depth), composite_part (it only occurs inside a composite), before_arrow and a polarity label.

    rejected lists formula-like fragments that are not terms: strings starting with an element, with a digit,
    that use a symbol that is not an element; templates with an R group (RCOOH); multiplied formulas (2H2O).
    composites lists expressions that join formulas into one loss ('CH2O+2H2O', 'RCOOH + NH3'); their parts
    become terms with composite_part set, the expression itself is not a term. isotopes lists isotope-labelled
    notations ('C10H10(15N)+'), which are not terms and whose parts are not read."""
    text = normalise_text(text)
    tokens, rejected, isotopes = [], [], []
    allow_leading = leading

    def accept(body, counts, z, radical, dot, label, original, start, end, leading=False):
        if body in NOTATIONS:
            return
        if counts is None:
            first = PART_RE.match(body)
            if re.search(r"\d", body) and first and first.group(1) in ATOMIC_NUMBER:
                rejected.append(body)
            elif "R" in body and TEMPLATE_RE.fullmatch(body):
                rest = re.sub(r"R\d*'?", "", body)
                if rest and parse_body(rest) is not None:
                    rejected.append(body)
            return
        if z == 0 and not leading:
            if len(counts) == 1 and body not in SINGLE_ELEMENT_NEUTRALS:
                return
            if not re.search(r"\d", body) and len(body) <= 3 and body not in SHORT_NEUTRALS:
                return
        if from_value and (PART_RE.match(body).group(1) != "C" or "H" not in counts):
            return
        tok = _make_token(body, counts, z, radical, dot, label, original)
        tok["start"], tok["end"] = start, end
        tokens.append(tok)

    masked = list(text)
    for m in ISOTOPE_RE.finditer(text):
        z, radical, dot, end = _read_charge(text, m.end())
        isotopes.append(text[m.start():end])
        for i in range(m.start(), end):
            masked[i] = " "
    for m in re.finditer(r"\[([^\[\]]*)\]", text):
        end = m.end()
        inner = m.group(1).strip()
        if re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", inner):
            z, radical, dot, end2 = _read_charge(text, end)
            if z:
                accept(inner, parse_body(inner), z, radical, dot, _label_after(text, end2), text[m.start():end2],
                       m.start(), end2)
                end = end2
        for i in range(m.start(), end):
            masked[i] = " "
    masked = "".join(masked)

    composites = [(m.start(), m.end(), m.group(0)) for m in COMPOSITE_RE.finditer(masked)
                  if all(_operand_is_formula(op) for op in re.split(r"\s*\+\s*", m.group(0)))]
    if not from_value:
        for m in MULTIPLIED_RE.finditer(masked):
            if parse_body(m.group(2)) is not None:
                rejected.append(m.group(0))

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
        leading = (allow_leading and not from_value and not masked[:start].strip()
                   and (not rest or rest[0] in "(;,"))
        # the label is read from the unmasked text, where ion notations such as [M-H]- are still visible
        accept(body, parse_body(body), z, radical, dot, _label_after(text, end2), masked[start:end2].strip(),
               start, end2, leading)
    tokens.sort(key=lambda t: t["start"])
    for tok in tokens:
        tok["role"], tok["role_evidence"] = precursor_role(text, tok["start"], tok["label"])
        if PRINTED_AS_RE.search(text, max(0, tok["start"] - 60), tok["start"]):
            tok["role"], tok["role_evidence"] = "printed_as", "printed as"
    for prev, tok in zip(tokens, tokens[1:]):
        # a list that continues a role: 'precursor m/z 910 (C36H44O22S2), 928 (C36H46O23S2)',
        # 'fragments printed as C18H13O8 (357), C17H11O7 (327)'
        if (prev["role"] in ("precursor", "printed_as") and not tok["role"]
                and LIST_GAP_RE.fullmatch(text, prev["end"], tok["start"])
                and _depth(text, prev["start"]) == _depth(text, tok["start"])):
            tok["role"], tok["role_evidence"] = prev["role"], "list"
    spans = [(m.end(), _scope_end(text, m.end())) for m in HYPOTHETICAL_RE.finditer(text)]
    for tok in tokens:
        if any(s <= tok["start"] < e for s, e in spans) or WOULD_AFTER_RE.match(text, tok["end"]):
            tok["role"], tok["role_evidence"] = "not_observed", "would be"
        tok["composite_part"] = any(s <= tok["start"] and tok["end"] <= e for s, e, _ in composites)
        tok["before_arrow"] = bool(ARROW_AFTER_RE.match(text, tok["end"]))
        tok["polarity"] = polarity_label(text, tok["start"], tok["end"])
    return tokens, rejected, [c for _, _, c in composites], isotopes


def _operand_is_formula(operand):
    """An operand of a composite: a formula, with an optional multiplier and R groups ('2H2O', 'RCOOH')."""
    body = re.sub(r"R\d*'?", "", re.sub(r"^\d+\s?", "", operand.strip()))
    return bool(body) and parse_body(body) is not None


# ------------------------------------------------------------------------------------------- mass check
# A decimal number; a full stop after it (end of a sentence) is allowed, a further digit group is not.
NUMBER_RE = re.compile(r"(?<![\d.])(\d+)\.(\d+)(?!\d)(?!\.\d)(?!\s*(?:%|ppm|eV|kV|V\b|min\b|mDa|mM|uM|µ|"
                       r"°|ms\b|s\b|h\b))")
# Text allowed between a formula and the number it labels: 'C9H7O3+ 163.0390', 'C9H7O3+ (calc 163.0390)',
# '163.0390 (C9H7O3+)'. A ':' ends the pairing ('hexose + NH3: 369.3505').
PAIR_BEFORE_RE = re.compile(r"[\s(=]*(?:(?:calc(?:d|ulated)?|computed|theor(?:etical)?|exact|m/z)\.?[\s(=]*)?")
PAIR_AFTER_RE = re.compile(r"\s*\(\s*")
PAIR_AFTER_END_RE = re.compile(r"\s*(?:[),;]|$)")


def printed_numbers(*texts):
    """Decimal numbers in the texts as (value, decimals, printed text, paired hill); intensities, ppm, eV and
    times are skipped. paired hill is the Hill formula written right next to the number (before it, or alone in
    parentheses after it) whose own mass lies within NEAR_WINDOW_DA of the number, or None. Tokens of one
    element (I, H2), composite parts (C6H12O6 in 'C6H12O6 + NH3') and leading words are never paired: a
    number next to a formula whose mass is far from it is not that formula's mass."""
    out = []
    for text in texts:
        text = normalise_text(text)
        if not text:
            continue
        tokens = [t for t in tokenise_formulas(text, leading=False)[0]
                  if t["mass"] is not None and len(PART_RE.findall(t["hill"])) > 1 and not t["composite_part"]]
        for m in NUMBER_RE.finditer(text):
            value = float(m.group(0))
            near = [t for t in tokens if abs(t["mass"] - value) <= NEAR_WINDOW_DA]
            paired = None
            for tok in sorted(near, key=lambda t: -t["end"]):
                if tok["end"] <= m.start() and PAIR_BEFORE_RE.fullmatch(text, tok["end"], m.start()):
                    paired = tok["hill"]
                    break
            if paired is None:
                for tok in near:
                    if (tok["start"] >= m.end() and PAIR_AFTER_RE.fullmatch(text, m.end(), tok["start"])
                            and PAIR_AFTER_END_RE.match(text, tok["end"])):
                        paired = tok["hill"]
                        break
            out.append((value, len(m.group(2)), m.group(0), paired))
    return out


def product_tolerance(mass):
    return max(PRODUCT_TOL_DA, PRODUCT_TOL_PPM * 1e-6 * mass)


def mass_check(mass, printed, fallback=(), hill_formula=None):
    """Compare a computed m/z (or neutral mass) with printed decimals (value, decimals, text[, paired hill]).

    Numbers within NEAR_WINDOW_DA of the mass are compared, except a number written next to a different
    formula (it is that formula's mass). The fallback numbers (the note) are used only when the primary ones
    (the value) hold none, so a misprint in the value is not hidden by a recomputed mass in the note. A number
    printed with d decimals matches when |printed - mass| <= tol + 0.5 * 10^-d. When none matches, the flag is
    printed_mismatch if any near number could resolve the tolerance, and lowres_printed if all have fewer than
    LOWRES_DECIMALS decimals and half their last digit is >= tol. A mismatch found in the note, which mostly
    holds masses recomputed by the extractor, is advisory: note_mismatch or note_lowres.
    Status: matched | printed_mismatch | lowres_printed | note_mismatch | note_lowres | no_printed_value |
    no_mass."""
    if mass is None:
        return {"status": "no_mass"}
    tol = product_tolerance(mass)

    def usable(x):
        paired = x[3] if len(x) > 3 else None
        return abs(x[0] - mass) <= NEAR_WINDOW_DA and (paired is None or hill_formula is None
                                                        or paired == hill_formula)

    near, origin = [], ""
    for origin, numbers in (("value", printed), ("note", fallback)):
        near = [x for x in numbers if usable(x)]
        if near:
            break
    if not near:
        return {"status": "no_printed_value", "tolerance": tol}
    nearest = min(near, key=lambda x: (abs(x[0] - mass), -x[1]))
    result = {"tolerance": tol, "printed_from": origin, "printed_near": sorted({x[2] for x in near}),
              "nearest": nearest[2], "nearest_decimals": nearest[1], "nearest_diff": nearest[0] - mass}
    if any(abs(x[0] - mass) <= tol + 0.5 * 10 ** -x[1] for x in near):
        result["status"] = "matched"
        return result
    resolving = [x for x in near if x[1] >= LOWRES_DECIMALS or 0.5 * 10 ** -x[1] < tol]
    status = "printed_mismatch" if resolving else "lowres_printed"
    if origin == "note":
        status = "note_mismatch" if resolving else "note_lowres"
    result["status"] = status
    return result


# --------------------------------------------------------------------------------------------- the build
RECORD_FIELDS = [
    "source", "id", "paper_id", "paper_key", "doi", "family_original", "family", "family_group", "scope",
    "family_override", "kind_original", "kind", "evidence_sense", "manual_queue", "ion_mode_original", "ion_mode",
    "ion_mode_inferred", "adduct_original", "adduct", "adduct_flags", "instrument_original", "instrument_class",
    "instrument_strat", "instrument_strat_computed", "instrument_rule", "instrument_activation",
    "instrument_override", "collision_regime", "collision_regime_basis", "in_source_mentioned", "standards_based",
    "formula_original", "formula_from", "formula_flags", "n_formula_tokens", "n_term_rows", "term_status",
    "class_scope", "value", "condition", "collision", "standards", "location", "note", "citation", "source_text",
    "primary_refs", "author_reviewed",
]
TERM_FIELDS = [
    "term_id", "record_id", "source", "paper_key", "family", "family_group", "scope", "kind", "kind_original",
    "ion_mode", "expanded_from_both", "polarity_label", "adduct", "expanded_from_adducts", "token_original",
    "token_label", "formula_from", "formula", "hill", "charge", "term_type", "role", "role_evidence",
    "odd_electron", "radical_printed", "mass", "unsupported_elements", "flags", "mass_check", "tolerance",
    "printed_from", "nearest_printed", "nearest_diff", "instrument_class", "instrument_strat", "collision_regime",
    "standards_based", "group_exclusion", "group_id", "support_id",
]
MASS_FLAG_FIELDS = ["term_id", "record_id", "source", "family", "kind", "formula", "term_type", "mass", "tolerance",
                    "flag", "advisory", "printed_from", "printed_near", "nearest_printed", "nearest_decimals",
                    "nearest_diff"]
MASS_FLAGS = ("printed_mismatch", "lowres_printed", "note_mismatch", "note_lowres")
ADVISORY_MASS_FLAGS = ("note_mismatch", "note_lowres")
MANUAL_FIELDS = ["source", "id", "paper_key", "family", "family_original", "kind_original", "reason", "ion_mode",
                 "adduct", "instrument_class", "standards_based", "value", "formula_original", "note", "citation"]
UNMAPPED_FIELDS = ["field", "value", "n_records", "record_ids"]
SUPPORT_FIELDS = ["support_id", "family", "family_group", "scope", "kind", "ion_mode", "formula", "term_type", "mass",
                  "n_groups", "adducts", "n_term_rows", "n_records", "n_papers", "n_standards_papers", "sources",
                  "instrument_strata", "collision_regimes", "kind_consistent", "paper_keys", "standards_paper_keys",
                  "group_ids"]


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
    used_stratum_overrides = set()
    term_overrides = maps.get("term_override", {})
    record_keys = {(label, r["id"]) for label, (_, recs) in inputs.items() for r in recs}
    used_term_overrides = set()
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
            override = maps.get("family_override", {}).get((label, rid))
            if override is not None:
                if override["family_original"] != family_original:
                    raise NormaliseError(f"family override for {label}:{rid} expects family "
                                         f"{override['family_original']!r}, the record has {family_original!r}")
                fm = override
            row.update(family_original=family_original, family=fm["family"], family_group=fm["family_group"],
                       scope=fm["scope"], family_override=override["reason"] if override else "")
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
            inst = classify_instrument(rec.get("instrument"), maps["instrument"], rec.get("collision"),
                                       rec.get("condition"))
            if inst is None:
                unmapped[("instrument", str(rec.get("instrument")))].append(rid)
                inst = ("unknown", "unknown", "unmapped", "")
            regime = collision_regime(inst, rec.get("instrument"), rec.get("collision"), rec.get("condition"),
                                      maps["instrument"])
            strat, instrument_override = inst[1], ""
            for source, pattern, o_strat, o_regime, reason in maps.get("stratum_override", ()):
                if source == label and pattern.search(rid):
                    strat, instrument_override = o_strat, reason
                    regime = (o_regime, "override", regime[2])
                    used_stratum_overrides.add(pattern.pattern)
                    break
            row.update(instrument_original=rec.get("instrument"), instrument_class=inst[0],
                       instrument_strat=strat, instrument_strat_computed=inst[1], instrument_rule=inst[2],
                       instrument_activation=inst[3], instrument_override=instrument_override,
                       collision_regime=regime[0], collision_regime_basis=regime[1], in_source_mentioned=regime[2])
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
            formula_flags = []
            if isinstance(formula_field, str) and formula_field.strip():
                term_text = formula_field
                tokens, rejected, composites, isotopes = tokenise_formulas(formula_field)
                formula_from = "formula"
                for fragment in rejected:
                    unmapped[("formula_fragment", fragment)].append(rid)
                for expression in composites:
                    unmapped[("formula_composite", expression)].append(rid)
                for notation in isotopes:
                    unmapped[("formula_isotope_label", notation)].append(rid)
                if rejected:
                    formula_flags.append("fragment_not_a_term")
                if composites:
                    formula_flags.append("composite_expression")
                if isotopes:
                    formula_flags.append("isotope_label_not_a_term")
                if tokens:
                    status = "has_terms"
                else:
                    status = "formula_unparsed"
                    unmapped[("formula", formula_field)].append(rid)
            else:
                term_text = rec.get("value") or ""
                tokens, _, _, isotopes = tokenise_formulas(term_text, from_value=True)
                for notation in isotopes:
                    unmapped[("formula_isotope_label", notation)].append(rid)
                if isotopes:
                    formula_flags.append("isotope_label_not_a_term")
                formula_from = "value" if tokens else ""
                status = "has_terms" if tokens else "no_formula"
            for tok in tokens:
                # in a neutral-loss record, the ion before an arrow is the precursor of the loss
                if km["kind"] == "neutral_loss" and tok["charge"] and tok["before_arrow"] and not tok["role"]:
                    tok["role"], tok["role_evidence"] = "precursor", "arrow"
            unique = {}
            for tok in tokens:
                if tok["formula"] in unique:
                    prev = unique[tok["formula"]]
                    if tok["label"] and tok["label"] not in prev["label"]:
                        prev["label"] = (prev["label"] + " | " + tok["label"]).strip(" |")
                    # a role or a composite part holds only when every occurrence says so; not_observed and
                    # printed_as occurrences of one formula keep the first of the two
                    prev["composite_part"] = prev["composite_part"] and tok["composite_part"]
                    if prev["role"] != tok["role"] and not (prev["role"] in UNOBSERVED_ROLES
                                                            and tok["role"] in UNOBSERVED_ROLES):
                        prev["role"], prev["role_evidence"] = "", ""
                    prev["polarity"] = merge_polarity(prev["polarity"], tok["polarity"])
                    continue
                unique[tok["formula"]] = dict(tok)
            for (o_source, o_id, o_formula), override in term_overrides.items():
                if (o_source, o_id) == (label, rid):
                    if o_formula not in unique:
                        raise NormaliseError(f"term override for {label}:{rid} names {o_formula}, which is not a "
                                             "formula token of the record")
                    unique[o_formula]["role"], unique[o_formula]["role_evidence"] = override["role"], "override"
                    used_term_overrides.add((o_source, o_id, o_formula))
            tokens = list(unique.values())
            for tok in tokens:
                if tok["role"] in UNOBSERVED_ROLES:
                    unmapped[(f"formula_{tok['role']}", tok["formula"])].append(rid)
            mentions_polarity = bool(POLARITY_WORD_RE.search(normalise_text(term_text)))
            printed_value = printed_numbers(rec.get("value"))
            printed_note = printed_numbers(rec.get("note"))
            n_rows = 0
            for k, tok in enumerate(tokens, 1):
                check = mass_check(tok["mass"], printed_value, printed_note, tok["hill"])
                polarity_flags = []
                if mode[0] != "both":
                    pols = [mode[0]]
                elif tok["charge"]:
                    pols = ["positive"] if tok["charge"] > 0 else ["negative"]
                elif tok["polarity"] in ("positive", "negative"):
                    pols = [tok["polarity"]]
                    polarity_flags.append("polarity_from_text")
                else:
                    pols = ["positive", "negative"]
                    if not tok["polarity"] and mentions_polarity:
                        polarity_flags.append("polarity_unverified")
                exclusion = []
                if tok["role"]:
                    exclusion.append(tok["role"])
                if km["kind"] in NEGATIVE_KINDS:
                    exclusion.append("role_unverified")
                if km["kind"] == "manual":
                    exclusion.append("manual")
                if tok["composite_part"]:
                    exclusion.append("composite_part")
                flags = []
                if tok["role"] != "precursor":
                    if km["kind"] == "neutral_loss" and tok["charge"]:
                        flags.append("kind_conflict")
                    if km["kind"] == "product_ion" and not tok["charge"]:
                        flags.append("kind_conflict")
                if km["kind"] == "neutral_loss" and not tok["charge"] and tok["odd_electron"]:
                    flags.append("odd_electron_neutral")
                if tok["unsupported_elements"]:
                    flags.append("unsupported_element")
                if tok["radical_mark_mismatch"]:
                    flags.append("radical_mark_mismatch")
                if check["status"] in MASS_FLAGS:
                    flags.append(check["status"])
                flags += polarity_flags
                flags += [f for f in (*ROLES, "role_unverified", "manual", "composite_part") if f in exclusion]
                for pol in pols:
                    if mode[0] == "both":
                        sign = "+" if pol == "positive" else "-"
                        pol_adducts = [a for a in adducts if adduct_sign(a) in (sign, "")] or ["unstated"]
                    else:
                        pol_adducts = adducts
                    pol_adducts = list(dict.fromkeys(pol_adducts))
                    pol_flags = list(flags)
                    if (pol == "positive" and tok["charge"] < 0) or (pol == "negative" and tok["charge"] > 0):
                        pol_flags.append("charge_mode_conflict")
                    for j, adduct in enumerate(pol_adducts, 1):
                        term_id = (f"{rid}#t{k}" + (f".{pol[:3]}" if mode[0] == "both" else "")
                                   + (f".a{j}" if len(pol_adducts) > 1 else ""))
                        terms.append({
                            "term_id": term_id, "record_id": rid, "source": label, "paper_key": row["paper_key"],
                            "doi": row["doi"], "family": row["family"], "family_group": row["family_group"],
                            "scope": row["scope"], "kind": row["kind"], "kind_original": kind_original,
                            "ion_mode": pol, "expanded_from_both": mode[0] == "both",
                            "polarity_label": tok["polarity"], "adduct": adduct,
                            "expanded_from_adducts": len(pol_adducts) > 1,
                            "token_original": tok["token_original"], "token_label": tok["label"],
                            "formula_from": formula_from, "formula": tok["formula"], "hill": tok["hill"],
                            "charge": tok["charge"], "term_type": "ion" if tok["charge"] else "neutral",
                            "role": tok["role"], "role_evidence": tok["role_evidence"],
                            "odd_electron": tok["odd_electron"], "radical_printed": tok["radical_printed"],
                            "mass": tok["mass"], "unsupported_elements": tok["unsupported_elements"],
                            "flags": pol_flags, "mass_check": check["status"], "check": check,
                            "tolerance": check.get("tolerance"), "printed_from": check.get("printed_from"),
                            "nearest_printed": check.get("nearest"), "nearest_diff": check.get("nearest_diff"),
                            "instrument_class": row["instrument_class"], "instrument_strat": row["instrument_strat"],
                            "collision_regime": row["collision_regime"],
                            "standards_based": row["standards_based"], "group_exclusion": exclusion,
                            "group_id": "", "support_id": "",
                            "secondary": json.loads(row["primary_refs"]) if row["primary_refs"] else [],
                            "citation_text": rec.get("citation") or rec.get("source") or "",
                        })
                        n_rows += 1
            row.update(formula_from=formula_from, formula_flags=formula_flags, n_formula_tokens=len(tokens),
                       n_term_rows=n_rows, term_status=status, manual_queue=row["kind"] == "manual")
            records.append(row)
            if row["kind"] == "manual":
                reason = "kind_other" if kind_original == "other" else "kind_unmapped"
                manual.append({**row, "reason": reason, "adduct": adducts})
    groups = group_terms(terms)
    support = support_rollup(terms, groups)
    overrides = {
        "stratum_override_rows_used": sorted(used_stratum_overrides),
        "stratum_override_rows_unused": sorted(p.pattern for _, p, _, _, _ in maps.get("stratum_override", ())
                                               if p.pattern not in used_stratum_overrides),
        "term_override_rows_used": len(used_term_overrides),
        "term_override_rows_without_record": sorted(":".join(k) for k in term_overrides
                                                    if (k[0], k[1]) not in record_keys),
    }
    return {"records": records, "terms": terms, "groups": groups, "support": support, "manual": manual,
            "unmapped": unmapped, "overrides": overrides}


def group_key(term):
    return (term["family"], term["kind"], term["ion_mode"], term["adduct"], term["formula"])


def support_key(term):
    """The term without its adduct: the roll-up that counts every paper supporting a term."""
    return (term["family"], term["kind"], term["ion_mode"], term["formula"])


def _short_id(prefix, key):
    return prefix + hashlib.sha256(json.dumps(key).encode("utf-8")).hexdigest()[:12]


def group_terms(terms):
    """Group the term rows that carry no group_exclusion by group_key; excluded rows keep group_id ''."""
    by_key = defaultdict(list)
    for t in terms:
        if not t["group_exclusion"]:
            by_key[group_key(t)].append(t)
    groups = []
    for key in sorted(by_key):
        members = by_key[key]
        gid = _short_id("g-", key)
        papers = sorted({t["paper_key"] for t in members if t["paper_key"]})
        std_papers = sorted({t["paper_key"] for t in members if t["paper_key"] and t["standards_based"] is True})
        flag_counts = Counter(f for t in members for f in t["flags"])
        n_conflict = sum(1 for t in members if "kind_conflict" in t["flags"])
        first = members[0]
        for t in members:
            t["group_id"] = gid
        groups.append({
            "group_id": gid, "family": key[0], "kind": key[1], "ion_mode": key[2], "adduct": key[3],
            "formula": key[4],
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
            "collision_regimes": dict(sorted(Counter(t["collision_regime"] for t in members).items())),
            # kind_consistent: no member is a charged neutral loss or an uncharged product ion
            "kind_consistent": n_conflict == 0, "n_kind_conflict_rows": n_conflict,
            "kind_original": sorted({str(t["kind_original"]) for t in members}),
            "expanded_from_both": sum(1 for t in members if t["expanded_from_both"]),
            "expanded_from_adducts": sum(1 for t in members if t["expanded_from_adducts"]),
            "flags": dict(sorted(flag_counts.items())),
        })
    return groups


def support_rollup(terms, groups):
    """Per-term support across adducts: the grouped term rows by (family, kind, ion_mode, formula), with the union
    of papers and standards-based papers over every adduct group, 'unstated' included. Sets support_id on the
    term rows and support_* on each group."""
    by_key = defaultdict(list)
    for t in terms:
        if t["group_id"]:
            by_key[support_key(t)].append(t)
    group_index = {g["group_id"]: g for g in groups}
    rows = []
    for key in sorted(by_key):
        members = by_key[key]
        sid = _short_id("s-", key)
        papers = sorted({t["paper_key"] for t in members if t["paper_key"]})
        std_papers = sorted({t["paper_key"] for t in members if t["paper_key"] and t["standards_based"] is True})
        gids = sorted({t["group_id"] for t in members})
        for t in members:
            t["support_id"] = sid
        for gid in gids:
            group_index[gid].update(support_id=sid, support_n_papers=len(papers),
                                    support_n_standards_papers=len(std_papers), support_n_groups=len(gids))
        first = members[0]
        rows.append({
            "support_id": sid, "family": key[0], "family_group": first["family_group"], "scope": first["scope"],
            "kind": key[1], "ion_mode": key[2], "formula": key[3], "term_type": first["term_type"],
            "mass": first["mass"], "n_groups": len(gids), "adducts": sorted({t["adduct"] for t in members}),
            "n_term_rows": len(members), "n_records": len({t["record_id"] for t in members}),
            "n_papers": len(papers), "n_standards_papers": len(std_papers),
            "sources": sorted({t["source"] for t in members}),
            "instrument_strata": dict(sorted(Counter(t["instrument_strat"] for t in members).items())),
            "collision_regimes": dict(sorted(Counter(t["collision_regime"] for t in members).items())),
            "kind_consistent": not any("kind_conflict" in t["flags"] for t in members),
            "paper_keys": papers, "standards_paper_keys": std_papers, "group_ids": gids,
        })
    return rows


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
        return {"n_records": len(recs), "n_term_rows": len(trms),
                "n_term_rows_grouped": sum(1 for t in trms if t["group_id"]), "n_groups": len(grps),
                "n_terms_supported": len({g.get("support_id") for g in grps}), "n_papers": len(papers),
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
    fields = ["family_group", "family", "scope", "n_records", "n_term_rows", "n_term_rows_grouped", "n_groups",
              "n_terms_supported", "n_papers", "n_standards_papers", "n_records_standards_based", "n_records_no_term"]
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
    write_tsv(out / "term_support.tsv", SUPPORT_FIELDS, result["support"])
    write_tsv(out / "manual_queue.tsv", MANUAL_FIELDS, result["manual"])
    flag_rows, seen = [], set()
    for t in terms:
        # one row per term token: the polarity and adduct copies of a token share its mass check
        key = (t["record_id"], t["formula"])
        if t["mass_check"] in MASS_FLAGS and key not in seen:
            seen.add(key)
            c = t["check"]
            flag_rows.append({**t, "flag": t["mass_check"], "advisory": t["mass_check"] in ADVISORY_MASS_FLAGS,
                              "printed_from": c["printed_from"], "printed_near": c["printed_near"],
                              "nearest_decimals": c["nearest_decimals"]})
    write_tsv(out / "mass_flags.tsv", MASS_FLAG_FIELDS, flag_rows)
    unmapped_rows = [{"field": f, "value": v, "n_records": len(set(ids)), "record_ids": sorted(set(ids))}
                     for (f, v), ids in sorted(result["unmapped"].items())]
    write_tsv(out / "unmapped.tsv", UNMAPPED_FIELDS, unmapped_rows)
    counts = {
        "n_records_in": reconciliation["n_input"], "n_records_out": len(records), "n_terms": len(terms),
        "n_term_tokens": len({(t["record_id"], t["formula"]) for t in terms}),
        "n_terms_grouped": sum(1 for t in terms if t["group_id"]),
        "terms_excluded_from_groups": dict(Counter(e for t in terms for e in t["group_exclusion"])),
        "n_groups": len(groups), "n_term_support_rows": len(result["support"]),
        "n_manual_queue": len(result["manual"]), "n_mass_flags": len(flag_rows),
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
        "instrument_class_by_strat": {f"{c}->{s}": n for (c, s), n in sorted(
            Counter((r["instrument_class"], r["instrument_strat"]) for r in records).items())},
        "instrument_overrides": sum(1 for r in records if r["instrument_override"]),
        "collision_regime": dict(Counter(r["collision_regime"] for r in records)),
        "collision_regime_by_strat": {f"{s}/{g}": n for (s, g), n in sorted(
            Counter((r["instrument_strat"], r["collision_regime"]) for r in records).items())},
        "collision_regime_basis": dict(Counter(r["collision_regime_basis"] for r in records)),
        "in_source_mentioned": sum(1 for r in records if r["in_source_mentioned"]),
        "family_overrides": dict(Counter(r["family"] for r in records if r["family_override"])),
        "formula_flags": dict(Counter(f for r in records for f in r["formula_flags"])),
        "term_flags": dict(Counter(f for t in terms for f in t["flags"])),
        "term_roles": dict(Counter(t["role"] or "term" for t in terms)),
        "mass_check": dict(Counter(t["mass_check"] for t in terms)),
        "terms_expanded_from_both": sum(1 for t in terms if t["expanded_from_both"]),
        "terms_expanded_from_adducts": sum(1 for t in terms if t["expanded_from_adducts"]),
        "adduct_tokens": dict(Counter(a for r in records for a in r["adduct"]).most_common()),
        "groups_by_scope": dict(Counter(g["scope"] for g in groups)),
        "groups_by_kind": dict(Counter(g["kind"] for g in groups)),
        "groups_kind_consistent": sum(1 for g in groups if g["kind_consistent"]),
        "groups_only_kind_conflict": sum(1 for g in groups if g["n_kind_conflict_rows"] == g["n_term_rows"]),
        "overrides": result["overrides"],
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
                 "mapping_files": {p.name: sha256_file(p) for p in MAP_FILES}},
        "environment": {"python": platform.python_version(), "implementation": platform.python_implementation(),
                        "platform": platform.system()},
        "parameters": {
            "element_mass": ELEMENT_MASS, "element_mass_source": "RDKit 2026.09.2 GetMostCommonIsotopeMass, which "
            "tools/ws8_p2/p2_diagnostics.py:73-77 calls at run time; values checked equal to RDKit on 2026-11-18",
            "electron_mass": ELECTRON_MASS,
            "electron_mass_source": "tools/ws8_p2/p2_diagnostics.py:77",
            "mz": "(sum of element masses - z * electron_mass) / |z|; neutrals: sum of element masses",
            "product_tolerance": {"da": PRODUCT_TOL_DA, "ppm": PRODUCT_TOL_PPM, "rule": "max(da, ppm * mass)"},
            "near_window_da": NEAR_WINDOW_DA, "printed_precision_allowance": "0.5 * 10^-decimals",
            "lowres_decimals": LOWRES_DECIMALS,
            "lowres_rule": "a non-matching value with fewer decimals whose half last digit is >= the tolerance",
            "mass_check_sources": "value first, note only when value has no number within the window; a note "
            "mismatch is advisory (note_mismatch, note_lowres); a number written next to a different formula is "
            "not compared",
            "group_key": ["family", "kind", "ion_mode", "adduct", "formula"],
            "support_key": ["family", "kind", "ion_mode", "formula"],
            "term_rows": "one per formula token, polarity (ion mode both) and adduct; in a both record an uncharged "
            "token takes the polarity the text gives it, else both (polarity_unverified when the text names "
            "polarities)",
            "group_exclusion": {"precursor": "the token is the precursor or intact molecule (label, nearby text, a "
                                "'precursor ...:' heading, the ion before an arrow in a neutral-loss record, or "
                                "term-override-p3.tsv)",
                                "not_observed": "the text names the formula as not observed ('would be', "
                                "'hypothetical')",
                                "printed_as": "a notation quoted as printed ('printed as X'), not the extractor's "
                                "formula",
                                "role_unverified": f"kinds {list(NEGATIVE_KINDS)}: which token is absent is not "
                                "assigned", "manual": "kind manual (stage-3 other): no kind until a reader assigns one",
                                "composite_part": "the token occurs only inside a composite expression"},
            "collision_regime": "in_source (in-source activation only, or only a source voltage named), trap, "
            "beam, ft_cell, mixed or unknown, from the collision field and the fragmenting instrument; "
            "stratum-override-p3.tsv sets stratum and regime for the records of one paper",
            "kind_consistent": "a group is kind_consistent when no member carries kind_conflict",
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
