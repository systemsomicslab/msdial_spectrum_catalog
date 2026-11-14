"""Compare MS-Emblator and CANOPUS predictions with strict reference-library matches, and check the report.

Run with an environment that has RDKit and pandas (structure SMILES are turned into InChIKeys); nothing
else is needed. The truth is the strict reference-library identity MS-DIAL assigned (a Level-2 match),
not an authentic standard, so every rate is "agreement with the strict reference match", never accuracy.

Inputs, all found from the batch run folder:

- plan.json, and for every benchmark export its manifest (*_manifest.local-vs21.tsv), read for the
  reference identity;
- per chunk: name_map.tsv, formula_score.csv, structure_score.csv, and SIRIUS's formula, structure and
  CANOPUS summaries under work/;
- the ClassyFire direct-parent table (InChIKey, Direct parent name, Direct parent ID);
- candidate structure databases, as <path>::<first-block InChIKey column>.

Outputs, in --out:

- benchmark_rows.local-vs21.tsv: one row per spectrum (raw, local; keeps its local-vs21 name);
- benchmark_summary.json and benchmark_summary.md: aggregates, which may be shared. When they name
  reference identities they cite the library by file name and sha256.

A chunk whose driver exit is not 0, or that reports a failed step (for example formula:msbuddy), is left
out of every rate and listed under excluded_chunks; its spectra keep their rows, with the reason in
`excluded`. A spectrum whose reference formula (charged or unparseable) or InChIKey (empty) cannot be
evaluated is left out of that metric's rates and counted.

The script ends with a self-check over the shareable outputs: absolute paths, MSP record text and account
details must not appear (`leak hits: 0`), and every shareable file that names a reference identity must
carry the library citation (`citation: ok`).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# WS6's tool and score identities (WS5 spec section 7); each metric names the one it was measured on.
TOOL_FORMULA = "msemblator.formula"
TOOL_STRUCTURE = "msemblator.structure"
TOOL_CANOPUS = "sirius.canopus"
# score_type is msemblator.<step>.rank_score.<ranker>, the ranker taken from the run records (WS0 ruling,
# 2026-11-13): MS-Emblator loads the first of xgbranker, catboostranker, lgbmranker it finds.
SCORE_TEMPLATE = "msemblator.{step}.rank_score.{ranker}"
POSTERIOR_MSEMBLATOR = "msemblator.softmax_top10"
POSTERIOR_CANOPUS = "sirius.canopus.posterior"

STRUCTURE_K = (1, 5, 10, 100)
FORMULA_K = (1, 3, 5, 10)
CANOPUS_P = (0.5, 0.7, 0.9)
MZ_BINS = ((0, 300, "<300"), (300, 650, "300-650"), (650, float("inf"), ">650"))
REFERENCE_FIELDS = ("reference_name", "reference_inchikey", "reference_formula")

# ------------------------------------------------------------------------------------ leak self-check

_ABSOLUTE_PATH = re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]|\\\\[A-Za-z0-9_.-]+\\|(?<![A-Za-z0-9.])/(?:Users|home|mnt)/")
_LIBRARY_FILE_PATH = re.compile(r"[\\/][^\s\\/]*\.(?:msp|lbm2|lbm)\b", re.IGNORECASE)
_MSP_TEXT = re.compile(r"^\s*(?:Num Peaks|NumPeaks)\s*:|^\s*\d+(?:\.\d+)?[\t ]+\d+(?:\.\d+)?\s*$",
                       re.IGNORECASE | re.MULTILINE)
_ACCOUNT = re.compile(r"(?<!hit)@|\bsid\b", re.IGNORECASE)  # every "@" except the report's own hit@k


def leak_hits(text: str) -> dict[str, int]:
    """Counts of what a shareable report must never hold."""
    return {
        "absolute_or_library_path": len(_ABSOLUTE_PATH.findall(text)) + len(_LIBRARY_FILE_PATH.findall(text)),
        "msp_record_text": len(_MSP_TEXT.findall(text)),
        "account_detail": len(_ACCOUNT.findall(text)),
    }


def names_reference(text: str, identities: list[str]) -> bool:
    return any(identity and identity in text for identity in identities)


def self_check(files: list[Path], identities: list[str],
               citation_marker: str | None) -> tuple[int, list[str], dict]:
    """(leak hit count, files naming a reference identity without the citation, hits per file).

    citation_marker is None when no reference library was recorded: then any file that names a reference
    identity is uncited, whatever fallback text it carries.
    """
    total, uncited, detail = 0, [], {}
    for path in files:
        text = path.read_text(encoding="utf-8")
        hits = leak_hits(text)
        detail[path.name] = hits
        total += sum(hits.values())
        if names_reference(text, identities) and (citation_marker is None or citation_marker not in text):
            uncited.append(path.name)
    return total, uncited, detail


# ------------------------------------------------------------------------------------------ helpers

def wilson(hits: int, n: int, z: float = 1.96) -> dict:
    if n == 0:
        return {"rate": None, "low": None, "high": None, "hits": hits, "n": n}
    p = hits / n
    centre = p + z * z / (2 * n)
    spread = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    denominator = 1 + z * z / n
    return {"rate": round(p, 4), "low": round((centre - spread) / denominator, 4),
            "high": round((centre + spread) / denominator, 4), "hits": hits, "n": n}


_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")


def hill(formula: str) -> str | None:
    """Hill-normalised neutral formula; None for an empty or charged/unparseable one."""
    formula = (formula or "").strip()
    if not formula or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula):
        return None
    counts: dict[str, int] = defaultdict(int)
    for element, number in _ELEMENT.findall(formula):
        counts[element] += int(number) if number else 1
    order = (["C", "H"] if "C" in counts else []) + sorted(e for e in counts if not ("C" in counts and e in "CH"))
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order if counts[e] > 0)


def first_block(inchikey: str | None) -> str:
    return (inchikey or "").strip().upper()[:14]


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def column(row: dict, *names: str) -> str:
    """The first present column among names (case-insensitive)."""
    lowered = {key.lower(): key for key in row}
    for name in names:
        if name.lower() in lowered:
            return row[lowered[name.lower()]] or ""
    return ""


def number(text: str) -> float | None:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


def mz_bin(mz: float | None) -> str:
    if mz is None:
        return "unknown"
    for low, high, label in MZ_BINS:
        if low <= mz < high:
            return label
    return "unknown"


def instrument_family(instrument: str | None) -> str:
    text = (instrument or "").lower()
    if any(word in text for word in ("orbitrap", "exactive", "q exactive", "fusion", "lumos", "exploris")):
        return "Orbitrap"
    if any(word in text for word in ("tripletof", "triple tof", "zenotof", "sciex")):
        return "TripleTOF"
    if any(word in text for word in ("tof", "qtof", "synapt", "xevo", "impact", "maxis", "6546", "6545", "6550",
                                     "lcms-9030", "9030")):
        return "QTOF"
    return "unknown" if not text else "other"


# ------------------------------------------------------------------------------- reference classes

def load_direct_parents(icd: Path) -> tuple[dict[str, tuple[str, str]], dict[str, set[tuple[str, str]]]]:
    full, by_block = {}, defaultdict(set)
    with icd.open("r", encoding="utf-8-sig", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3 or not parts[0]:
                continue
            key, name, identifier = parts[0].strip().upper(), parts[1].strip(), parts[2].strip()
            full[key] = (name, identifier)
            by_block[key[:14]].add((name, identifier))
    return full, by_block


def reference_parent(inchikey: str, full: dict, by_block: dict) -> tuple[str, str, str] | None:
    key = (inchikey or "").strip().upper()
    if key in full:
        return (*full[key], "full_inchikey")
    parents = by_block.get(key[:14], set())
    if len(parents) == 1:
        return (*next(iter(parents)), "first_block_unique")
    return None


def load_db_blocks(specs: list[str]) -> tuple[set[str], list[dict]]:
    blocks, described = set(), []
    for spec in specs:
        path_text, _, column_name = spec.rpartition("::")
        path = Path(path_text)
        delimiter = "|" if path.read_text(encoding="utf-8", errors="replace")[:4096].count("|") > 3 else "\t"
        count = 0
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                value = first_block(row.get(column_name))
                if value:
                    blocks.add(value)
                    count += 1
        described.append({"file": path.name, "column": column_name, "rows": count})
    return blocks, described


# ---------------------------------------------------------------------------------------- readers

class Chunk:
    """One chunk's outputs, keyed by exported NAME."""

    def __init__(self, folder: Path):
        self.folder = folder
        self.record = json.loads((folder / "msemblator_run.json").read_text(encoding="utf-8")) \
            if (folder / "msemblator_run.json").is_file() else None
        self.name_by_id = {}
        if (folder / "name_map.tsv").is_file():
            self.name_by_id = {row["msemblator_id"]: row["name"] for row in read_tsv(folder / "name_map.tsv")}
            self.structure_step = {row["name"]: row["structure_step"] for row in read_tsv(folder / "name_map.tsv")}
        else:
            self.structure_step = {}

    def rows(self, name: str) -> list[dict]:
        path = self.folder / name
        return read_csv(path) if path.is_file() else []

    def sirius_table(self, relative: str) -> list[dict]:
        path = self.folder / relative
        return read_tsv(path) if path.is_file() else []

    def name_of_feature(self, row: dict) -> str | None:
        """SIRIUS's feature ID for .ms input ends in _<msemblator id>; map it back to the exported NAME."""
        for key in ("mappingFeatureId", "featureId", "id"):
            value = column(row, key)
            if value:
                match = re.search(r"(?:^|_)(\d+)$", value)
                if match and match.group(1) in self.name_by_id:
                    return self.name_by_id[match.group(1)]
        return None


def by_name(rows: list[dict]) -> dict[str, list[dict]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[column(row, "filename")].append(row)
    return grouped


def smiles_to_block(smiles: str, cache: dict) -> str:
    if smiles in cache:
        return cache[smiles]
    from rdkit import Chem, RDLogger  # noqa: PLC0415 (only the evaluation needs RDKit)

    RDLogger.DisableLog("rdApp.*")
    block = ""
    try:
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
        if molecule is not None:
            block = first_block(Chem.MolToInchiKey(molecule))
    except Exception:  # an unparseable SMILES ranks as no match
        block = ""
    cache[smiles] = block
    return block


def chunk_exclusion(record: dict) -> str:
    """Why a chunk's spectra stay out of the main numbers ('' when they count): a non-zero driver exit, or a
    step that failed (for example formula:msbuddy), which leaves the ranking without that tool's evidence.
    Counting such spectra as misses would measure the failure, not the prediction."""
    exit_code = record.get("exit_code")
    steps = record.get("failed_steps") or []
    if exit_code != 0:
        return f"driver exit {exit_code}" + (f"; failed steps {', '.join(steps)}" if steps else "")
    if steps:
        return f"degraded: failed steps {', '.join(steps)}"
    return ""


def min_rank(rows: list[dict], matches) -> float:
    ranks = [number(column(row, "rank")) for row in rows if matches(row)]
    ranks = [rank for rank in ranks if rank is not None]
    return min(ranks) if ranks else math.inf


# ------------------------------------------------------------------------------------- evaluation

def evaluate(run: Path, icd: Path, databases: list[str]) -> tuple[list[dict], dict]:
    plan = json.loads((run / "plan.json").read_text(encoding="utf-8"))
    units = {unit["unit_key"]: unit for unit in plan["units"]}
    full, by_block = load_direct_parents(icd)
    db_blocks, db_described = load_db_blocks(databases)
    smiles_cache: dict[str, str] = {}
    rows: list[dict] = []
    not_run: list[str] = []
    excluded: list[dict] = []
    rankers_seen: dict[str, set] = defaultdict(set)
    canopus_columns: dict[str, list[str]] = {}
    for export in plan["exports"]:
        if export["selection"] != "benchmark":
            continue
        manifest = {row["name"]: row for row in read_tsv(run / export["manifest"])}
        chunk_ids = set(export["chunks"])
        for chunk in plan["chunks"]:
            if chunk["chunk_id"] not in chunk_ids:
                continue
            outputs = Chunk(run / chunk["output_dir"])
            formula_rows = by_name(outputs.rows("formula_score.csv"))
            structure_rows = by_name(outputs.rows("structure_score.csv"))
            sirius_formula = defaultdict(list)
            for row in outputs.sirius_table("work/formula/sirius/output/formula_identifications_top-100.tsv"):
                name = outputs.name_of_feature(row)
                if name:
                    sirius_formula[name].append(row)
            sirius_structure = defaultdict(list)
            for row in outputs.sirius_table("work/structure/sirius_output/structure_identifications_top-100.tsv"):
                name = outputs.name_of_feature(row)
                if name:
                    sirius_structure[name].append(row)
            canopus = {}
            for relative in (outputs.record or {}).get("canopus_summary_files", []):
                table = outputs.sirius_table(f"work/structure/sirius_output/{relative}")
                if table:
                    canopus_columns[relative] = list(table[0].keys())
                if "structure" in relative or not table:
                    continue  # the formula-level summary carries one row per feature's top formula
                for row in table:
                    name = outputs.name_of_feature(row)
                    if name and name not in canopus:
                        canopus[name] = row
            exit_code = (outputs.record or {}).get("exit_code")
            for step, ranker in ((outputs.record or {}).get("rankers") or {}).items():
                rankers_seen[step].add((ranker.get("ranker"), ranker.get("file"), ranker.get("sha256")))
            chunk_names = [line[5:].strip() for line in
                           (run / chunk["chunk"]).read_text(encoding="utf-8").splitlines()
                           if line.upper().startswith("NAME:")]
            if outputs.record is None:
                not_run.extend(chunk_names)
                continue
            exclusion = chunk_exclusion(outputs.record)
            if exclusion:
                excluded.append({"chunk_id": chunk["chunk_id"], "exit_code": exit_code,
                                 "failed_steps": outputs.record.get("failed_steps") or [],
                                 "spectra": len(chunk_names), "reason": exclusion})
            for name in chunk_names:
                spectrum = manifest[name]
                reference_formula = hill(spectrum.get("reference_formula", ""))
                reference_block = first_block(spectrum.get("reference_inchikey"))
                f_rows, s_rows = formula_rows.get(name, []), structure_rows.get(name, [])
                rank_formula = min_rank(f_rows, lambda r: hill(column(r, "formula")) == reference_formula
                                        and reference_formula is not None)
                rank_structure = min_rank(s_rows, lambda r: smiles_to_block(column(r, "SMILES", "Canonical_SMILES"),
                                                                            smiles_cache) == reference_block
                                          and bool(reference_block))
                top = min(s_rows, key=lambda r: number(column(r, "rank")) or math.inf) if s_rows else {}
                sirius_formula_rank = min(
                    (number(column(r, "formulaRank", "rank")) or math.inf for r in sirius_formula.get(name, [])
                     if hill(column(r, "molecularFormula")) == reference_formula and reference_formula),
                    default=math.inf)
                sirius_structure_rank = min(
                    (number(column(r, "structurePerIdRank", "rank")) or math.inf for r in sirius_structure.get(name, [])
                     if first_block(column(r, "InChIkey2D", "InChIKey2D", "inchikey2D")) == reference_block
                     and reference_block),
                    default=math.inf)
                parent = reference_parent(spectrum.get("reference_inchikey", ""), full, by_block)
                canopus_row = canopus.get(name, {})
                most_specific = column(canopus_row, "ClassyFire#most specific class")
                level_terms = []
                for level in ("superclass", "class", "subclass", "level 5", "most specific class"):
                    term = column(canopus_row, f"ClassyFire#{level}")
                    probability = number(column(canopus_row, f"ClassyFire#{level} Probability",
                                                f"ClassyFire#{level} probability"))
                    if term:
                        level_terms.append((term, probability))
                all_classes = [t.strip() for t in column(canopus_row, "ClassyFire#all classifications").split(";")
                               if t.strip()]
                parent_name = (parent[0] if parent else "").casefold()
                unit = units.get(spectrum.get("unit_id"), {})
                precursor_mz = number(spectrum.get("precursor_mz"))
                rows.append({
                    "name": name,
                    "unit_key": spectrum.get("unit_id"),
                    "alignment_id": spectrum.get("alignment_id"),
                    "chunk_id": chunk["chunk_id"],
                    "chunk_exit_code": exit_code,
                    "excluded": exclusion,
                    "formula_evaluable": reference_formula is not None,
                    "structure_evaluable": bool(reference_block),
                    "ion_mode": spectrum.get("ion_mode") or unit.get("ion_mode") or "unknown",
                    "instrument_family": instrument_family(unit.get("instrument")),
                    "precursor_mz": precursor_mz,
                    "mz_bin": mz_bin(precursor_mz),
                    "reference_name": spectrum.get("reference_name", ""),
                    "reference_inchikey": spectrum.get("reference_inchikey", ""),
                    "reference_formula": spectrum.get("reference_formula", ""),
                    "reference_in_db": reference_block in db_blocks if reference_block else False,
                    "formula_returned": bool(f_rows),
                    "structure_step": outputs.structure_step.get(name, ""),
                    "structure_returned": bool(s_rows),
                    "rank_formula": rank_formula,
                    "rank_structure": rank_structure,
                    "top1_structure_score": number(column(top, "rank_score")),
                    "top1_structure_margin": number(column(top, "margin_score")),
                    "top1_structure_probability": number(column(top, "softmax_score", "probability",
                                                                "softmax_probability")),
                    "sirius_formula_rank": sirius_formula_rank,
                    "sirius_structure_rank": sirius_structure_rank,
                    "reference_direct_parent": parent[0] if parent else "",
                    "reference_direct_parent_id": parent[1] if parent else "",
                    "reference_parent_source": parent[2] if parent else "absent",
                    "canopus_returned": bool(canopus_row),
                    "canopus_most_specific": most_specific,
                    "canopus_c1": bool(parent) and bool(most_specific) and most_specific.casefold() == parent_name,
                    **{f"canopus_c2_p{p}": bool(parent) and any(
                        term.casefold() == parent_name and (prob is not None and prob >= p)
                        for term, prob in level_terms) for p in CANOPUS_P},
                    "canopus_in_all_classifications": bool(parent) and parent_name in
                    {t.casefold() for t in all_classes},
                })
    context = {"plan": plan, "databases": db_described, "canopus_columns": canopus_columns,
               "not_run": len(not_run), "excluded_chunks": excluded,
               "rankers": {step: sorted(values, key=str) for step, values in rankers_seen.items()},
               "icd": {"file": icd.name, "sha256": sha256(icd), "rows": len(full)}}
    return rows, context


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rates(rows: list[dict], field: str, ks: tuple[int, ...]) -> dict:
    return {f"hit@{k}": wilson(sum(r[field] <= k for r in rows), len(rows)) for k in ks}


def denominators(rows: list[dict]) -> dict[str, list[dict]]:
    return {
        "N_sent": rows,
        "N_formula_returned": [r for r in rows if r["formula_returned"]],
        "N_structure_returned": [r for r in rows if r["structure_returned"]],
        "N_ref_in_db": [r for r in rows if r["reference_in_db"]],
    }


def score_type(context: dict, step: str) -> str:
    """The step's score identity; 'unknown' or 'mixed' when the records do not name exactly one ranker."""
    seen = {ranker for ranker, _, _ in context.get("rankers", {}).get(step, [])}
    if len(seen) == 1:
        return SCORE_TEMPLATE.format(step=step, ranker=next(iter(seen)))
    return SCORE_TEMPLATE.format(step=step, ranker="mixed" if seen else "unknown")


def summarise(rows: list[dict], context: dict, benchmark_id: str, citation: dict) -> dict:
    # Spectra of failed or degraded chunks stay out of every rate (they are counted and listed instead), and
    # a spectrum whose reference formula or InChIKey cannot be evaluated stays out of that metric's rates.
    included = [r for r in rows if not r["excluded"]]
    structure_rows = [r for r in included if r["structure_evaluable"]]
    formula_rows = [r for r in included if r["formula_evaluable"]]
    structure_groups, formula_groups = denominators(structure_rows), denominators(formula_rows)
    score_formula, score_structure = score_type(context, "formula"), score_type(context, "structure")
    structure = {name: rates(group, "rank_structure", STRUCTURE_K) for name, group in structure_groups.items()}
    formula = {name: rates(group, "rank_formula", FORMULA_K) for name, group in formula_groups.items()}
    baselines = {
        "sirius_formula (SIRIUS formula_identifications_top-100.tsv, formulaRank)":
            {name: rates(group, "sirius_formula_rank", FORMULA_K) for name, group in formula_groups.items()},
        "sirius_csi_fingerid (SIRIUS structure_identifications_top-100.tsv, structurePerIdRank, InChIkey2D)":
            {name: rates(group, "sirius_structure_rank", STRUCTURE_K) for name, group in structure_groups.items()},
        "msfinder": "not evaluated: raw MS-FINDER outputs are not parsed by this version",
        "metfrag": "not evaluated: raw MetFrag outputs are not parsed by this version",
    }
    with_parent = [r for r in included if r["reference_direct_parent"]]
    formula_top1 = [r for r in with_parent if r["rank_formula"] == 1]

    def hit_rate(subset: list[dict], field: str, k: int) -> dict:
        return wilson(sum(r[field] <= k for r in subset), len(subset))

    def canopus_block(subset: list[dict]) -> dict:
        block = {"C1_most_specific_equals_direct_parent": wilson(sum(r["canopus_c1"] for r in subset), len(subset))}
        for p in CANOPUS_P:
            block[f"C2_direct_parent_at_posterior_ge_{p}"] = wilson(sum(r[f"canopus_c2_p{p}"] for r in subset),
                                                                  len(subset))
        block["direct_parent_in_all_classifications_no_posterior"] = wilson(
            sum(r["canopus_in_all_classifications"] for r in subset), len(subset))
        return block

    strata = {}
    for key in ("unit_key", "ion_mode", "instrument_family", "mz_bin"):
        values = sorted({r[key] for r in included if r[key] is not None}, key=str)
        strata[key] = {str(v): {"n": sum(r[key] == v for r in included),
                                "structure_hit@1": hit_rate([r for r in structure_rows if r[key] == v],
                                                            "rank_structure", 1),
                                "structure_hit@5": hit_rate([r for r in structure_rows if r[key] == v],
                                                            "rank_structure", 5),
                                "formula_hit@1": hit_rate([r for r in formula_rows if r[key] == v], "rank_formula", 1)}
                         for v in values}

    def unit_average(metric: str) -> float | None:
        values = [s[metric]["rate"] for s in strata["unit_key"].values() if s[metric]["rate"] is not None]
        return round(sum(values) / len(values), 4) if values else None

    parents = defaultdict(list)
    for row in with_parent:
        parents[row["reference_direct_parent"]].append(row)
    strata["reference_direct_parent_n_ge_5"] = {
        name: {"n": len(group),
               "structure_hit@1": hit_rate([r for r in group if r["structure_evaluable"]], "rank_structure", 1),
               "canopus_C1": wilson(sum(r["canopus_c1"] for r in group), len(group))}
        for name, group in sorted(parents.items()) if len(group) >= 5}
    margins = sorted(r["top1_structure_margin"] for r in structure_rows if r["top1_structure_margin"] is not None)
    closure = []
    for quantile in (0.0, 0.25, 0.5, 0.75, 0.9):
        if not margins:
            break
        cut = margins[min(len(margins) - 1, int(quantile * len(margins)))]
        subset = [r for r in structure_rows
                  if r["top1_structure_margin"] is not None and r["top1_structure_margin"] >= cut]
        closure.append({"top1_margin_ge": round(cut, 4), "n": len(subset),
                        **{f"hit@{k}": wilson(sum(r["rank_structure"] <= k for r in subset), len(subset))
                           for k in (1, 5)}})
    return {
        "benchmark_id": benchmark_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "truth": "agreement with the strict reference-library match (MS-DIAL Level-2), not accuracy",
        "reference_library_citation": citation,
        "run_id": context["plan"]["run_id"],
        "plan_sha256": context.get("plan_sha256"),
        "rdkit_version": context.get("rdkit_version"),
        "icd": context["icd"],
        "candidate_databases": context["databases"],
        "candidate_database_note": "N_ref_in_db checks the MS-FINDER and MetFrag text databases only; coverage of "
                                   "the SIRIUS custom database is not checked",
        "canopus_summary_columns": context["canopus_columns"],
        "counts": {
            "benchmark_spectra_evaluated": len(included),
            **{f"structure_{name}": len(group) for name, group in structure_groups.items()},
            **{f"formula_{name}": len(group) for name, group in formula_groups.items()},
            "structure_reference_unevaluable": sum(not r["structure_evaluable"] for r in included),
            "formula_reference_unevaluable": sum(not r["formula_evaluable"] for r in included),
            "benchmark_spectra_in_excluded_chunks": len(rows) - len(included),
            "excluded_chunks": len(context.get("excluded_chunks", [])),
            "benchmark_spectra_in_chunks_not_run": context.get("not_run", 0),
            "with_reference_direct_parent": len(with_parent),
            "icd_coverage": wilson(len(with_parent), len(included)),
            "skipped_no_formula": sum(r["structure_step"] == "skipped_no_formula" for r in included)},
        "excluded_chunks": context.get("excluded_chunks", []),
        "exclusion_rule": "a chunk whose driver exit is not 0, or that reports a failed step, is left out of every "
                          "rate; a spectrum whose reference formula (charged or unparseable) or InChIKey (empty) "
                          "cannot be evaluated is left out of that metric's rates",
        "rankers": {step: [{"ranker": r, "file": f, "sha256": h} for r, f, h in values]
                    for step, values in context.get("rankers", {}).items()},
        "structure": {"tool": TOOL_STRUCTURE, "score_type": score_structure,
                      "posterior_type": POSTERIOR_MSEMBLATOR, "match": "first InChIKey block (RDKit from SMILES)",
                      "denominators": structure},
        "formula": {"tool": TOOL_FORMULA, "score_type": score_formula, "posterior_type": POSTERIOR_MSEMBLATOR,
                    "match": "Hill-normalised neutral formula, minimum rank across adduct groups",
                    "denominators": formula},
        "baselines": baselines,
        "canopus": {"tool": TOOL_CANOPUS, "posterior_type": POSTERIOR_CANOPUS,
                    "reference": "ClassyFire direct parent from the ICD (full InChIKey, else a unique first block)",
                    "all": canopus_block(with_parent), "formula_top1": canopus_block(formula_top1),
                    "note": "CANOPUS runs on MS-Emblator's top-1 formula, so the formula_top1 block is the "
                            "informative one; class/superclass agreement needs the ChemOnt lineage (Q5)"},
        "unit_macro_average": {
            "units": sum(s["n"] > 0 for s in strata["unit_key"].values()),
            "structure_hit@1": unit_average("structure_hit@1"),
            "structure_hit@5": unit_average("structure_hit@5"),
            "formula_hit@1": unit_average("formula_hit@1")},
        "strata": strata,
        "closure_by_top1_margin": {"tool": TOOL_STRUCTURE, "score_type": score_structure,
                                   "margin": "rank_score(top1) - rank_score(top2)", "rows": closure},
        "caveats": [
            "Q3: CANOPUS is conditional on MS-Emblator's top-1 formula; class agreement with a wrong formula is "
            "uninformative.",
            "Q4: the truth is a Level-2 library match, not an authentic standard, and CSI:FingerID/CANOPUS "
            "training data overlap public libraries contained in the reference library; read the numbers as "
            "optimistic agreement.",
        ],
    }


def markdown(summary: dict, cited_rows: list[dict] | None) -> str:
    citation = summary["reference_library_citation"]
    lines = [f"# MS-Emblator benchmark {summary['benchmark_id']}", "",
             f"Truth: {summary['truth']}. Reference library: {citation['text']}.", "",
             "| Metric | Tool / score | Denominator | n | rate | 95% CI |", "|---|---|---|---|---|---|"]
    for kind in ("structure", "formula"):
        block = summary[kind]
        for denominator, values in block["denominators"].items():
            for metric, value in values.items():
                if value["n"]:
                    lines.append(f"| {kind} {metric} | {block['tool']} / {block['score_type']} | {denominator} | "
                                 f"{value['n']} | {value['rate']} | {value['low']}-{value['high']} |")
    for subset in ("all", "formula_top1"):
        for metric, value in summary["canopus"][subset].items():
            if value["n"]:
                lines.append(f"| CANOPUS {metric} ({subset}) | {summary['canopus']['tool']} / "
                             f"{summary['canopus']['posterior_type']} | with direct parent | {value['n']} | "
                             f"{value['rate']} | {value['low']}-{value['high']} |")
    lines += ["", "Counts: " + ", ".join(f"{k} {v}" for k, v in summary["counts"].items() if not isinstance(v, dict)),
              "", f"Exclusions: {summary['exclusion_rule']}."]
    if summary["excluded_chunks"]:
        lines += ["", "| Excluded chunk | spectra | reason |", "|---|---|---|"]
        lines += [f"| {c['chunk_id']} | {c['spectra']} | {c['reason']} |" for c in summary["excluded_chunks"]]
    lines += ["", "Caveats:"] + [f"- {c}" for c in summary["caveats"]]
    if cited_rows:
        lines += ["", f"Per-spectrum rows (reference identities from {citation['text']}):", "",
                  "| spectrum | reference | first block | rank formula | rank structure | CANOPUS most specific | "
                  "reference direct parent |", "|---|---|---|---|---|---|---|"]
        for row in cited_rows:
            lines.append(f"| {row['name']} | {row['reference_name']} | {first_block(row['reference_inchikey'])} | "
                         f"{row['rank_formula']} | {row['rank_structure']} | {row['canopus_most_specific']} | "
                         f"{row['reference_direct_parent']} |")
    return "\n".join(lines) + "\n"


def library_citation(run: Path) -> dict:
    """The strict reference libraries as run, by file name and sha256, from the units' campaign records."""
    plan = json.loads((run / "plan.json").read_text(encoding="utf-8"))
    libraries = {}
    for unit in plan["units"]:
        output = Path(unit["output_dir"])
        for record_path in (output.parent / "campaign-record.json", output / "campaign-record.json"):
            if not record_path.is_file():
                continue
            record = json.loads(record_path.read_text(encoding="utf-8"))
            for library in record.get("libraries_as_run", []):
                name = Path(str(library.get("name", ""))).name
                if name.lower().endswith(".msp") and "private" in name.lower():
                    libraries[name] = library.get("sha256")
    text = "; ".join(f"{name} (sha256 {digest})" for name, digest in sorted(libraries.items())) \
        or "reference library not recorded"
    return {"libraries": [{"file": n, "sha256": d} for n, d in sorted(libraries.items())], "text": text}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", required=True, type=Path, help="The batch run folder (holds plan.json)")
    parser.add_argument("--icd", required=True, type=Path)
    parser.add_argument("--db", action="append", default=[], metavar="PATH::COLUMN")
    parser.add_argument("--out", type=Path, help="Default: <run>/benchmark")
    parser.add_argument("--benchmark-id", required=True)
    parser.add_argument("--per-spectrum-rows", action="store_true",
                        help="Add per-spectrum rows with reference identities to the shareable report")
    args = parser.parse_args(argv)
    out = args.out or args.run / "benchmark"
    out.mkdir(parents=True, exist_ok=True)
    rows, context = evaluate(args.run, args.icd, args.db)
    context["plan_sha256"] = sha256(args.run / "plan.json")
    try:
        import rdkit  # noqa: PLC0415

        context["rdkit_version"] = rdkit.__version__
    except ImportError:
        context["rdkit_version"] = None
    citation = library_citation(args.run)
    fields = list(rows[0].keys()) if rows else ["name"]
    with (out / "benchmark_rows.local-vs21.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = summarise(rows, context, args.benchmark_id, citation)
    summary_json = out / "benchmark_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    summary_md = out / "benchmark_summary.md"
    summary_md.write_text(markdown(summary, rows if args.per_spectrum_rows else None), encoding="utf-8")
    identities = sorted({value for row in rows for value in (row["reference_name"], row["reference_inchikey"])
                         if value and len(value) > 3})
    marker = citation["text"] if citation["libraries"] else None  # the fallback text is not a citation
    hits, uncited, detail = self_check([summary_json, summary_md], identities, marker)
    print(f"benchmark rows: {len(rows)} ({summary['counts']['benchmark_spectra_evaluated']} evaluated, "
          f"{summary['counts']['benchmark_spectra_in_excluded_chunks']} in excluded chunks)")
    print(f"summary sha256: {sha256(summary_json)}")
    print(f"leak hits: {hits}" + ("" if not hits else f" {detail}"))
    print("citation: ok" if not uncited else f"citation: missing in {', '.join(uncited)}")
    return 0 if not hits and not uncited else 1


if __name__ == "__main__":
    sys.exit(main())
