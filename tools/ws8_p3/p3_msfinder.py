"""WS8 P3 MS-FINDER-level measurement of the author-reviewed diagnostic-ion rules: option off versus on.

Stages: select (representatives of the frozen rule-level run on which a reviewed rule fired or was veto-blocked in
the base pass and whose truth first block is in MS-FINDER's structure DB, plus a seeded control without a reviewed
outcome), predict (MS-FINDER console `predict`, rules off and on at each weight, BelowNormal, at most 8 workers),
measure (structure rank of the truth, top-1/top-5, boost-wrong rate and ChemOnt class-level accuracy, per family and
instrument class, record and cluster level) and all. A VS20 run picks the weight; a VS21 run (--vs21) uses only the
weight a completed VS20 selection run recorded. Optional chemistry (pyarrow) stays outside src/.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import queue
import random
import re
import shutil
import statistics
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from tools.ws8_p2 import p2_diagnostics as p2  # noqa: E402  (also pins the thread variables to 1)
from tools.ws8_p3 import p3_measure as p3m  # noqa: E402
from tools.ws8_p3 import rule_format as rf  # noqa: E402

OUT_ROOT = p3m.OUT_ROOT
RUN_NAME = r"p3ms(smoke)?-(\d{8})-([a-f0-9]{7,40})"
# The frozen rule-level run whose base pass decides the selection (prepared.*, fires.local-vs21.parquet).
FROZEN_RUN = OUT_ROOT / "p3cs-20261124-5a8fc37a"
EXPECTED_RULESET_VERSION = "0.2.0"
EXPECTED_RULESET_SHA256 = "ecf7988dea90992eb7233fc856ce2cb19a8afe4f1f5683628a4b35a459c93569"
WEIGHTS = (0.5, 1.0, 2.0, 4.0)
SMOKE_WEIGHTS = (1.0,)
CONTROL_N = 200
SEED = p3m.SEED
VS21_MIN_CELL = 5
CORES_PER_WORKER = 2  # 8 workers x 2 logical processors = the 16-thread ceiling
SUCCESS_EXIT = 1  # PredictProcess.Run returns 1 on success and -1 on error
CODE_FILES = (HERE / "p3_msfinder.py", HERE / "p3_measure.py", HERE / "rule_format.py",
              ROOT / "tools/ws8_p2/p2_diagnostics.py")
# Resources the console's structure search needs beside it (FileStorageUtility.GetResourcesPath); the console
# writes a MessagePack cache "<stem>_bin<ext>" beside each on first use, so the first job runs alone.
RESOURCES = ("MsfinderStructureDB-VS15.esd", "InchikeyClassyfireDB-VS5.icd", "MINEs-StructureDB-vs1.msd",
             "MsfinderFormulaDB-VS13.efd")
BIN_CACHES = ("MsfinderStructureDB-VS15_bin.esd", "MINEs-StructureDB-vs1_bin.msd", "MsfinderFormulaDB-VS13_bin.efd")
# MS-FINDER predict method: defaults except the tolerances (Da, as the rule-level console INI), the local structure DB
# only (no PubChem request, no MINEs), and the diagnostic-rule block appended per variant.
METHOD_LINES = ("MassToleranceType=Da", "Ms1Tolerance=0.01", "Ms2Tolerance=0.01", "IsTmsMeoxDerivative=False",
                "TryTopNMolecularFormulaSearch=5", "StructureMaximumReportNumber=100", "IsUseEiFragmentDB=False",
                "MinesNeverUse=True", "MinesOnlyUseForNecessary=False", "MinesAllways=False",
                "PubChemNeverUse=True", "PubChemOnlyUseForNecessary=False", "PubChemAllways=False")
RANKING_SCOPE = "Reviewed"
CLASS_DEPTH = 3  # ChemOnt: Chemical entities (0) > kingdom (1) > superclass (2) > class (3)
FILE_RE = re.compile(r"rec(\d+)(?:\.local-vs21)?")
ALL = "all"


def suffix(source):
    return p3m.suffix(source)


def variant_name(weight):
    return "off" if weight is None else "w" + p3m.num(weight)


def method_ini(weight):
    lines = list(METHOD_LINES)
    if weight is not None:
        lines += ["DiagnosticRuleEvaluation=TRUE", "DiagnosticRuleRankingWeight=" + p3m.num(weight),
                  "DiagnosticRuleRankingScope=" + RANKING_SCOPE]
    return "\r\n".join(lines) + "\r\n"


def console_class(instrument_type):
    """The console's instrument class for an MSP INSTRUMENTTYPE (blank is unknown; else P2's bins)."""
    return "unknown" if not (instrument_type or "").strip() else p2.instrument_bin(instrument_type)


def family_of(rule_id):
    return rule_id.split(".")[0]


# ------------------------------------------------------------------------------------------- libraries
def read_short_keys(esd):
    """First blocks of MS-FINDER's structure DB (the text .esd, column "Short InChIKey")."""
    with open(esd, encoding="utf-8-sig", errors="replace", newline="") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        col = header.index("Short InChIKey")
        return {parts[col] for parts in (line.rstrip("\r\n").split("\t") for line in f) if len(parts) > col and parts[col]}


def read_icd(icd):
    """MS-FINDER's InChIKey -> direct parent lookup, first row per first block (as SetClassyfireOntologies)."""
    by_block, name_to_id = {}, {}
    with open(icd, encoding="utf-8-sig", errors="replace", newline="") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3 or not parts[0]:
                continue
            block, name, cid = parts[0].split("-")[0], parts[1], parts[2]
            by_block.setdefault(block, (name, cid))
            if re.fullmatch(r"CHEMONTID:\d{7}", cid):
                name_to_id.setdefault(name, cid)
    return by_block, name_to_id


def read_obo(path):
    """ChemOnt terms: id -> first is_a parent, and name -> id."""
    parents, names, cur = {}, {}, None
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {}
            elif cur is not None and line.startswith("id: "):
                cur["id"] = line[4:].strip()
            elif cur is not None and line.startswith("name: "):
                names.setdefault(line[6:].strip(), cur.get("id"))
            elif cur is not None and line.startswith("is_a: ") and cur.get("id") and cur["id"] not in parents:
                parents[cur["id"]] = line[6:].split("!")[0].strip()
    return parents, names


def class_level(term, parents, depth=CLASS_DEPTH):
    """The ancestor of a ChemOnt term at the class level, or None above it or when unknown."""
    if not term:
        return None
    path, seen = [term], {term}
    while path[-1] in parents:
        nxt = parents[path[-1]]
        if nxt in seen:
            return None
        path.append(nxt)
        seen.add(nxt)
    path.reverse()  # root first
    return path[depth] if len(path) > depth else None


# ------------------------------------------------------------------------------------------- selection
def reviewed_rules(ruleset_dir):
    res = rf.load_files(rf.read_directory(ruleset_dir))
    if not res.ok:
        raise RuntimeError("Rule set does not load")
    if res.ruleset_sha256 != EXPECTED_RULESET_SHA256 or res.rule_set_version != EXPECTED_RULESET_VERSION:
        raise RuntimeError("Rule set is not the pinned 0.2.0 set")
    reviewed = {rid: dict(family=family_of(rid), motif_id=info.doc["claim"]["motif_id"], version=info.doc["version"])
                for rid, info in sorted(res.rules.items()) if info.doc.get("author_reviewed") is True}
    return reviewed, res


def select_records(reps, outcomes, any_fired, db_keys, reviewed, seed, control_n=CONTROL_N, limit=None):
    """reps: [{record_id, key, mode, cluster_id, instrument_type}]; outcomes: record_id -> [(rule_id, fired, veto)]
    of the base pass. Returns (rows, counts); rows carry no truth key."""
    counts = Counter(representatives=len(reps))
    selected, pool = [], []
    for r in sorted(reps, key=lambda x: x["record_id"]):
        rid = r["record_id"]
        hits = sorted((rule, "fired" if fired else "veto_blocked") for rule, fired, veto in outcomes.get(rid, ())
                      if rule in reviewed and (fired or veto))
        in_db = r["key"] in db_keys
        counts["in_db"] += in_db
        counts["reviewed_outcome"] += bool(hits)
        if not in_db:
            counts["reviewed_outcome_not_in_db"] += bool(hits)
            continue
        base = dict(record_id=rid, mode=r["mode"], cluster_id=r["cluster_id"],
                    instrument_class=console_class(r.get("instrument_type")))
        if hits:
            if len(hits) > 1:
                counts["selected_multi_rule"] += 1
            selected.append(dict(base, set="selected", rule_id=hits[0][0], family=reviewed[hits[0][0]]["family"],
                                 outcome=hits[0][1]))
        else:
            counts["in_db_no_reviewed_outcome"] += 1
            counts["in_db_no_rule_fired"] += rid not in any_fired
            pool.append(dict(base, set="control", rule_id="", family="control", outcome="none"))
    rng = random.Random(seed)
    if limit is not None:
        selected = sorted(rng.sample(selected, min(limit, len(selected))), key=lambda x: x["record_id"])
        control_n = min(control_n, max(1, limit // 4))
    control = sorted(rng.sample(pool, min(control_n, len(pool))), key=lambda x: x["record_id"])
    counts.update(selected=len(selected), control=len(control), control_pool=len(pool))
    for r in selected:
        counts[f"selected:{r['family']}:{r['outcome']}"] += 1
        counts[f"selected_instrument:{r['instrument_class']}"] += 1
    return selected + control, dict(sorted(counts.items()))


def write_query_msp(path, r):
    """One query per file, as MS-FINDER predict reads a folder: NAME carries only the record id."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"NAME: rec{r['record_id']}\nPRECURSORMZ: {p3m.num(r['precursor'])}\nPRECURSORTYPE: {r['adduct']}\n")
        f.write("IONMODE: " + ("Positive" if r["mode"] == "pos" else "Negative") + "\n")
        if r.get("instrument_type"):
            f.write(f"INSTRUMENTTYPE: {r['instrument_type']}\n")
        if r.get("collision_energy"):
            f.write(f"COLLISIONENERGY: {r['collision_energy']}\n")
        f.write(f"SPECTRUMTYPE: Centroid\nNum Peaks: {len(r['mz'])}\n")
        for mz, intensity in zip(r["mz"], r["intensity"]):
            f.write(f"{p3m.num(mz)}\t{p3m.num(intensity)}\n")


def load_frozen(frozen, source, reviewed):
    """Representatives of one source and their base-pass outcomes (the code reads VS21 files; nothing is printed)."""
    import pyarrow.parquet as pq
    cols = ["record_id", "key", "mode", "cluster_id", "instrument_type"]
    reps = pq.read_table(p2.row_path(Path(frozen), "prepared", source), columns=cols).to_pylist()
    fires = pq.read_table(Path(frozen) / "fires.local-vs21.parquet",
                          columns=["record_id", "rule_id", "fired", "veto_blocked"],
                          filters=[("source", "=", source), ("set", "=", "representatives"), ("pass_id", "=", "base")])
    outcomes, any_fired = defaultdict(list), set()
    for r in fires.to_pylist():
        if r["fired"]:
            any_fired.add(r["record_id"])
        if r["rule_id"] in reviewed:
            outcomes[r["record_id"]].append((r["rule_id"], r["fired"], r["veto_blocked"]))
    return reps, outcomes, any_fired


def select(args, run):
    reviewed, _ = reviewed_rules(args.ruleset_dir)
    res_dir = Path(args.console).parent / "Resources"
    db_keys = read_short_keys(res_dir / RESOURCES[0])
    reps, outcomes, any_fired = load_frozen(args.frozen_run, args.source, reviewed)
    rows, counts = select_records(reps, outcomes, any_fired, db_keys, reviewed, args.seed, limit=args.limit)
    counts["structure_db_first_blocks"] = len(db_keys)
    sfx = suffix(args.source)
    p3m.write_rows(args.out / f"selection{sfx}.tsv", rows,
                   ["record_id", "set", "family", "rule_id", "outcome", "mode", "instrument_class", "cluster_id"])
    import pyarrow.parquet as pq
    ids = [r["record_id"] for r in rows]
    spectra = pq.read_table(p2.row_path(Path(args.frozen_run), "prepared", args.source),
                            columns=["record_id", "precursor", "adduct", "mode", "instrument_type", "collision_energy",
                                     "mz", "intensity"], filters=[("record_id", "in", ids)]).to_pylist()
    msp_dir = args.out / "msp" / (args.source + (".local-vs21.d" if sfx else ""))
    msp_dir.mkdir(parents=True, exist_ok=True)
    for r in spectra:
        write_query_msp(msp_dir / f"rec{r['record_id']}{sfx}.msp", r)
    if len(spectra) != len(rows):
        raise RuntimeError("Selected records missing from the prepared table")
    run["selection"] = dict(source=args.source, reviewed_rules=reviewed, counts=counts)
    print("select", args.source, {k: v for k, v in counts.items()}, flush=True)


# ------------------------------------------------------------------------------------------- predict
def shards_of(ids, n):
    blocks = p3m.shard_blocks(len(ids), n)
    return [ids[lo:hi] for lo, hi in blocks if hi > lo]


def set_affinity(proc, slot, cores):
    if sys.platform != "win32":
        return None
    import ctypes
    mask = ((1 << cores) - 1) << (slot * cores)
    if not ctypes.windll.kernel32.SetProcessAffinityMask(int(proc._handle), mask):
        raise RuntimeError("SetProcessAffinityMask failed")
    return mask


def validate_job(out_dir, weight, n_records):
    """One structure result; for a rules-on job the diagnostic manifest names the pinned 0.2.0 set, the weight, the
    scope and no errors or refusals."""
    results = sorted(Path(out_dir).glob("Structure result-*.txt"))
    if len(results) != 1:
        raise RuntimeError("Expected one structure result in " + Path(out_dir).name)
    info = dict(result=results[0].name, result_sha256=p2.sha256(results[0]))
    manifests = sorted(Path(out_dir).glob("Diagnostic rule run-*.json"))
    if weight is None:
        if manifests:
            raise RuntimeError("Rules-off job wrote a diagnostic manifest")
        return info
    if len(manifests) != 1:
        raise RuntimeError("Expected one diagnostic manifest")
    m = json.loads(manifests[0].read_text(encoding="utf-8"))
    rs, st, c = m["rule_set"], m["settings"], m["counts"]
    if rs["version"] != EXPECTED_RULESET_VERSION or rs["sha256"] != EXPECTED_RULESET_SHA256 or not rs["pinned"]:
        raise RuntimeError("Console rule set is not the pinned 0.2.0 set")
    if not st["evaluation"] or float(st["ranking_weight"]) != float(weight) or st["ranking_scope"].lower() != RANKING_SCOPE.lower():
        raise RuntimeError("Console settings disagree with the variant")
    if c["errors"] or c["refused"] or c["queries"] != n_records:
        raise RuntimeError("Console reported errors, refusals or a query count mismatch")
    return dict(info, manifest=manifests[0].name, manifest_sha256=p2.sha256(manifests[0]),
                build_identity=m["build_identity"], engine_semantics=m["engine_semantics"],
                counts={k: c[k] for k in ("queries", "evaluated", "fired", "veto_blocked", "none", "skipped")},
                timing=m.get("timing"))


def run_job(console, job, slots, cores):
    """One console process per (variant, shard) at BelowNormal on its own cores. A job folder is never removed: a
    folder without job.json is a failed or running job and is refused."""
    d = job["dir"]
    done = d / "job.json"
    if done.exists():
        info = json.loads(done.read_text(encoding="utf-8"))
        validate_job(d / "out", job["weight"], len(job["ids"]))
        return dict(info, reused=True)
    if d.exists():
        raise RuntimeError(f"Partial job folder {d.name} exists; start a new run id")
    (d / "input").mkdir(parents=True)
    for p in job["msps"]:
        shutil.copyfile(p, d / "input" / p.name)
    (d / "method.ini").write_text(method_ini(job["weight"]), encoding="ascii", newline="")
    slot = slots.get()
    try:
        start = time.perf_counter()
        with open(d / "console.log", "w", encoding="utf-8") as log:
            proc = subprocess.Popen([str(console), "predict", "-i", str(d / "input"), "-m", str(d / "method.ini"),
                                     "-o", str(d / "out")], stdout=log, stderr=subprocess.STDOUT,
                                    creationflags=p3m.BELOW_NORMAL | p3m.NO_WINDOW if sys.platform == "win32" else 0)
            set_affinity(proc, slot, cores)
            code = proc.wait()
        wall = time.perf_counter() - start
    finally:
        slots.put(slot)
    if code != SUCCESS_EXIT:
        raise RuntimeError(f"predict exited {code} in {d.name}")
    info = dict(validate_job(d / "out", job["weight"], len(job["ids"])), variant=job["variant"], shard=job["shard"],
                n=len(job["ids"]), wall_seconds=round(wall, 3), exit_code=code)
    p2.write_json(done, info)
    return info


def plan_jobs(args, weights, ids):
    sfx = suffix(args.source)
    msp_dir = args.out / "msp" / (args.source + (".local-vs21.d" if sfx else ""))
    jobs = []
    for w in (None,) + tuple(weights):
        for s, block in enumerate(shards_of(ids, args.shards)):
            jobs.append(dict(weight=w, variant=variant_name(w), shard=s, ids=block,
                             dir=args.out / "predict" / variant_name(w) / f"shard-{s:02d}{'.local-vs21.d' if sfx else ''}",
                             msps=[msp_dir / f"rec{i}{sfx}.msp" for i in block]))
    return jobs


def predict(args, run):
    sel = p3m.read_rows(args.out / f"selection{suffix(args.source)}.tsv")
    ids = [int(r["record_id"]) for r in sel]
    jobs = plan_jobs(args, args.weights, ids)
    slots = queue.Queue()
    for i in range(args.workers):
        slots.put(i)
    res_dir = Path(args.console).parent / "Resources"
    results = {}
    first = [] if all((res_dir / b).is_file() for b in BIN_CACHES) else jobs[:1]
    for job in first:  # the console's MessagePack caches are written on first use: never by two processes at once
        results[(job["variant"], job["shard"])] = run_job(args.console, job, slots, args.cores_per_worker)
    with ThreadPoolExecutor(args.workers) as pool:
        futures = {(j["variant"], j["shard"]): pool.submit(run_job, args.console, j, slots, args.cores_per_worker)
                   for j in jobs if (j["variant"], j["shard"]) not in results}
        for k, f in futures.items():
            results[k] = f.result()
    identities = {(r["build_identity"], r["engine_semantics"]) for r in results.values() if "build_identity" in r}
    if len(identities) > 1:
        raise RuntimeError("Jobs disagree on the console build identity")
    walls = [r["wall_seconds"] for r in results.values() if "wall_seconds" in r]
    n_spectra = sum(r["n"] for r in results.values() if "wall_seconds" in r)
    run["predict"] = dict(variants=sorted({j["variant"] for j in jobs}), shards=len({j["shard"] for j in jobs}),
                          jobs={f"{v}/{s}": r for (v, s), r in sorted(results.items())},
                          build_identity=next(iter(identities))[0] if identities else None,
                          caches_written_by_this_run=bool(first),
                          process_seconds_per_spectrum=round(sum(walls) / n_spectra, 3) if n_spectra else None,
                          bin_caches={b: p2.sha256(res_dir / b) for b in BIN_CACHES if (res_dir / b).is_file()})
    print("predict", run["predict"]["variants"], "jobs", len(results), "s/spectrum",
          run["predict"]["process_seconds_per_spectrum"], flush=True)


# ------------------------------------------------------------------------------------------- measure
def parse_structure_result(path):
    """record_id -> candidates in the written order: (first block, ontology, rule term, base rank)."""
    out = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            m = FILE_RE.fullmatch(row["File name"])
            if not m:
                raise RuntimeError("Unexpected query file name in the structure result")
            out[int(m.group(1))].append(dict(block=(row.get("InChIKey") or "").split("-")[0], ontology=row.get("Ontology") or "",
                                             term=float(row.get("Diagnostic rule term") or 0),
                                             base_rank=int(row["Base rank"]) if row.get("Base rank") else None))
    return dict(out)


def parse_diagnostic_result(path):
    """record_id -> the predict run's own outcomes of the reviewed rules (claim eligibility, instrument class)."""
    out = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            m = FILE_RE.fullmatch(Path(row["File name"]).stem)
            if m:
                out[int(m.group(1))].append(row)
    return dict(out)


def rank_of(cands, truth):
    for i, c in enumerate(cands, 1):
        if c["block"] == truth:
            return i
    return None


def boost_wrong(off, on, truth):
    """The rules raised a wrong candidate above the true structure (or above the off-run top-1 that was correct):
    some wrong candidate is below the truth in the off run and above it in the on run."""
    r_off, r_on = rank_of(off, truth), rank_of(on, truth)
    if r_off is None or r_on is None:
        return False
    below = {c["block"] for c in off[r_off:]} - {truth}
    above = {c["block"] for c in on[:r_on - 1]} - {truth}
    return bool(below & above)


def spectrum_row(rec, truth, truth_class, cands, off_cands, parents, name_to_id, onto_names):
    rank = rank_of(cands, truth)
    top = cands[0] if cands else None
    top_id = (name_to_id.get(top["ontology"]) or onto_names.get(top["ontology"])) if top else None
    top_class = class_level(top_id, parents)
    row = dict(rec, rank=rank, found=rank is not None, top1=rank == 1, top5=rank is not None and rank <= 5,
               n_candidates=len(cands), class_defined=truth_class is not None,
               class_correct=truth_class is not None and top_class == truth_class)
    if off_cands is not None:
        row["boost_wrong"] = boost_wrong(off_cands, cands, truth)
        row["order_changed"] = [c["block"] for c in cands] != [c["block"] for c in off_cands]
        row["rank_changed"] = rank != rank_of(off_cands, truth)
        # The on run's Base rank must reproduce the off run's order (deterministic search).
        base = sorted(cands, key=lambda c: (c["base_rank"] is None, c["base_rank"] or 0))
        row["base_rank_mismatch"] = [c["block"] for c in base] != [c["block"] for c in off_cands]
    return row


RATE_FIELDS = ("top1", "top5", "found", "class_correct", "boost_wrong", "rank_changed")


def cell(rows, on):
    n = len(rows)
    found = [r["rank"] for r in rows if r["rank"] is not None]
    c = dict(n=n, top1=sum(r["top1"] for r in rows), top5=sum(r["top5"] for r in rows), found=len(found),
             not_found=n - len(found), rank_le10=sum(1 for x in found if x <= 10),
             median_rank_found=statistics.median(found) if found else None,
             class_n=sum(r["class_defined"] for r in rows),
             class_correct=sum(r["class_correct"] for r in rows if r["class_defined"]))
    if on:
        c.update(boost_wrong=sum(r["boost_wrong"] for r in rows), rank_changed=sum(r["rank_changed"] for r in rows),
                 order_changed=sum(r["order_changed"] for r in rows))
    c.update(top1_rate=_rate(c["top1"], n), top5_rate=_rate(c["top5"], n), not_found_rate=_rate(c["not_found"], n),
             class_rate=_rate(c["class_correct"], c["class_n"]),
             boost_wrong_rate=_rate(c["boost_wrong"], n) if on else None)
    # Cluster level: every cluster weighs one (the mean over clusters of the within-cluster rate).
    by_cluster = defaultdict(list)
    for r in rows:
        by_cluster[r["cluster_id"]].append(r)
    c["n_clusters"] = len(by_cluster)

    def cmean(field, cond=lambda r: True):
        vals = [sum(r[field] for r in g if cond(r)) / k for g in by_cluster.values()
                for k in [sum(1 for r in g if cond(r))] if k]
        return (round(sum(vals) / len(vals), 6), len(vals)) if vals else (None, 0)
    c["c_top1_rate"], _ = cmean("top1")
    c["c_top5_rate"], _ = cmean("top5")
    nf, _ = cmean("found")
    c["c_not_found_rate"] = None if nf is None else round(1 - nf, 6)
    c["c_class_rate"], c["c_class_n"] = cmean("class_correct", lambda r: r["class_defined"])
    c["c_boost_wrong_rate"] = cmean("boost_wrong")[0] if on else None
    return c


def _rate(k, n):
    return round(k / n, 6) if n else None


CELL_FIELDS = ["source", "variant", "weight", "set", "family", "instrument_class", "n", "n_clusters", "top1", "top1_rate",
               "top5", "top5_rate", "found", "not_found", "not_found_rate", "rank_le10", "median_rank_found", "class_n",
               "class_correct", "class_rate", "boost_wrong", "boost_wrong_rate", "rank_changed", "order_changed",
               "c_top1_rate", "c_top5_rate", "c_not_found_rate", "c_class_n", "c_class_rate", "c_boost_wrong_rate",
               "suppressed"]
SUPPRESS_BY = {"n": ("top1", "top1_rate", "top5", "top5_rate", "found", "not_found", "not_found_rate", "rank_le10",
                     "median_rank_found", "boost_wrong", "boost_wrong_rate", "rank_changed", "order_changed"),
               "class_n": ("class_correct", "class_rate"),
               "n_clusters": ("c_top1_rate", "c_top5_rate", "c_not_found_rate", "c_boost_wrong_rate"),
               "c_class_n": ("c_class_rate",)}


def suppress_cell(c, source, n_min=VS21_MIN_CELL):
    """VS21 only: every value whose denominator is below n_min is blanked (the denominators themselves stay)."""
    c = dict(c)
    blanked = []
    if source == "vs21":
        for den, fields in SUPPRESS_BY.items():
            if (c.get(den) or 0) < n_min:
                for f in fields:
                    if c.get(f) is not None:
                        c[f] = None
                        blanked.append(f)
    c["suppressed"] = ";".join(blanked)
    return c


def aggregate(rows_by_variant, source, weights):
    cells = []
    for w in (None,) + tuple(weights):
        v = variant_name(w)
        rows = rows_by_variant[v]
        for s in ("selected", "control"):
            srows = [r for r in rows if r["set"] == s]
            for fam in [ALL] + sorted({r["family"] for r in srows}):
                for ic in [ALL] + sorted({r["instrument_class"] for r in srows}):
                    sub = [r for r in srows if fam in (ALL, r["family"]) and ic in (ALL, r["instrument_class"])]
                    if not sub:
                        continue
                    c = dict(source=source, variant=v, weight=w, set=s, family=fam, instrument_class=ic, **cell(sub, w is not None))
                    cells.append(suppress_cell(c, source))
    return cells


def choose_weight(top1_by_variant, weights):
    """The weight with the most top-1 hits on the selected spectra (ties to the smaller weight), only if it beats
    rules off; else None and the VS21 confirmatory run is skipped."""
    off = top1_by_variant["off"]
    best = None
    for w in sorted(weights):
        k = top1_by_variant[variant_name(w)]
        if best is None or k > best[1]:
            best = (w, k)
    if best is None or best[1] <= off:
        return None
    return best[0]


def load_truth(frozen, source, ids):
    import pyarrow.parquet as pq
    t = pq.read_table(p2.row_path(Path(frozen), "prepared", source), columns=["record_id", "key"],
                      filters=[("record_id", "in", list(ids))]).to_pylist()
    return {r["record_id"]: r["key"] for r in t}


def measure(args, run):
    sfx = suffix(args.source)
    sel = p3m.read_rows(args.out / f"selection{sfx}.tsv")
    recs = {int(r["record_id"]): dict(record_id=int(r["record_id"]), set=r["set"], family=r["family"], rule_id=r["rule_id"],
                                      instrument_class=r["instrument_class"], cluster_id=r["cluster_id"]) for r in sel}
    truth = load_truth(args.frozen_run, args.source, recs)
    res_dir = Path(args.console).parent / "Resources"
    icd_by_block, name_to_id = read_icd(res_dir / RESOURCES[1])
    parents, onto_names = read_obo(args.chemont_obo)
    truth_class = {i: class_level(icd_by_block.get(k, (None, None))[1], parents) for i, k in truth.items()}
    cands, diag = {}, {}
    for w in (None,) + tuple(args.weights):
        v = variant_name(w)
        cands[v] = {}
        for d in sorted((args.out / "predict" / v).iterdir()):
            out = d / "out"
            cands[v].update(parse_structure_result(next(out.glob("Structure result-*.txt"))))
            if w is not None:
                diag.setdefault(v, {}).update(parse_diagnostic_result(next(out.glob("Diagnostic rule result-*.txt"))))
    rows_by_variant, checks = {}, Counter()
    for w in (None,) + tuple(args.weights):
        v = variant_name(w)
        rows = []
        for i, rec in sorted(recs.items()):
            row = spectrum_row(rec, truth[i], truth_class[i], cands[v].get(i, []),
                               None if w is None else cands["off"].get(i, []), parents, name_to_id, onto_names)
            if w is not None:
                outs = [o for o in diag[v].get(i, []) if o["Rule id"] in run["selection"]["reviewed_rules"]
                        and o["Outcome"] in ("fired", "veto_blocked")]
                row["predict_reviewed_outcome"] = bool(outs)
                row["predict_claim_eligible"] = any(o["Claim eligible"] == "true" for o in outs)
                row["predict_instrument_class"] = (diag[v].get(i) or [{}])[0].get("Instrument class", "")
                checks[f"{v}:instrument_class_mismatch"] += row["predict_instrument_class"] not in ("", rec["instrument_class"])
                checks[f"{v}:{rec['set']}:predict_reviewed_outcome"] += row["predict_reviewed_outcome"]
                checks[f"{v}:{rec['set']}:predict_claim_eligible"] += row["predict_claim_eligible"]
                checks[f"{v}:base_rank_mismatch"] += row["base_rank_mismatch"]
                if rec["set"] == "control":
                    checks[f"{v}:control_order_changed"] += row["order_changed"]
            checks[f"{v}:no_candidates"] += row["n_candidates"] == 0
            rows.append(row)
        rows_by_variant[v] = rows
    cells = aggregate(rows_by_variant, args.source, args.weights)
    p3m.write_rows(args.out / f"msfinder_cells{sfx}.csv", cells, CELL_FIELDS)
    spectrum_fields = ["record_id", "set", "family", "rule_id", "instrument_class", "cluster_id", "variant", "rank", "found",
                       "top1", "top5", "n_candidates", "class_defined", "class_correct", "boost_wrong", "rank_changed",
                       "order_changed", "base_rank_mismatch", "predict_reviewed_outcome", "predict_claim_eligible",
                       "predict_instrument_class"]
    p3m.write_rows(args.out / f"per_spectrum{sfx}.tsv",
                   [dict(r, variant=v) for v, rows in rows_by_variant.items() for r in rows], spectrum_fields)
    for f in spectrum_fields + CELL_FIELDS:
        if re.search(p3m.FORBIDDEN_COLUMN, f):
            raise RuntimeError("Forbidden column in an output: " + f)
    top1 = {v: sum(r["top1"] for r in rows if r["set"] == "selected") for v, rows in rows_by_variant.items()}
    result = dict(source=args.source, smoke=args.limit is not None, weights=list(args.weights), top1_selected=top1,
                  n_selected=sum(1 for r in recs.values() if r["set"] == "selected"),
                  checks=dict(sorted(checks.items())),
                  truth_class_defined=sum(c is not None for c in truth_class.values()))
    if args.source == "vs20":
        chosen = choose_weight(top1, args.weights)
        result.update(chosen_weight=chosen, vs21_confirmatory="skipped: no weight beats off" if chosen is None else "run",
                      binding=args.limit is None and tuple(args.weights) == WEIGHTS,
                      rule="most top-1 (first-block) hits on the VS20 selected spectra; ties to the smaller weight; "
                           "must exceed rules off")
    else:
        result["weight"] = args.weights[0]
        if result["n_selected"] < VS21_MIN_CELL:
            result["top1_selected"] = None
    p2.write_json(args.out / f"selection_result{sfx}.json", result)
    run["measure"] = dict(cells=len(cells), checks=result["checks"],
                          chosen_weight=result.get("chosen_weight"), truth_class_defined=result["truth_class_defined"])
    write_report(args.out, run, cells, result, sfx)
    print("measure", args.source, {k: result[k] for k in ("top1_selected", "n_selected")}, flush=True)


def write_report(out, run, cells, result, sfx):
    lines = [f"# WS8 P3 MS-FINDER off/on measurement: {run['run_id']}", "",
             f"Source: {result['source']}. Rule set {EXPECTED_RULESET_VERSION} ({EXPECTED_RULESET_SHA256[:12]}), scope "
             f"{RANKING_SCOPE}. Weights: {', '.join(p3m.num(w) for w in result['weights'])}.", "",
             "## Selection counts", ""]
    lines += [f"- {k}: {v}" for k, v in run["selection"]["counts"].items()]
    lines += ["", "## Selected spectra, all families and instrument classes", "",
              "| variant | n | top-1 | top-5 | not found | class acc. | boost-wrong | cluster top-1 |", "|---|---|---|---|---|---|---|---|"]
    for c in cells:
        if c["set"] == "selected" and c["family"] == ALL and c["instrument_class"] == ALL:
            lines.append(f"| {c['variant']} | {c['n']} | {c['top1_rate']} | {c['top5_rate']} | {c['not_found']} | "
                         f"{c['class_rate']} | {c['boost_wrong_rate']} | {c['c_top1_rate']} |")
    lines += ["", "## Checks", ""] + [f"- {k}: {v}" for k, v in result["checks"].items()]
    if result["source"] == "vs20":
        lines += ["", f"Chosen weight: {result['chosen_weight']} (binding: {result['binding']}); VS21: {result['vs21_confirmatory']}."]
    lines += ["", "Full cells: msfinder_cells" + sfx + ".csv (VS21 values with a denominator below 5 are blanked)."]
    (out / f"report{sfx}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ------------------------------------------------------------------------------------------- guards and main
def check_out(out, smoke, commit, out_root=OUT_ROOT):
    out = Path(out).resolve()
    if out.parent != Path(out_root).resolve():
        raise ValueError("--out must be directly under " + Path(out_root).as_posix())
    m = re.fullmatch(RUN_NAME, out.name)
    if not m:
        raise ValueError("--out must be named p3ms-<yyyymmdd>-<sha> or p3mssmoke-<yyyymmdd>-<sha>")
    if bool(m.group(1)) != smoke:
        raise ValueError("A smoke run (p3mssmoke-) needs --limit and a full run (p3ms-) must not have it")
    if not commit.startswith(m.group(3)):
        raise ValueError("The run id's sha must be a prefix of the code commit")
    return out


def identity_hashes(console, ruleset_dir):
    res_dir = Path(console).parent / "Resources"
    return dict(console_sha256=p2.sha256(console), console_files=p3m.console_files(console),
                method_sha256=hashlib.sha256(method_ini(None).encode("ascii")).hexdigest(),
                resources={n: p2.sha256(res_dir / n) for n in RESOURCES},
                ruleset_files=p3m.ruleset_file_hashes(ruleset_dir))


def vs21_guard(args, out_root=OUT_ROOT):
    """VS21 runs only with --vs21 and --selection-run: a completed, binding VS20 p3ms- run whose recorded chosen weight is
    not None and whose console, resources, method and rule set equal the ones now. Returns the weight."""
    if not getattr(args, "vs21", False):
        raise PermissionError("VS21 needs the explicit --vs21 flag")
    sel = getattr(args, "selection_run", None)
    if not sel:
        raise PermissionError("VS21 needs --selection-run (a completed VS20 p3ms- run)")
    sel = Path(sel).resolve()
    m = re.fullmatch(RUN_NAME, sel.name)
    if sel.parent != Path(out_root).resolve() or not m or m.group(1):
        raise PermissionError("--selection-run must be a full p3ms- run directly under " + Path(out_root).as_posix())
    run_path, res_path = sel / "run.json", sel / "selection_result.json"
    if not run_path.is_file() or not res_path.is_file():
        raise PermissionError("Selection run has no run.json or selection_result.json")
    run = json.loads(run_path.read_text(encoding="utf-8"))
    res = json.loads(res_path.read_text(encoding="utf-8"))
    if run.get("status") != "complete" or run.get("parameters", {}).get("source") != "vs20":
        raise PermissionError("Selection run is not a completed VS20 run")
    if run.get("files", {}).get("selection_result.json") != p2.sha256(res_path):
        raise PermissionError("selection_result.json differs from the one the run recorded")
    if not res.get("binding") or res.get("smoke"):
        raise PermissionError("Selection result is not binding")
    if res.get("chosen_weight") is None:
        raise PermissionError("No weight beat rules off on VS20: the VS21 confirmatory run is skipped")
    now = identity_hashes(args.console, args.ruleset_dir)
    diffs = [k for k in now if run.get("identity", {}).get(k) != now[k]]
    if diffs:
        raise PermissionError("Selection run identity mismatch, VS21 refused: " + ", ".join(diffs))
    return float(res["chosen_weight"])


def detach(argv, out):
    """WMI launch outside the caller's Windows job, as p3_measure.detach, for this script."""
    command = [sys.executable, "-B", str(Path(__file__).resolve())] + [x for x in argv if x != "--detach"] + ["--detached-child"]
    env = dict(os.environ, WS8_P3_COMMAND=subprocess.list2cmdline(command), WS8_P3_CWD=str(ROOT))
    script = ("$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$env:WS8_P3_COMMAND; "
              "CurrentDirectory=$env:WS8_P3_CWD}; if ($r.ReturnValue -ne 0) { throw 'WMI detached launch failed' }; Write-Output $r.ProcessId")
    return subprocess.check_output(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", script],
                                   env=env, text=True, creationflags=p3m.NO_WINDOW).strip()


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("stage", choices=("select", "predict", "measure", "all"))
    parser.add_argument("--out", type=Path, required=True, help="run directory p3ms[smoke]-<yyyymmdd>-<sha> under " + OUT_ROOT.as_posix())
    parser.add_argument("--source", choices=("vs20", "vs21"), default="vs20")
    parser.add_argument("--console", type=Path, required=True, help="MsfinderConsoleApp.exe (0.2.0 build, Resources beside it)")
    parser.add_argument("--ruleset-dir", type=Path, required=True, help="the RuleSet folder of the console's source (0.2.0)")
    parser.add_argument("--chemont-obo", type=Path, required=True, help="ChemOnt_2_1.obo")
    parser.add_argument("--frozen-run", type=Path, default=FROZEN_RUN)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--shards", type=int, default=8, help="shards per variant")
    parser.add_argument("--cores-per-worker", type=int, default=CORES_PER_WORKER)
    parser.add_argument("--limit", type=int, help="smoke: number of selected spectra (control: limit // 4)")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--vs21", action="store_true", help="allow the VS21 source (needs --selection-run)")
    parser.add_argument("--selection-run", type=Path, help="completed VS20 p3ms- run whose chosen weight VS21 uses")
    parser.add_argument("--detach", action="store_true", help="launch outside the caller's Windows job (WMI), log to run.log")
    parser.add_argument("--detached-child", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 8 or args.workers * args.cores_per_worker > 16:
        parser.error("Use 1..8 workers and at most 16 logical processors in all")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.source == "vs21":
        if args.limit is not None:
            parser.error("VS21 has no smoke")
        try:
            args.weights = (vs21_guard(args),)
        except PermissionError as exc:
            parser.error(str(exc))
    else:
        if args.vs21 or args.selection_run:
            parser.error("--vs21 and --selection-run belong to --source vs21")
        args.weights = SMOKE_WEIGHTS if args.limit is not None else WEIGHTS
    return parser, args


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser, args = parse_args(argv)
    commit = p3m.git("rev-parse", "HEAD")
    try:
        args.out = check_out(args.out, args.limit is not None, commit)
    except ValueError as exc:
        parser.error(str(exc))
    if p3m.git("status", "--porcelain", "--untracked-files=no"):
        parser.error("Run from a clean, committed pin")
    if not p3m.git("branch", "-r", "--contains", "HEAD"):
        parser.error("Run from a pushed commit")
    args.out.mkdir(parents=True, exist_ok=True)
    if args.detach:
        print("Detached P3 MS-FINDER process:", detach(argv, args.out), flush=True)
        return
    if args.detached_child:
        sys.stdout = open(args.out / "run.log", "a", encoding="utf-8", buffering=1)
        sys.stderr = sys.stdout
    run_path = args.out / "run.json"
    run = json.loads(run_path.read_text(encoding="utf-8")) if run_path.exists() else {}
    if run.get("status") == "complete":
        parser.error("Completed run is immutable; use a new code commit/run ID")
    params = dict(source=args.source, weights=list(args.weights), limit=args.limit, seed=args.seed, workers=args.workers,
                  shards=args.shards, cores_per_worker=args.cores_per_worker, control_n=CONTROL_N,
                  ranking_scope=RANKING_SCOPE, method_ini=method_ini(None), priority="BelowNormal",
                  selection_run=str(args.selection_run) if args.selection_run else None,
                  frozen_run=Path(args.frozen_run).name, vs21_min_cell=VS21_MIN_CELL)
    if run.get("parameters") not in (None, params):
        parser.error("Stage arguments must match the run manifest")
    p2.priority()
    frozen = Path(args.frozen_run)
    sfx = suffix(args.source)
    inputs = {f"frozen/{n}": p2.sha256(frozen / n) for n in (f"prepared{sfx}.parquet", "fires.local-vs21.parquet")}
    inputs["frozen/run.json"] = p2.sha256(frozen / "run.json")
    inputs["chemont_obo"] = p2.sha256(args.chemont_obo)
    import pyarrow
    run.update(label="ws8-p3-msfinder-off-on", run_id=args.out.name, code_commit=commit, pin_path=str(ROOT),
               interpreter=sys.executable, versions=dict(python=sys.version, pyarrow=pyarrow.__version__),
               code_sha256={p.relative_to(ROOT).as_posix(): p2.sha256(p) for p in CODE_FILES},
               console_path=str(args.console), identity=identity_hashes(args.console, args.ruleset_dir),
               ruleset=dict(version=EXPECTED_RULESET_VERSION, sha256=EXPECTED_RULESET_SHA256), inputs=inputs,
               parameters=params, process=dict(detached_child=args.detached_child, **p2.process_context()), status="running")
    run.setdefault("timings_seconds", {})
    p2.write_json(run_path, run)
    stages = ("select", "predict", "measure") if args.stage == "all" else (args.stage,)
    try:
        for stage in stages:
            start = time.perf_counter()
            globals()[stage](args, run)
            run["timings_seconds"][stage] = round(time.perf_counter() - start, 3)
            p2.write_json(run_path, run)
        if args.stage in ("measure", "all"):
            run["status"] = "complete"
        run["files"] = {p.relative_to(args.out).as_posix(): p2.sha256(p) for p in sorted(args.out.glob("*"))
                        if p.is_file() and p.name not in ("run.json", "run.log")}
        p2.write_json(run_path, run)
        print("completed", args.stage, run["timings_seconds"], flush=True)
    except Exception as exc:
        run.update(status="failed", error_type=type(exc).__name__)
        p2.write_json(run_path, run)
        print("FAILED", type(exc).__name__, flush=True)
        import traceback
        for frame in traceback.extract_tb(exc.__traceback__):
            print(Path(frame.filename).name, frame.lineno, frame.name, flush=True)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
