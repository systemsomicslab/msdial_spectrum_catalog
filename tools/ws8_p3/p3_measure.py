"""WS8 P3 rule-level measurement of the MS-FINDER diagnostic-ion rules (VS20 development pass).

Stages: prepare (record units, truth, clusters, allowlisted MSP), fire (the C# console's diagrules on the
representatives and on the synthetic chimeric set), measure (join, counts, intervals, gate, aggregates), all,
and freeze (the hashes the VS21 pass must match). Optional chemistry (RDKit, pyarrow) stays outside src/.
"""
import argparse
from fractions import Fraction
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from decimal import Decimal
import bisect
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from tools.ws8_p2 import p2_diagnostics as p2  # noqa: E402  (also pins the thread variables to 1)
from tools.ws8_p3 import rule_format as rf  # noqa: E402

# Beside P2's run root (tools/ws8_p3 files carry no drive paths); P2's VS20 run holds the exclusions to match.
OUT_ROOT = p2.OUT_ROOT.parent / "ws8-p3"
P2_RUN = p2.OUT_ROOT / "p2-20261113-c1d5dc3" / "run.json"
RUN_NAME = r"p3cs(smoke)?-(\d{8})-([a-f0-9]{7,40})"
FREEZE_NAME = r"freeze-\d{8}-[a-f0-9]{7,40}\.json"
EXPECTED_RULESET_SHA256 = "af4a7ebf2732b30b9610c7e410ba34f50f176b3fc0155d5a4f9a28fb413ab2eb"
EXPECTED_RULESET_VERSION = "0.1.0"
# VS21 is identified by name and sha256 only (AGENTS.md); its files are opened only behind vs21_guard.
VS21 = [("pos", "MSMS-Private-pos-VS21.msp", "0fa0a96bdf1ebbd86349883cb96bb1ff9e327eb6ba7650ee3434056aacb4804f"),
        ("neg", "MSMS-Private-neg-VS21.msp", "7e2581ab31461b5a9102158a6eaf73d0a1702b8564dd75a9281c9b86dc4d132f")]
PROVENANCE_TABLE = HERE / "provenance-p3.tsv"
PROVENANCE_FIELDS = ("comment",)
PROVENANCE_BUCKETS = ("massbank", "gnps", "mona", "respect", "riken", "other", "unknown")
CODE_FILES = (HERE / "p3_measure.py", HERE / "rule_format.py", ROOT / "tools/ws8_p2/p2_diagnostics.py",
              ROOT / "tools/ws8_p1/mcs_core.py", PROVENANCE_TABLE)

# The console's passes (DiagrulesVariants.AllPasses), in run order.
PASSES = ("base", "tol_product_0.003", "tol_product_0.01", "floor_x0.5", "floor_x2", "k_minus_1",
          "loss_ref_hypothesis", "loss_ref_observed", "decoy_plus_5", "decoy_plus_9", "decoy_minus_7",
          "tol_precursor_0.005", "mz_only_curated")
DECOY_PASSES = ("decoy_plus_5", "decoy_plus_9", "decoy_minus_7")


def decoy_fdr_within(decoy_fires, real_fires, limit):
    """The mean over the shifts of decoy fires / real fires is at most ``limit``, decided exactly.

    A float mean of three ratios can land a hair above a bound it meets exactly (6/120 = 0.05), so the
    comparison is made on rationals: sum(decoy) <= limit * n_shifts * real.
    """
    return Fraction(sum(decoy_fires)) <= Fraction(str(limit)) * len(decoy_fires) * real_fires
PRIMARY_PASS = {"adduct_ambiguous": "base", "mz_only": "mz_only_curated"}
FULL_STRATA_PASSES = ("base", "mz_only_curated")
CHIMERIC_LEVELS = (0.1, 0.3, 0.5)
CHIMERIC_WINDOW_DA = 12.5
SEED = 20261118
BOOTSTRAP = 1000
# Pre-registered criteria (MsdialWorkbench PR #908); they may change before the VS21 pass, never after.
CRITERIA = dict(counting_unit="cluster", n_min=100, wilson_target=0.95, source_n_min=30, source_precision=0.95,
                prior_fraction=0.5, prior_target=0.90, decoy_fdr_max=0.05, decoy_shifts_da=[5, 9, -7],
                occupancy_max=0.05, synthetic_precision=0.95, chimeric_levels=list(CHIMERIC_LEVELS),
                interval="wilson95", lr_interval="log95", vs21_min_cell=5,
                stress_grid_priors=[0.5, 0.2, 0.1, 0.05], stress_grid_prevalence_fractions=[0.5, 0.25])
CONSOLE_INI = "MassToleranceType=Da\r\nMs2Tolerance=0.01\r\nIsTmsMeoxDerivative=False\r\n"
BELOW_NORMAL, NO_WINDOW = 0x4000, 0x08000000
FORBIDDEN_COLUMN = r"(^|_)(key|first.?block|smiles|name)($|_)"
NOT_DONE = [
    "in-source and adduct recast synthetic set (no insource-p3.tsv exists; it is a scientific input for the author)",
    "cross-library VS20/VS21 cluster matching by spectral cosine >= 0.9 and precursor within 0.01 Da",
    "all-spectra sensitivity variant weighted 1/(spectra per compound) (only representatives are fired)",
    "as_applied precision with conflict suppression and nesting (needs the passing rules first)",
    "MS-FINDER-level measurement (option off versus on)",
    "runtime measurement on the WS5 development set",
    "conditional mining (mined_conditional.csv)",
    "precondition decoy (report-only)",
    "figures/",
]


# ------------------------------------------------------------------------------------------- small helpers
def num(x):
    """Fixed-point text of a float's shortest repr (the console's peak grammar has no negative exponent)."""
    return format(Decimal(repr(float(x))), "f")


def write_rows(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n", delimiter="\t" if str(path).endswith(".tsv") else ",")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t" if str(path).endswith(".tsv") else ","))


def suffix(source):
    return ".local-vs21" if source == "vs21" else ""


def any_vs21(sources):
    return "vs21" in sources


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


# ------------------------------------------------------------------------------------- reader and provenance
def load_provenance(path=PROVENANCE_TABLE):
    """Declared provenance regexes: (order, field, pattern, outcome) with outcome a bucket or in_silico."""
    rows = sorted(read_rows(path), key=lambda r: int(r["order"]))
    out = []
    for r in rows:
        if r["field"] not in PROVENANCE_FIELDS:
            raise ValueError("Provenance field is not allowlisted: " + r["field"])
        if r["outcome"] not in PROVENANCE_BUCKETS + ("in_silico",) or r["outcome"] in ("other", "unknown"):
            raise ValueError("Unknown provenance outcome: " + r["outcome"])
        out.append((r["field"], re.compile(r["pattern"]), r["outcome"]))
    return out


def provenance_of(fields, table):
    """Only the bucket and the in-silico flag leave this function; the field text is never kept."""
    present = any(fields.get(f, "").strip() for f in PROVENANCE_FIELDS)
    in_silico = any(o == "in_silico" and p.search(fields.get(f, "")) for f, p, o in table)
    bucket = next((o for f, p, o in table if o != "in_silico" and p.search(fields.get(f, ""))), None)
    return dict(bucket=bucket or ("other" if present else "unknown"), in_silico=bool(in_silico))


def iter_msp_p3(path, table):
    """P2's iter_msp (same allowlist, record splitting and peak parsing; a test checks equality), plus the
    provenance outcome of the declared allowlisted source fields. The field text is reduced to a bucket at once."""
    allowed = {"precursormz", "precursortype", "ionmode", "instrumenttype", "collisionenergy", "smiles"}
    meta, peaks, decimals, extra, in_peaks = {}, [], [], {}, False
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                if meta:
                    yield meta, peaks, decimals, provenance_of(extra, table)
                meta, peaks, decimals, extra, in_peaks = {}, [], [], {}, False
            elif in_peaks:
                try:
                    tokens = line.split()
                    mz, intensity = float(tokens[0]), float(tokens[1])
                    if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
                        peaks.append((mz, intensity))
                        mantissa = tokens[0].lower().split("e")[0]
                        decimals.append(len(mantissa.split(".")[1]) if "." in mantissa else 0)
                except (ValueError, IndexError):
                    pass
            elif ":" in line:
                k, v = line.split(":", 1)
                k = k.lower().strip()
                if k in allowed:
                    meta[k] = v.strip()
                elif k in PROVENANCE_FIELDS:
                    extra[k] = v.strip()
                if k in ("num peaks", "numpeaks"):
                    in_peaks = True
        if meta:
            yield meta, peaks, decimals, provenance_of(extra, table)


# ------------------------------------------------------------------------------------------- rule-set tables
def read_ruleset_tables(ruleset_dir):
    d = Path(ruleset_dir)
    motifs = json.loads((d / "tables/motifs.json").read_text(encoding="utf-8"))["motifs"]
    adducts = json.loads((d / "tables/adducts.json").read_text(encoding="utf-8"))["adducts"]
    return [dict(motif_id=m["motif_id"], smarts=m["smarts"]) for m in motifs], adducts


def alias_map(adducts):
    """The console's DiagrulesAdductAliases: canonical name and aliases, ASCII spaces stripped, first row wins."""
    out = {}
    for row in adducts:
        for spelling in [row["adduct"]] + list(row["aliases"]):
            out.setdefault(spelling.replace(" ", ""), (row["adduct"], int(row["charge"])))
    return out


def adduct_charge(raw, mapped):
    if mapped:
        return abs(mapped[1])
    m = re.search(r"\](\d+)[+-]$", raw)
    return int(m.group(1)) if m else 1


def ruleset_file_hashes(ruleset_dir):
    d = Path(ruleset_dir)
    return {p.relative_to(d).as_posix(): p2.sha256(p) for p in sorted(d.rglob("*")) if p.is_file()}


# ------------------------------------------------------------------------------------------------- prepare
def _chem_init(pairs, set_priority):
    if set_priority:
        p2.priority()
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    p2.QUERIES = [(mid, Chem.MolFromSmarts(s)) for mid, s in pairs]
    if any(q is None for _, q in p2.QUERIES):
        raise ValueError("Invalid motif SMARTS")


def compute_chemistry(smiles, motif_rows, workers):
    """P1 standardisation, first block, Hill formula and the shipped motif SMARTS (P2's chemistry())."""
    pairs = [(r["motif_id"], r["smarts"]) for r in motif_rows]
    if workers <= 1:
        _chem_init(pairs, False)
        return {s: p2.chemistry(s) for s in smiles}
    with mp.get_context("spawn").Pool(workers, initializer=_chem_init, initargs=(pairs, True)) as pool:
        return dict(zip(smiles, pool.imap(p2.chemistry, smiles, chunksize=64)))


RAW_FIELDS = ("record_id", "raw_smiles", "mode", "adduct", "precursor", "instrument", "energy",
              "instrument_type", "collision_energy", "bucket", "in_silico", "mz", "intensity")


def _schemas():
    import pyarrow as pa
    raw = pa.schema([("record_id", pa.int64()), ("raw_smiles", pa.string()), ("mode", pa.string()),
                     ("adduct", pa.string()), ("precursor", pa.float64()), ("instrument", pa.string()),
                     ("energy", pa.string()), ("instrument_type", pa.string()), ("collision_energy", pa.string()),
                     ("bucket", pa.string()), ("in_silico", pa.bool_()),
                     ("mz", pa.list_(pa.float64())), ("intensity", pa.list_(pa.float64()))])
    prepared = pa.schema([("record_id", pa.int64()), ("source", pa.string()), ("key", pa.string()),
                          ("formula", pa.string()), ("motifs", pa.list_(pa.string())), ("cluster_id", pa.string()),
                          ("mode", pa.string()), ("adduct", pa.string()), ("adduct_mapped", pa.bool_()),
                          ("precursor", pa.float64()), ("instrument", pa.string()), ("energy", pa.string()),
                          ("instrument_type", pa.string()), ("collision_energy", pa.string()),
                          ("bucket", pa.string()), ("spectrum_count", pa.int64()), ("msp_file", pa.string()),
                          ("msp_index", pa.int64()), ("mz", pa.list_(pa.float64())),
                          ("intensity", pa.list_(pa.float64()))])
    return raw, prepared


def prepare_source(source, inputs, out, limit, workers, motif_rows, aliases, table, shards, public_keys=None):
    """One library source. inputs: [(mode, path, expected sha256 or None)].

    The P2 counters (same names, same order of checks) go to `counts`; P3's own exclusions go to `extra`.
    Returns (counts, extra, stats, public keys, input hashes)."""
    raw_schema, prepared_schema = _schemas()
    out = Path(out)
    counts, extra, hashes = Counter(), Counter(), []
    raw_path = p2.row_path(out, "p3raw", source)
    sink = p2.ParquetRows(raw_path, raw_schema)
    raw_smiles = set()
    for mode, path, expected in inputs:
        digest = p2.sha256(path)
        if expected and digest != expected:
            raise RuntimeError("Input checksum mismatch")
        hashes.append(dict(name=Path(path).name, sha256=digest))
        for ordinal, (meta, peaks, decimals, prov) in enumerate(iter_msp_p3(path, table)):
            prefix = source + ":" + mode + ":"
            counts[prefix + "input"] += 1
            raw_smiles.add(meta.get("smiles", ""))
            if limit and ordinal % max(1, round(1 / limit)) != 0:
                counts[prefix + "smoke_unsampled"] += 1
                continue
            if not peaks:
                counts[prefix + "no_peaks"] += 1
                continue
            if p2.nominal_mass(peaks, decimals):
                counts[prefix + "nominal_mass"] += 1
                continue
            try:
                precursor = float(meta.get("precursormz", ""))
                if not math.isfinite(precursor) or precursor <= 0:
                    raise ValueError
            except ValueError:
                counts[prefix + "precursor"] += 1
                continue
            if p2.ion_mode(meta.get("ionmode", ""), mode) != mode:
                counts[prefix + "mode_mismatch"] += 1
                continue
            adduct = meta.get("precursortype", "").replace(" ", "")
            if not adduct:
                counts[prefix + "missing_precursor_type"] += 1
                continue
            sink.add(dict(record_id=ordinal * 2 + (mode == "neg"), raw_smiles=meta.get("smiles", ""), mode=mode,
                          adduct=adduct, precursor=precursor, instrument=p2.instrument_bin(meta.get("instrumenttype", "")),
                          energy=p2.energy_bin(meta.get("collisionenergy", "")),
                          instrument_type=meta.get("instrumenttype", ""), collision_energy=meta.get("collisionenergy", ""),
                          bucket=prov["bucket"], in_silico=prov["in_silico"],
                          mz=[p[0] for p in peaks], intensity=[p[1] for p in peaks]))
    sink.close()
    chem = compute_chemistry(sorted(raw_smiles), motif_rows, workers)
    for result in chem.values():
        if "rejected" in result:
            counts[source + ":unique_smiles:" + result["rejected"]] += 1
    if source == "vs20":
        public_keys = {r["key"] for r in chem.values() if "key" in r}
    reps, ns = {}, Counter()

    def unit_of(row, c):
        mapped = aliases.get(row["adduct"])
        return (c["key"], row["mode"], mapped[0] if mapped else row["adduct"]), mapped

    for row in p2.parquet_rows(raw_path):
        c = chem[row["raw_smiles"]]
        if "rejected" in c:
            counts[source + ":structure_rejected_records"] += 1
            continue
        if source == "vs21" and c["key"] in public_keys:
            counts[source + ":vs20_overlap_records"] += 1
            continue
        unit, mapped = unit_of(row, c)
        if adduct_charge(row["adduct"], mapped) >= 2:
            extra[source + ":charge_ge2_records"] += 1
            continue
        if row["in_silico"]:
            extra[source + ":in_silico_records"] += 1
            continue
        if not mapped:
            extra[source + ":adduct_unmapped_records"] += 1
        extra[source + ":bucket:" + row["bucket"]] += 1
        ns[unit] += 1
        candidate = (len(row["mz"]), -row["record_id"])
        if unit not in reps or candidate > reps[unit]:
            reps[unit] = candidate
    chosen = {-v[1]: u for u, v in reps.items()}
    rows = []
    for row in p2.parquet_rows(raw_path):
        unit = chosen.get(row["record_id"])
        if unit is None:
            continue
        c = chem[row["raw_smiles"]]
        mapped = aliases.get(row["adduct"])
        rows.append(dict(record_id=row["record_id"], source=source, key=c["key"], formula=c["formula"],
                         motifs=sorted(c["motifs"]), mode=row["mode"], adduct=unit[2], adduct_mapped=bool(mapped),
                         precursor=row["precursor"], instrument=row["instrument"], energy=row["energy"],
                         instrument_type=row["instrument_type"], collision_energy=row["collision_energy"],
                         bucket=row["bucket"], spectrum_count=ns[unit], mz=row["mz"], intensity=row["intensity"]))
    raw_path.unlink()
    rows.sort(key=lambda r: r["record_id"])
    cluster_keys = sorted({(r["formula"], r["adduct"], tuple(r["motifs"])) for r in rows})
    cluster_ids = {k: f"{source}-c{i:06d}" for i, k in enumerate(cluster_keys)}
    msp_dir = out / "msp" / source
    msp_dir.mkdir(parents=True, exist_ok=True)
    blocks = shard_blocks(len(rows), shards)
    index = []
    for s, (lo, hi) in enumerate(blocks):
        name = f"representatives-{s:02d}{suffix(source)}.msp"
        for i, r in enumerate(rows[lo:hi]):
            r["msp_file"], r["msp_index"] = name, i
            r["cluster_id"] = cluster_ids[(r["formula"], r["adduct"], tuple(r["motifs"]))]
            index.append(dict(file_name=name, record_index=i, record_id=r["record_id"]))
        write_msp(msp_dir / name, rows[lo:hi])
    write_rows(msp_dir / f"representatives-index{suffix(source)}.tsv", index, ["file_name", "record_index", "record_id"])
    sink = p2.ParquetRows(p2.row_path(out, "prepared", source), prepared_schema)
    for r in rows:
        sink.add(r)
    sink.close()
    stats = dict(representatives=len(rows), compounds=len({r["key"] for r in rows}), clusters=len(cluster_ids),
                 spectra=sum(ns.values()), shards=len(blocks))
    return counts, extra, stats, public_keys, hashes


def shard_blocks(n, shards):
    shards = max(1, min(shards, n)) if n else 1
    size = -(-n // shards) if n else 0
    return [(i * size, min(n, (i + 1) * size)) for i in range(shards) if i * size < n] or [(0, 0)]


def write_msp(path, rows):
    """The allowlisted MSP for diagrules: NAME carries only the record id; the console reads the six keys."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for r in rows:
            f.write(f"NAME: rec{r['record_id']}\n")
            f.write(f"PRECURSORMZ: {num(r['precursor'])}\n")
            f.write(f"PRECURSORTYPE: {r['adduct']}\n")
            f.write("IONMODE: " + ("Positive" if r["mode"] == "pos" else "Negative") + "\n")
            if r.get("instrument_type"):
                f.write(f"INSTRUMENTTYPE: {r['instrument_type']}\n")
            if r.get("collision_energy"):
                f.write(f"COLLISIONENERGY: {r['collision_energy']}\n")
            f.write(f"Num Peaks: {len(r['mz'])}\n")
            for mz, intensity in zip(r["mz"], r["intensity"]):
                f.write(f"{num(mz)}\t{num(intensity)}\n")
            f.write("\n")


def p2_parity(counts, p2_exclusions, source):
    """P3's P2-named counters must equal P2's run.json exclusions for the same source; a mismatch is a hard stop."""
    mine = {k: v for k, v in counts.items() if k.startswith(source + ":")}
    theirs = {k: v for k, v in p2_exclusions.items() if k.startswith(source + ":")}
    diffs = {k: [mine.get(k, 0), theirs.get(k, 0)] for k in sorted(set(mine) | set(theirs)) if mine.get(k, 0) != theirs.get(k, 0)}
    if diffs:
        raise RuntimeError("P2 exclusion parity mismatch: " + json.dumps(diffs, sort_keys=True))
    return dict(result="match", keys=len(theirs))


def library_inputs(source):
    if source == "vs20":
        return [(mode, p2.LIB_ROOT / rel, sha) for mode, (rel, sha) in zip(("pos", "neg"), p2.PUBLIC)]
    return [(mode, next(p2.LIB_ROOT.rglob(name)), sha) for mode, name, sha in VS21]


def prepare(args, run, inputs=None, p2_run=P2_RUN):
    motif_rows, adducts = read_ruleset_tables(args.ruleset_dir)
    aliases, table = alias_map(adducts), load_provenance()
    run["inputs"] = [dict(name=p.name, sha256=p2.sha256(p)) for p in (PROVENANCE_TABLE,)]
    run["inputs"] += [dict(name="RuleSet/" + rel, sha256=h) for rel, h in ruleset_file_hashes(args.ruleset_dir).items()
                      if rel in ("tables/motifs.json", "tables/adducts.json")]
    run["prepared"], run["exclusions"], run["p3_exclusions"], run["p2_parity"] = {}, {}, {}, {}
    public_keys = None
    for source in args.sources:
        if source == "vs21":
            vs21_guard(args)
        src_inputs = (inputs or {}).get(source) or library_inputs(source)
        counts, extra, stats, public_keys, hashes = prepare_source(
            source, src_inputs, args.out, args.limit, args.workers, motif_rows, aliases, table, args.workers, public_keys)
        run["inputs"] += hashes
        run["exclusions"].update(dict(sorted(counts.items())))
        run["p3_exclusions"].update(dict(sorted(extra.items())))
        run["prepared"][source] = stats
        if args.limit or p2_run is None:
            run["p2_parity"][source] = dict(result="not_checked", reason="smoke sample" if args.limit else "no P2 run")
        else:
            p2_doc = json.loads(Path(p2_run).read_text(encoding="utf-8"))
            run["inputs"].append(dict(name="P2 run.json", sha256=p2.sha256(p2_run)))
            run["p2_parity"][source] = p2_parity(counts, p2_doc["exclusions"], source)
        buckets = {k.split(":")[-1]: v for k, v in extra.items() if ":bucket:" in k}
        run["prepared"][source]["provenance_unavailable"] = set(buckets) <= {"unknown"}
        print("prepare", source, stats, flush=True)


# -------------------------------------------------------------------------------------------------- fire
def build_chimeric(reps, seed, levels=CHIMERIC_LEVELS, window=CHIMERIC_WINDOW_DA, source="vs20"):
    """Each representative mixed with one seeded same-mode partner of another compound within +-window Da, scaled
    so the partner's highest peak is level x the primary's base peak (precursor +-1 Da excluded). Truth: primary."""
    by_mode = defaultdict(list)
    for r in reps:
        by_mode[r["mode"]].append(r)
    mzs = {}
    for mode in by_mode:
        by_mode[mode].sort(key=lambda r: (r["precursor"], r["record_id"]))
        mzs[mode] = [x["precursor"] for x in by_mode[mode]]
    out, no_partner, no_base = [], 0, 0
    for r in sorted(reps, key=lambda r: r["record_id"]):
        pool, mz = by_mode[r["mode"]], mzs[r["mode"]]
        lo = bisect.bisect_left(mz, r["precursor"] - window)
        hi = bisect.bisect_right(mz, r["precursor"] + window)
        rng = random.Random(f"{seed}|chimeric|{source}|{r['record_id']}")
        partner = None
        for _ in range(32):
            if hi <= lo:
                break
            c = pool[rng.randrange(lo, hi)]
            if c["key"] != r["key"]:
                partner = c
                break
        if partner is None:
            others = [c for c in pool[lo:hi] if c["key"] != r["key"]]
            partner = rng.choice(others) if others else None
        if partner is None:
            no_partner += 1
            continue
        base = max((i for m, i in zip(r["mz"], r["intensity"]) if abs(m - r["precursor"]) > 1.0), default=0.0)
        top = max(partner["intensity"], default=0.0)
        if base <= 0 or top <= 0:
            no_base += 1
            continue
        for level in levels:
            scale = level * base / top
            peaks = sorted(list(zip(r["mz"], r["intensity"])) + [(m, i * scale) for m, i in zip(partner["mz"], partner["intensity"])])
            out.append(dict(record_id=r["record_id"], partner_record_id=partner["record_id"], level=level,
                            mode=r["mode"], precursor=r["precursor"], adduct=r["adduct"],
                            instrument_type=r.get("instrument_type", ""), collision_energy=r.get("collision_energy", ""),
                            mz=[p[0] for p in peaks], intensity=[p[1] for p in peaks]))
    return out, dict(no_partner=no_partner, no_base=no_base, records=len(out))


def write_chimeric(out, source, rows, shards):
    msp_dir = Path(out) / "msp" / source
    msp_dir.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda r: (r["level"], r["record_id"]))
    index = []
    for s, (lo, hi) in enumerate(shard_blocks(len(rows), shards)):
        name = f"chimeric-{s:02d}{suffix(source)}.msp"
        for i, r in enumerate(rows[lo:hi]):
            index.append(dict(file_name=name, record_index=i, level=r["level"], record_id=r["record_id"],
                              partner_record_id=r["partner_record_id"]))
        write_msp(msp_dir / name, rows[lo:hi])
    write_rows(msp_dir / f"chimeric-index{suffix(source)}.tsv", index,
               ["file_name", "record_index", "level", "record_id", "partner_record_id"])


def console_files(console):
    """sha256 of every .exe, .dll and .config in the console's bin folder (recursive, dot folders skipped)."""
    root = Path(console).resolve().parent
    return {q.relative_to(root).as_posix(): p2.sha256(q) for q in sorted(root.rglob("*"))
            if q.is_file() and q.suffix.lower() in (".exe", ".dll", ".config")
            and not any(part.startswith(".") for part in q.relative_to(root).parts)}


def run_console(console, msp, ini, outdir, n_records=None, local_vs21=False):
    """One console run per shard at BelowNormal; a shard folder without a manifest is this run's own partial
    output and is replaced. A complete one is reused only when its manifest validates against the shard."""
    outdir = Path(outdir)
    if outdir.exists() and list(outdir.glob("Diagrules run-*.json")):
        try:
            validate_shard(outdir, msp, n_records, local_vs21)
        except RuntimeError as exc:
            raise RuntimeError(f"Refusing to reuse shard {outdir.name}: {exc}") from None
        return "reused"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    flags = BELOW_NORMAL | NO_WINDOW if sys.platform == "win32" else 0
    with open(outdir / "console.log", "w", encoding="utf-8") as log:
        code = subprocess.Popen([str(console), "diagrules", "-i", str(msp), "-m", str(ini), "-o", str(outdir)],
                                stdout=log, stderr=subprocess.STDOUT, creationflags=flags).wait()
    if code != 0:
        raise RuntimeError(f"diagrules exited {code} on {Path(msp).name}")
    return "ran"


def validate_shard(outdir, msp, n_records, local_vs21):
    outdir = Path(outdir)
    manifests = sorted(outdir.glob("Diagrules run-*.json"))
    if len(manifests) != 1:
        raise RuntimeError("Expected one diagrules manifest in " + outdir.name)
    m = json.loads(manifests[0].read_text(encoding="utf-8"))
    rs = m["rule_set"]
    if rs["sha256"] != EXPECTED_RULESET_SHA256 or not rs["pinned"] or rs["version"] != EXPECTED_RULESET_VERSION:
        raise RuntimeError("Console rule set is not the pinned 0.1.0 set")
    if bool(m["local_vs21"]) != local_vs21:
        raise RuntimeError("Console local_vs21 marking disagrees with the source")
    if len(m["inputs"]) != 1 or m["inputs"][0]["sha256"] != p2.sha256(msp) or             (n_records is not None and m["inputs"][0]["records"] != n_records):
        raise RuntimeError("Console manifest input disagrees with the shard MSP")
    if set(m["passes"]) != set(PASSES):
        raise RuntimeError("Console passes differ from the expected 13")
    for p in PASSES:
        if not (outdir / m["passes"][p]["evidence"]).is_file():
            raise RuntimeError("Missing evidence file for pass " + p)
    return dict(manifest=manifests[0].name, sha256=p2.sha256(manifests[0]), engine_semantics=m["engine_semantics"],
                build_identity=m["build_identity"],
                refusals={p: m["passes"][p]["refusals"] for p in PASSES if m["passes"][p]["refusals"]},
                counts={p: m["passes"][p]["counts"] for p in PASSES})


def check_shards(results, frozen_identity=None):
    """Every shard ran the same console build and engine, without record errors; for VS21 the frozen one."""
    ids = {(r["build_identity"], r["engine_semantics"]) for r in results.values()}
    if len(ids) > 1:
        raise RuntimeError("Shards disagree on the console build identity or engine semantics")
    if frozen_identity is not None and ids and ids != {(frozen_identity["build_identity"], frozen_identity["engine_semantics"])}:
        raise RuntimeError("Console build identity or engine semantics differ from the freeze")
    errors = sorted(f"{k}:{p}" for k, r in results.items() for p, c in r["counts"].items() if c["errors"] > 0)
    if errors:
        raise RuntimeError("Console reported record errors in " + ", ".join(errors[:5]))
    return dict(build_identity=next(iter(ids))[0] if ids else None, engine_semantics=next(iter(ids))[1] if ids else None)


def console_identity(console):
    """build_identity and engine_semantics of the console, from one synthetic record in a temporary folder."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        msp, ini, out = Path(d) / "probe.msp", Path(d) / "probe.ini", Path(d) / "out"
        write_msp(msp, [dict(record_id=0, precursor=200.1, adduct="[M+H]+", mode="pos", mz=[100.1], intensity=[10.0])])
        ini.write_text(CONSOLE_INI, encoding="ascii", newline="")
        run_console(console, msp, ini, out)
        r = validate_shard(out, msp, 1, False)
        return dict(build_identity=r["build_identity"], engine_semantics=r["engine_semantics"])


def shard_jobs(out, source, kind):
    msp_dir = Path(out) / "msp" / source
    index = read_rows(msp_dir / f"{kind}-index{suffix(source)}.tsv")
    per_file = Counter(r["file_name"] for r in index)
    return [(msp_dir / name, Path(out) / "fire" / source / kind / Path(name).stem, n) for name, n in sorted(per_file.items())]


def fire(args, run):
    console = Path(args.console)
    ini = args.out / "diagrules.ini"
    if not ini.exists():
        ini.write_text(CONSOLE_INI, encoding="ascii", newline="")
    run["console"] = dict(path=str(console), sha256=p2.sha256(console), ini_sha256=p2.sha256(ini),
                          files=console_files(console))
    run["chimeric"] = {}
    jobs = []
    for source in args.sources:
        if source == "vs21":
            vs21_guard(args)
        reps = list(p2.parquet_rows(p2.row_path(args.out, "prepared", source)))
        rows, stats = build_chimeric(reps, args.seed, source=source)
        write_chimeric(args.out, source, rows, args.workers)
        run["chimeric"][source] = stats
        for kind in ("representatives", "chimeric"):
            jobs += [(source, kind, *j) for j in shard_jobs(args.out, source, kind)]
    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(run_console, console, msp, ini, outdir, n, source == "vs21"): (source, kind, msp, outdir, n)
                   for source, kind, msp, outdir, n in jobs}
        for fut, (source, kind, msp, outdir, n) in futures.items():
            state = fut.result()
            results[f"{source}/{kind}/{outdir.name}"] = dict(state=state, **validate_shard(outdir, msp, n, source == "vs21"))
    run["fire"] = results
    frozen = None
    if any_vs21(args.sources):
        frozen = json.loads(Path(args.freeze).read_text(encoding="utf-8"))["console_identity"]
    run["console"]["identity"] = check_shards(results, frozen)
    print("fire", len(results), "shards", flush=True)


# ----------------------------------------------------------------------------------------------- evidence
def read_evidence(shard_dir, msp_name, index_ids, passes=PASSES):
    """{pass: [(uid, status, instrument_class, [(rule_id, version, fired, n_hyp, veto, adduct)])]} in index order.
    index_ids: uid per record_index. Every index row must have exactly one line per pass, in order."""
    shard_dir = Path(shard_dir)
    out = {}
    for p in passes:
        files = sorted(shard_dir.glob(f"Diagrules evidence-{p}-*.jsonl"))
        files = [f for f in files if re.fullmatch(rf"Diagrules evidence-{re.escape(p)}-\d+(\.local-vs21)?\.jsonl", f.name)]
        if len(files) != 1:
            raise RuntimeError(f"Expected one evidence file for pass {p} in {shard_dir.name}")
        recs = []
        with open(files[0], encoding="utf-8") as f:
            for i, line in enumerate(f):
                e = json.loads(line)
                if e["file_name"] != msp_name or e["record_index"] != i or i >= len(index_ids):
                    raise RuntimeError(f"Evidence does not join the index ({p}, line {i})")
                res = e["result"] or {}
                outcomes = [(o["rule_id"], o["version"], bool(o["fired"]), o["n_hypotheses"], bool(o["veto_blocked"]),
                             (o["chosen"] or {}).get("adduct")) for o in res.get("outcomes", []) if o["n_hypotheses"] > 0]
                recs.append((index_ids[i], e["status"], res.get("instrument_class"), outcomes))
        if len(recs) != len(index_ids):
            raise RuntimeError(f"Evidence has {len(recs)} lines for {len(index_ids)} index rows ({p})")
        out[p] = recs
    return out


def collect(evidence_by_shard, versions=None):
    """by_rule[pass][rule_id] = [(uid, fired, n_hyp, veto, adduct)]; notev[pass] = uids not evaluated."""
    by_rule = defaultdict(lambda: defaultdict(list))
    notev = defaultdict(set)
    inst = {}
    for ev in evidence_by_shard:
        for p, recs in ev.items():
            for uid, status, iclass, outcomes in recs:
                if status != "evaluated":
                    notev[p].add(uid)
                if p == "base":
                    inst[uid] = iclass
                for rid, ver, fired, nh, veto, adduct in outcomes:
                    if versions is not None and versions.get(rid) not in (None, ver):
                        raise RuntimeError("Evidence rule version differs from the rule file: " + rid)
                    by_rule[p][rid].append((uid, fired, nh, veto, adduct))
    return by_rule, notev, inst


# ------------------------------------------------------------------------------------------------ statistics
def stress(lr_lo, prior):
    if lr_lo is None or prior is None:
        return None
    return prior * lr_lo / (1 - prior + prior * lr_lo)


def percentile_pair(values, lo=0.025, hi=0.975):
    if not values:
        return None, None
    v = sorted(values)
    n = len(v)
    return v[min(n - 1, int(math.floor(lo * n)))], v[max(0, int(math.ceil(hi * n)) - 1)]


def cluster_bootstrap(fired_types, n_clusters, replicates, rng):
    """Resample n_clusters clusters with replacement (multinomial by conditional binomials over cluster types).
    fired_types: Counter{(truth, fired_units): count} of fired clusters; the rest of n_clusters never fired.
    Returns record-level and cluster-level precision intervals and the count of undefined replicates."""
    types = sorted(fired_types.items())
    if not types or n_clusters <= 0:
        return dict(boot_lo=None, boot_hi=None, c_boot_lo=None, c_boot_hi=None, boot_undefined=replicates if replicates else 0)
    rec, clu, undefined = [], [], 0
    for _ in range(replicates):
        remaining_n, remaining_c = n_clusters, n_clusters
        tp_u = fired_u = tp_c = fired_c = 0
        for (truth, units), c in types:
            x = rng.binomialvariate(remaining_n, c / remaining_c) if remaining_n and remaining_c else 0
            remaining_n -= x
            remaining_c -= c
            fired_c += x
            fired_u += x * units
            if truth:
                tp_c += x
                tp_u += x * units
        if fired_c == 0:
            undefined += 1
            continue
        rec.append(tp_u / fired_u)
        clu.append(tp_c / fired_c)
    (a, b), (c, d) = percentile_pair(rec), percentile_pair(clu)
    return dict(boot_lo=a, boot_hi=b, c_boot_lo=c, c_boot_hi=d, boot_undefined=undefined)


def axis_occupancy(precursors_pda, tol_pda):
    """Fraction of the rule's m/z axis within tolerance of any hypothesis precursor. The axis is the nominal-mass
    span of the series, [floor(min), floor(max) + 1) Da, so a single hypothesis spans one dalton."""
    if not precursors_pda:
        return dict(occupancy=None, n_hypotheses=0, lo_da=None, hi_da=None)
    pda = rf.PDA_PER_DA
    lo, hi = (min(precursors_pda) // pda) * pda, (max(precursors_pda) // pda + 1) * pda
    covered, end = 0, lo
    for a, b in sorted((max(lo, p - tol_pda), min(hi, p + tol_pda)) for p in precursors_pda):
        if b > end:
            covered += b - max(a, end)
            end = b
    return dict(occupancy=covered / (hi - lo), n_hypotheses=len(precursors_pda), lo_da=lo // pda, hi_da=hi // pda)


# --------------------------------------------------------------------------------------------- measurement
def cells_of(rec, multi):
    sources = ("pooled", rec["source"]) if multi else ("pooled",)
    return [(s, i, e) for s in sources for i in ("all", rec["instrument"]) for e in ("all", rec["energy"])]


def universe(truth_by_mode, mode, motif, multi, sources):
    U = defaultdict(lambda: dict(n=0, pos=0, cl=set(), pcl=set()))
    for s in ("pooled",) + (tuple(sources) if multi else ()):
        U[(s, "all", "all")]
    for uid, rec in truth_by_mode.get(mode, ()):
        t = motif in rec["motifs"]
        for cell in cells_of(rec, multi):
            u = U[cell]
            u["n"] += 1
            u["cl"].add(rec["cluster"])
            if t:
                u["pos"] += 1
                u["pcl"].add(rec["cluster"])
    return U


def conditional(entries, truth, motif, multi, cell_fn=None):
    C = defaultdict(lambda: dict(cond=0, cond_pos=0, fires=0, tp=0, ccl=set(), cpcl=set(), fired=Counter(), uids=set()))
    for uid, fired, *_ in entries:
        rec = truth[uid[:2]]
        t = motif in rec["motifs"]
        for cell in (cell_fn(rec) if cell_fn else cells_of(rec, multi)):
            c = C[cell]
            c["cond"] += 1
            c["ccl"].add(rec["cluster"])
            c["uids"].add(uid)
            if t:
                c["cond_pos"] += 1
                c["cpcl"].add(rec["cluster"])
            if fired:
                c["fires"] += 1
                c["tp"] += t
                c["fired"][rec["cluster"]] += 1
    return C


def cell_row(u, c, notev, replicates, rng):
    """Record- and cluster-level counts and statistics of one cell; u universe, c conditional accumulator."""
    n, pos = u["n"], u["pos"]
    cond, cond_pos, fires, tp = c["cond"], c["cond_pos"], c["fires"], c["tp"]
    fp = fires - tp
    lo, hi = p2.wilson(tp, fires)
    rec = p2.statistics(tp, fp, pos - tp, n - pos - fp)
    c_univ, c_pos, c_cond, c_cond_pos = len(u["cl"]), len(u["pcl"]), len(c["ccl"]), len(c["cpcl"])
    c_fires = len(c["fired"])
    c_tp = sum(1 for k in c["fired"] if k in u["pcl"])
    c_fp = c_fires - c_tp
    clo, chi = p2.wilson(c_tp, c_fires)
    cst = p2.statistics(c_tp, c_fp, c_pos - c_tp, c_univ - c_pos - c_fp)
    ccond = p2.statistics(c_tp, c_fp, c_cond_pos - c_tp, c_cond - c_cond_pos - c_fp)
    row = dict(n_universe=n, n_not_evaluated=notev, n_positive=pos, n_cond=cond, n_cond_positive=cond_pos,
               fires=fires, tp=tp, fp=fp, precision=tp / fires if fires else None, wilson_lo=lo if fires else None,
               wilson_hi=hi if fires else None, recall_cond=tp / cond_pos if cond_pos else None,
               recall_universe=tp / pos if pos else None, coverage=cond_pos / pos if pos else None,
               lr_plus=rec["lr_plus"], lr_lo=rec["lr_lo"], prior_stress=rec["prior_stress_precision"],
               c_universe=c_univ, c_positive=c_pos, c_cond=c_cond, c_cond_positive=c_cond_pos, c_fires=c_fires,
               c_tp=c_tp, c_fp=c_fp, c_precision=c_tp / c_fires if c_fires else None,
               c_wilson_lo=clo if c_fires else None, c_wilson_hi=chi if c_fires else None,
               c_recall_cond=c_tp / c_cond_pos if c_cond_pos else None, c_coverage=c_cond_pos / c_pos if c_pos else None,
               c_lr_plus=cst["lr_plus"], c_lr_lo=cst["lr_lo"], c_prevalence=cst["prevalence"],
               c_prior_stress=cst["prior_stress_precision"], c_cond_prevalence=ccond["prevalence"],
               c_cond_lr_lo=ccond["lr_lo"])
    for pi in CRITERIA["stress_grid_priors"]:
        row[f"stress_pi_{pi}"] = stress(ccond["lr_lo"], pi)
    for frac in CRITERIA["stress_grid_prevalence_fractions"]:
        row[f"stress_cond_prev_x{frac}"] = stress(ccond["lr_lo"], frac * ccond["prevalence"] if ccond["prevalence"] is not None else None)
    if rng is not None:
        types = Counter((k in u["pcl"], units) for k, units in c["fired"].items())
        row.update(cluster_bootstrap(types, c_cond, replicates, rng))
    else:
        row.update(boot_lo=None, boot_hi=None, c_boot_lo=None, c_boot_hi=None, boot_undefined=None)
    return row


EMPTY_C = dict(cond=0, cond_pos=0, fires=0, tp=0, ccl=set(), cpcl=set(), fired=Counter(), uids=set())


def mode_of_pass(p):
    return "mz_only" if p == "mz_only_curated" else "adduct_ambiguous"


def measure_core(rules, truth, by_rule, notev, sources, seed=SEED, replicates=BOOTSTRAP, chimeric=None, refused=()):
    """Pure measurement. rules: dicts (rule_id, version, ion_mode, motif_id, competitors, measure, source_checked,
    primary_mode, origin, axis{precondition mode: axis_occupancy()}). truth: {(source, record_id): dict(mode,
    motifs, cluster, instrument, energy, bucket, source)}. by_rule[pass][rule_id] = [(uid, fired, ...)] for the
    representatives; chimeric = dict(by_rule={level: by_rule}, universe={level: [uids]}, stats)."""
    multi = len(sources) > 1
    truth_by_mode = defaultdict(list)
    for uid in sorted(truth):
        truth_by_mode[truth[uid]["mode"]].append((uid, truth[uid]))
    cache = {}
    precision_rows, competitor_rows, synthetic_rows, occupancy_rows, gates = [], [], [], [], {}
    for rule in sorted(rules, key=lambda r: r["rule_id"]):
        rid, mode, motif = rule["rule_id"], rule["ion_mode"], rule["motif_id"]
        ids = dict(rule_id=rid, version=rule["version"], ion_mode=mode, motif_id=motif)
        eligible = rule["measure"] and rule["source_checked"] and rid not in refused
        if not eligible:
            reason = "refused" if rid in refused else "measure_false" if not rule["measure"] else "source_unchecked"
            gates[rid] = dict(ids, status="not_measured", reason=reason, origin=rule["origin"])
            continue
        if (mode, motif) not in cache:
            cache[(mode, motif)] = universe(truth_by_mode, mode, motif, multi, sources)
        U = cache[(mode, motif)]
        cells = {}
        for p in PASSES:
            C = conditional(by_rule.get(p, {}).get(rid, []), truth, motif, multi)
            nv = Counter(cell for uid in notev.get(p, ()) if truth[uid[:2]]["mode"] == mode
                         for cell in cells_of(truth[uid[:2]], multi))
            for cell in sorted(U):
                if p not in FULL_STRATA_PASSES and cell[1:] != ("all", "all"):
                    continue
                rng = random.Random(f"{seed}|boot|{rid}|{p}|{'|'.join(cell)}")
                row = dict(ids, pass_id=p, precondition_mode=mode_of_pass(p), source=cell[0], instrument=cell[1],
                           energy=cell[2], **cell_row(U[cell], C.get(cell, EMPTY_C), nv.get(cell, 0), replicates, rng))
                if p == "base":
                    # Decoy fires inside the base conditional population of the cell, in records and in clusters;
                    # the gate counts clusters (PR 908), the record-level FDR is report-only.
                    base_uids = C.get(cell, EMPTY_C)["uids"]
                    fdrs, cfdrs = [], []
                    for d in DECOY_PASSES:
                        hits = [e[0] for e in by_rule.get(d, {}).get(rid, []) if e[1] and e[0] in base_uids]
                        k, ck = len(hits), len({truth[u[:2]]["cluster"] for u in hits})
                        row["decoy_fires_" + d[6:]], row["c_decoy_fires_" + d[6:]] = k, ck
                        fdrs.append(k / row["fires"] if row["fires"] else None)
                        cfdrs.append(ck / row["c_fires"] if row["c_fires"] else None)
                    row["decoy_fdr"] = sum(fdrs) / len(fdrs) if row["fires"] else None
                    row["c_decoy_fdr"] = sum(cfdrs) / len(cfdrs) if row["c_fires"] else None
                cells[(p,) + cell] = row
                precision_rows.append(row)
        primary = PRIMARY_PASS[rule["primary_mode"]]
        gate_row = cells[(primary, "pooled", "all", "all")]
        decoy_row = cells[("base", "pooled", "all", "all")]
        # Occupancy, both definitions, for both precondition modes.
        occ = {}
        for pm, p in PRIMARY_PASS.items():
            r = cells[(p, "pooled", "all", "all")]
            neg, neg_cond = r["n_universe"] - r["n_positive"], r["n_cond"] - r["n_cond_positive"]
            cneg, cneg_cond = r["c_universe"] - r["c_positive"], r["c_cond"] - r["c_cond_positive"]
            ax = rule["axis"].get(pm) or {}
            orow = dict(ids, precondition_mode=pm, axis_occupancy=ax.get("occupancy"), axis_lo_da=ax.get("lo_da"),
                        axis_hi_da=ax.get("hi_da"), n_hypotheses=ax.get("n_hypotheses"), n_positive=r["n_positive"],
                        c_positive=r["c_positive"], n_negative=neg,
                        n_negative_cond=neg_cond, record_occupancy=neg_cond / neg if neg else None,
                        c_negative=cneg, c_negative_cond=cneg_cond, c_record_occupancy=cneg_cond / cneg if cneg else None)
            orow["non_conditional"] = any(v is not None and v > CRITERIA["occupancy_max"]
                                          for v in (orow["axis_occupancy"], orow["record_occupancy"]))
            occ[pm] = orow
            occupancy_rows.append(orow)
        # Competitors, primary pass, pooled.
        entries = by_rule.get(primary, {}).get(rid, [])
        for cm in rule["competitors"]:
            sel = [e for e in entries if cm in truth[e[0][:2]]["motifs"] and motif not in truth[e[0][:2]]["motifs"]]
            n, k = len(sel), sum(1 for e in sel if e[1])
            cn = {truth[e[0][:2]]["cluster"] for e in sel}
            ck = {truth[e[0][:2]]["cluster"] for e in sel if e[1]}
            lo, hi = p2.wilson(k, n)
            competitor_rows.append(dict(ids, pass_id=primary, competitor_motif_id=cm,
                                        n_universe_competitor=sum(1 for _, r in truth_by_mode.get(mode, ())
                                                                  if cm in r["motifs"] and motif not in r["motifs"]),
                                        n_cond=n, fires=k, fpr=k / n if n else None, fpr_wilson_hi=hi if n else None,
                                        c_cond=len(cn), c_fires=len(ck), c_fpr=len(ck) / len(cn) if cn else None))
        # Synthetic chimeric set, primary pass, by the primary record's truth.
        synth = {}
        for level in CHIMERIC_LEVELS:
            if chimeric is None:
                synth[level] = None
                continue
            ents = chimeric["by_rule"].get(level, {}).get(primary, {}).get(rid, [])
            uni = [u for u in chimeric["universe"].get(level, ()) if truth[u[:2]]["mode"] == mode]
            fires = [e for e in ents if e[1]]
            tp = sum(1 for e in fires if motif in truth[e[0][:2]]["motifs"])
            cf = {truth[e[0][:2]]["cluster"]: motif in truth[e[0][:2]]["motifs"] for e in fires}
            lo, hi = p2.wilson(tp, len(fires))
            clo, chi = p2.wilson(sum(cf.values()), len(cf))
            nv = sum(1 for u in chimeric.get("notev", {}).get(level, {}).get(primary, ()) if truth[u[:2]]["mode"] == mode)
            # The gate reads the cluster-level precision (PR 908); the record-level values are report-only.
            srow = dict(ids, pass_id=primary, level=level, n_universe=len(uni), n_not_evaluated=nv, n_cond=len(ents),
                        fires=len(fires), tp=tp, fp=len(fires) - tp, precision=tp / len(fires) if fires else None,
                        wilson_lo=lo if fires else None, wilson_hi=hi if fires else None, c_fires=len(cf),
                        c_tp=sum(cf.values()), c_fp=len(cf) - sum(cf.values()),
                        c_precision=sum(cf.values()) / len(cf) if cf else None,
                        c_wilson_lo=clo if cf else None, c_wilson_hi=chi if cf else None)
            synth[level] = srow
            synthetic_rows.append(srow)
        by_source = {s: cells[(primary, s, "all", "all")] for s in sources} if multi else {sources[0]: gate_row}
        by_prov = {}
        prov_c = conditional(entries, truth, motif, multi, cell_fn=lambda rec: [rec["bucket"]])
        for b in sorted(prov_c):
            c = prov_c[b]
            c_tp = sum(1 for k in c["fired"] if truth_cluster_positive(k, U))
            by_prov[b] = dict(fires=c["fires"], tp=c["tp"], precision=c["tp"] / c["fires"] if c["fires"] else None,
                              c_fires=len(c["fired"]), c_tp=c_tp,
                              c_precision=c_tp / len(c["fired"]) if c["fired"] else None)
        gates[rid] = dict(ids, origin=rule["origin"], primary_precondition_mode=rule["primary_mode"], gate_pass=primary,
                          **decide(gate_row, decoy_row, by_source, occ[rule["primary_mode"]], synth, rule["origin"]),
                          by_provenance=by_prov,
                          provenance_flag=any(v["c_fires"] >= CRITERIA["source_n_min"] and v["c_precision"] < CRITERIA["source_precision"]
                                              for v in by_prov.values()))
    return dict(precision=precision_rows, competitors=competitor_rows, synthetic=synthetic_rows,
                occupancy=occupancy_rows, gates=gates)


def truth_cluster_positive(cluster, U):
    return cluster in U[("pooled", "all", "all")]["pcl"]


def decide(gate_row, decoy_row, by_source, occ, synth, origin):
    """The pre-registered gate on the primary-mode pooled cell, counting clusters."""
    K = CRITERIA
    n = gate_row["c_fires"]

    def at_least(v, t):  # None: not evaluable (no fires), never a pass
        return None if v is None else v >= t

    chim_ok = None if any(s is None or s["n_universe"] == 0 for s in synth.values()) else \
        all(s["c_precision"] is None or s["c_precision"] >= K["synthetic_precision"] for s in synth.values())
    crit = {
        "n_clusters": dict(value=n, threshold=K["n_min"], ok=n >= K["n_min"]),
        "cluster_wilson_lo": dict(value=gate_row["c_wilson_lo"], threshold=K["wilson_target"],
                                  ok=at_least(gate_row["c_wilson_lo"], K["wilson_target"])),
        "per_source_precision": dict(value={s: dict(c_fires=r["c_fires"], c_precision=r["c_precision"]) for s, r in by_source.items()},
                                     threshold=K["source_precision"], min_n=K["source_n_min"],
                                     ok=all(r["c_fires"] < K["source_n_min"] or r["c_precision"] >= K["source_precision"]
                                            for r in by_source.values())),
        "prior_stress": dict(value=gate_row["c_prior_stress"], threshold=K["prior_target"],
                             ok=at_least(gate_row["c_prior_stress"], K["prior_target"])),
        "decoy_fdr": dict(value=decoy_row["c_decoy_fdr"], unit="cluster", record_level=decoy_row["decoy_fdr"],
                          threshold=K["decoy_fdr_max"],
                          ok=None if decoy_row["c_decoy_fdr"] is None else decoy_fdr_within(
                              [decoy_row["c_decoy_fires_" + d[6:]] for d in DECOY_PASSES], decoy_row["c_fires"],
                              K["decoy_fdr_max"])),
        "occupancy": dict(value=dict(axis=occ["axis_occupancy"], records=occ["record_occupancy"]),
                          threshold=K["occupancy_max"], ok=not occ["non_conditional"]),
        "synthetic_chimeric": dict(value={str(lv): (s["c_precision"] if s else None) for lv, s in synth.items()},
                                   unit="cluster", record_level={str(lv): (s["precision"] if s else None) for lv, s in synth.items()},
                                   threshold=K["synthetic_precision"], ok=chim_ok),
        "synthetic_insource": dict(value=None, ok=None, pending="in-source/adduct recast set not built in this phase"),
    }
    if origin == "tuned_vs20":
        crit["heldout"] = dict(ok=None, pending="tuned_vs20: must also pass on VS21 alone")
    elif origin == "tuned_pooled":
        crit["heldout"] = dict(ok=None, pending="tuned_pooled: can pass only on an external set")
    else:
        crit["heldout"] = dict(ok=None, pending="as_published: the gate is on VS20 and VS21 pooled")
    if crit["synthetic_chimeric"]["ok"] is None:
        crit["synthetic_chimeric"]["pending"] = "chimeric set not measured"
    if not crit["n_clusters"]["ok"]:
        status = "not_estimable"
    elif any(c["ok"] is False for c in crit.values()):
        status = "fail"
    elif any(c["ok"] is None for c in crit.values()):
        status = "incomplete"
    else:
        status = "pass"
    return dict(status=status, criteria=crit, non_conditional=occ["non_conditional"],
                failed=[k for k, c in crit.items() if c["ok"] is False],
                pending=[k for k, c in crit.items() if c["ok"] is None])


# ------------------------------------------------------------------------------------------- aggregates out
PRECISION_FIELDS = ["rule_id", "version", "ion_mode", "motif_id", "pass_id", "precondition_mode", "source",
                    "instrument", "energy", "n_universe", "n_not_evaluated", "n_positive", "n_cond", "n_cond_positive",
                    "fires", "tp", "fp", "precision", "wilson_lo", "wilson_hi", "recall_cond", "recall_universe",
                    "coverage", "lr_plus", "lr_lo", "prior_stress", "c_universe", "c_positive", "c_cond",
                    "c_cond_positive", "c_fires", "c_tp", "c_fp", "c_precision", "c_wilson_lo", "c_wilson_hi",
                    "c_recall_cond", "c_coverage", "c_lr_plus", "c_lr_lo", "c_prevalence", "c_prior_stress",
                    "c_cond_prevalence", "c_cond_lr_lo"] + \
                   [f"stress_pi_{pi}" for pi in CRITERIA["stress_grid_priors"]] + \
                   [f"stress_cond_prev_x{f}" for f in CRITERIA["stress_grid_prevalence_fractions"]] + \
                   ["boot_lo", "boot_hi", "c_boot_lo", "c_boot_hi", "boot_undefined",
                    "decoy_fires_plus_5", "decoy_fires_plus_9", "decoy_fires_minus_7", "decoy_fdr",
                    "c_decoy_fires_plus_5", "c_decoy_fires_plus_9", "c_decoy_fires_minus_7", "c_decoy_fdr"]
COMPETITOR_FIELDS = ["rule_id", "version", "ion_mode", "motif_id", "pass_id", "competitor_motif_id",
                     "n_universe_competitor", "n_cond", "fires", "fpr", "fpr_wilson_hi", "c_cond", "c_fires", "c_fpr"]
SYNTHETIC_FIELDS = ["rule_id", "version", "ion_mode", "motif_id", "pass_id", "level", "n_universe", "n_not_evaluated",
                    "n_cond", "fires", "tp", "fp", "precision", "wilson_lo", "wilson_hi", "c_fires", "c_tp", "c_fp",
                    "c_precision", "c_wilson_lo", "c_wilson_hi"]
OCCUPANCY_FIELDS = ["rule_id", "version", "ion_mode", "motif_id", "precondition_mode", "axis_occupancy", "axis_lo_da",
                    "axis_hi_da", "n_hypotheses", "n_positive", "n_negative", "n_negative_cond", "record_occupancy",
                    "c_positive", "c_negative", "c_negative_cond", "c_record_occupancy", "non_conditional"]
# Rates and the count that is their denominator, for VS21 small-cell suppression.
RATE_DENOMINATOR = {"precision": "fires", "wilson_lo": "fires", "wilson_hi": "fires", "boot_lo": "fires",
                    "boot_hi": "fires", "decoy_fdr": "fires", "recall_cond": "n_cond_positive",
                    "recall_universe": "n_positive", "coverage": "n_positive", "c_precision": "c_fires",
                    "c_wilson_lo": "c_fires", "c_wilson_hi": "c_fires", "c_boot_lo": "c_fires", "c_boot_hi": "c_fires",
                    "c_recall_cond": "c_cond_positive", "c_coverage": "c_positive", "fpr": "n_cond",
                    "fpr_wilson_hi": "n_cond", "c_fpr": "c_cond", "record_occupancy": "n_negative",
                    "c_record_occupancy": "c_negative"}
ID_FIELDS = {"rule_id", "version", "ion_mode", "motif_id", "pass_id", "precondition_mode", "source", "instrument",
             "energy", "competitor_motif_id", "level", "axis_lo_da", "axis_hi_da", "n_hypotheses", "axis_occupancy",
             "non_conditional"}


C_RATES = "every c_* value that is not a count"
DECOYS = ("plus_5", "plus_9", "minus_7")
# Complementary sets (counts, rates): one small count blanks the whole set, so no member can be recovered by
# subtraction from the others (for example c_fp = c_fires - c_tp). A count shared by two sets links them.
SUPPRESSION_GROUPS = [
    (("fires", "tp", "fp"), ("precision", "wilson_lo", "wilson_hi", "boot_lo", "boot_hi", "recall_cond",
                             "recall_universe", "fpr", "fpr_wilson_hi")),
    (("c_fires", "c_tp", "c_fp"), C_RATES),
    (("n_universe", "n_positive"), ("coverage", "recall_universe")),
    (("n_cond", "n_cond_positive"), ("recall_cond", "fpr", "fpr_wilson_hi")),
    (("c_universe", "c_positive"), C_RATES),
    (("c_cond", "c_cond_positive"), C_RATES),
    (tuple("decoy_fires_" + d for d in DECOYS), ("decoy_fdr",)),
    (tuple("c_decoy_fires_" + d for d in DECOYS), ("c_decoy_fdr",)),
    (("n_positive", "n_negative", "n_negative_cond"), ("record_occupancy",)),
    (("c_positive", "c_negative", "c_negative_cond"), C_RATES),
]
WHOLE_TABLE = ("lr_", "c_lr_", "prior_stress", "c_prior_stress", "stress_", "c_cond_lr", "c_prevalence", "c_cond_prevalence")


def is_count(k, v):
    return k not in ID_FIELDS and isinstance(v, int) and not isinstance(v, bool)


def suppress(row, n_min=CRITERIA["vs21_min_cell"]):
    """VS21 small cells: a count in 1..n_min-1 is blanked with every count and rate of its complementary sets
    (closed over shared counts), every rate whose denominator is below n_min is blanked, and likelihood ratios and
    stress values (whole-table functions) are blanked when any count of the row is."""
    out = dict(row)
    blank = {k for k, v in row.items() if is_count(k, v) and 0 < v < n_min}
    changed = True
    while changed:
        changed = False
        for counts, rates in SUPPRESSION_GROUPS:
            present = [k for k in counts if k in row]
            if present and blank & set(present):
                members = set(present) | ({k for k in row if k.startswith("c_") and not is_count(k, row[k])}
                                          if rates is C_RATES else {k for k in rates if k in row})
                if not members <= blank:
                    blank |= members
                    changed = True
    for k in blank:
        out[k] = None
    for k, d in RATE_DENOMINATOR.items():
        if k in out and isinstance(row.get(d), int) and row[d] < n_min:
            out[k] = None
    if any(is_count(k, row[k]) for k in blank):
        for k in out:
            if k.startswith(WHOLE_TABLE):
                out[k] = None
    return out


def suppress_counts(d, n_min=CRITERIA["vs21_min_cell"]):
    """A flat {label: count} for a VS21 report: counts in 1..n_min-1 are blanked; when exactly one is, the next
    smallest positive count is blanked too, so a total elsewhere cannot give the blanked one back."""
    out = dict(d)
    small = [k for k, v in d.items() if is_count(k, v) and 0 < v < n_min]
    for k in small:
        out[k] = None
    if len(small) == 1:
        rest = sorted((v, k) for k, v in d.items() if is_count(k, v) and v > 0 and k not in small)
        if rest:
            out[rest[0][1]] = None
    return out


def write_aggregates(out, result, run_id, rules, local_vs21, extra=None):
    out = Path(out)
    sfx = ".local-vs21" if local_vs21 else ""
    tr = suppress if local_vs21 else (lambda r: r)
    rows = {k: [tr(r) for r in result[k]] for k in ("precision", "competitors", "synthetic", "occupancy")}
    write_rows(out / f"rules_precision{sfx}.csv", rows["precision"], PRECISION_FIELDS)
    write_rows(out / f"rules_competitors{sfx}.csv", rows["competitors"], COMPETITOR_FIELDS)
    write_rows(out / f"rules_synthetic{sfx}.csv", rows["synthetic"], SYNTHETIC_FIELDS)
    write_rows(out / f"rules_occupancy{sfx}.csv", rows["occupancy"], OCCUPANCY_FIELDS)
    gate_rows = {(r["rule_id"], r["pass_id"]): r for r in rows["precision"]
                 if (r["source"], r["instrument"], r["energy"]) == ("pooled", "all", "all")}
    by_rule = []
    for rule in sorted(rules, key=lambda r: r["rule_id"]):
        g = result["gates"][rule["rule_id"]]
        entry = dict(g, parameters={k: rule.get(k) for k in ("product_tolerance", "precursor_tolerance_da", "convention",
                                                              "loss_reference", "semantic_sha256", "rule_file_sha256")})
        if g["status"] != "not_measured":
            entry["passes"] = {p: gate_rows[(rule["rule_id"], p)] for p in PASSES}
            entry["synthetic"] = [r for r in rows["synthetic"] if r["rule_id"] == rule["rule_id"]]
            entry["occupancy"] = [r for r in rows["occupancy"] if r["rule_id"] == rule["rule_id"]]
            entry["competitors"] = [r for r in rows["competitors"] if r["rule_id"] == rule["rule_id"]]
            if local_vs21:
                entry["criteria"] = {k: {kk: vv for kk, vv in v.items() if kk != "value"} for k, v in entry["criteria"].items()}
                entry["by_provenance"] = {b: suppress(v) for b, v in entry["by_provenance"].items()}
        by_rule.append(entry)
    doc = dict(label="ws8-p3-class-rules", run_id=run_id, criteria=CRITERIA, ruleset_sha256=EXPECTED_RULESET_SHA256,
               ruleset_version=EXPECTED_RULESET_VERSION, rules=by_rule,
               as_applied=dict(status="not_done", reason="needs the passing rules and their joint application"),
               status_counts=dict(sorted(Counter(g["status"] for g in result["gates"].values()).items())),
               not_done=NOT_DONE, **(extra or {}))
    p2.write_json(out / f"rules_precision{sfx}.json", doc)
    proposal = dict(label="ws8-p3-review-proposal", run_id=run_id,
                    passed=[dict(rule_id=g["rule_id"], version=g["version"], ion_mode=g["ion_mode"])
                            for g in by_rule if g["status"] == "pass"],
                    incomplete_measured_criteria_met=[dict(rule_id=g["rule_id"], version=g["version"], ion_mode=g["ion_mode"],
                                                           pending=g["pending"]) for g in by_rule if g["status"] == "incomplete"],
                    note="Only `passed` met every pre-registered criterion. `incomplete` rules met every criterion this "
                         "run measures and are NOT eligible for author_reviewed until the pending criteria are measured.")
    p2.write_json(out / f"review_proposal-{run_id}{sfx}.json", proposal)
    return doc


def audit_aggregates(out):
    """No identifier-like value or forbidden column in any aggregate (P2-A5's scan)."""
    for p in sorted(Path(out).iterdir()):
        if p.suffix in (".csv", ".json", ".md") and p.name != "run.json":
            content = p.read_text(encoding="utf-8")
            if re.search(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", content):
                raise AssertionError("Aggregate contains an identifier: " + p.name)
            if p.suffix == ".csv":
                columns = next(csv.reader(content.splitlines()), [])
                if any(re.search(FORBIDDEN_COLUMN, c, re.I) for c in columns):
                    raise AssertionError("Forbidden aggregate column in " + p.name)
    return "pass"


def skip_counts(fire_results):
    """Per (source, set, pass), summed over shards: the console's skipped and error counts and skip reasons."""
    out = {}
    for key, shard in sorted(fire_results.items()):
        source, kind = key.split("/")[:2]
        for p in PASSES:
            c = shard["counts"].get(p)
            if c is None:
                continue
            agg = out.setdefault(f"{source}/{kind}/{p}", Counter())
            agg["skipped"] += c["skipped"]
            agg["errors"] += c["errors"]
            for reason, n in c.get("skipped_by_reason", {}).items():
                agg["skipped:" + reason] += n
    return {k: dict(sorted(v.items())) for k, v in out.items()}


def write_report(out, doc, run, local_vs21):
    sfx = ".local-vs21" if local_vs21 else ""
    sc = suppress_counts if local_vs21 else (lambda d: dict(d))

    def fmt(v):
        return "NA" if v is None else f"{v:.4g}" if isinstance(v, float) else str(v)
    lines = [f"# WS8 P3 rule-level measurement ({'VS20 + VS21' if local_vs21 else 'VS20 development pass'})", "",
             f"Run `{doc['run_id']}`. Rule set {EXPECTED_RULESET_VERSION} (`{EXPECTED_RULESET_SHA256}`). "
             "Library-only measurement; counts use clusters (truth Hill formula, curated precursor type, truth motif "
             "set within a source) unless a column says records.", "",
             "Status: `pass` met every pre-registered criterion; `fail` failed at least one measured criterion; "
             "`not_estimable` has fewer than 100 fired clusters; `incomplete` met every criterion this run measures "
             "while others are pending (the in-source recast set and the VS21 pass are not part of this phase); "
             "`not_measured` is outside measure:true and source_checked:true (or refused by the console).", "",
             "Status counts: " + json.dumps(doc["status_counts"], sort_keys=True), "",
             "Decoy FDR and chimeric precision are cluster-level, as the gate reads them; the record-level values are "
             "in rules_precision.csv and rules_synthetic.csv (report only).", "",
             "| Rule | v | Mode | Primary pass | Fired clusters | Cluster precision | Cluster Wilson lower | Prior stress | Decoy FDR | Occupancy axis / records | Chimeric precision 10/30/50 % | Status | Failed | Pending |",
             "|---|---:|---|---|---:|---:|---:|---:|---:|---|---|---|---|---|"]
    for g in doc["rules"]:
        if g["status"] == "not_measured":
            continue
        row = g["passes"][g["gate_pass"]]
        decoy = g["passes"]["base"]["c_decoy_fdr"]
        occ = next((o for o in g["occupancy"] if o["precondition_mode"] == g["primary_precondition_mode"]), {})
        chim = "/".join(fmt(s["c_precision"]) for s in g["synthetic"]) or "NA"
        lines.append("| " + " | ".join([g["rule_id"], str(g["version"]), g["ion_mode"], g["gate_pass"], fmt(row["c_fires"]),
                                        fmt(row["c_precision"]), fmt(row["c_wilson_lo"]), fmt(row["c_prior_stress"]),
                                        fmt(decoy), fmt(occ.get("axis_occupancy")) + " / " + fmt(occ.get("record_occupancy")),
                                        chim, g["status"], ",".join(g["failed"]) or "-", ",".join(g["pending"]) or "-"]) + " |")
    not_measured = Counter(g["reason"] for g in doc["rules"] if g["status"] == "not_measured")
    lines += ["", "Not measured: " + json.dumps(dict(sorted(not_measured.items())), sort_keys=True), ""]
    prov = sc(run.get("p3_exclusions", {}))
    raw_prov = run.get("p3_exclusions", {})
    lines += ["## Circularity control", ""]
    for source, stats in sorted(run.get("prepared", {}).items()):
        buckets = {k.split(":")[-1]: v for k, v in prov.items() if k.startswith(source + ":bucket:")}
        raw_buckets = {k.split(":")[-1]: v for k, v in raw_prov.items() if k.startswith(source + ":bucket:")}
        total = sum(raw_buckets.values())
        if stats.get("provenance_unavailable"):
            lines.append(f"- {source}: `provenance_unavailable`. No record could be bucketed; no record was excluded on guessed provenance.")
        else:
            lines.append(f"- {source}: provenance buckets (record counts) " + json.dumps(buckets, sort_keys=True) + ".")
        if total and raw_buckets.get("unknown", 0) / total > 0.5:
            lines.append(f"- {source}: `unknown` exceeds 50 %, so label circularity is uncontrolled here and relies on the external set (P3-Q22).")
        insil = raw_prov.get(source + ":in_silico_records", 0)
        shown = prov.get(source + ":in_silico_records", 0)
        if shown is None or any(v is None for v in buckets.values()):
            lines.append(f"- {source}: in-silico detector counts suppressed (VS21 cells below 5).")
        else:
            lines.append(f"- {source}: the declared in-silico detector removed {insil} records "
                         f"({insil / (insil + total):.4%} of the eligible records)." if insil + total else f"- {source}: no eligible records.")
    flagged = [g["rule_id"] for g in doc["rules"] if g.get("provenance_flag")]
    lines += ["- Rules with a provenance bucket of >= 30 fired clusters below 0.95 cluster precision (report only): "
              + (", ".join(flagged) if flagged else "none") + ".", "",
              "## Exclusions and checks", "", "P2-named counters (must equal P2's run.json for the same library on a full run):",
              "", "```json", json.dumps(sc(run.get("exclusions", {})), indent=2, sort_keys=True), "```", "",
              "P2 parity: " + json.dumps(run.get("p2_parity", {}), sort_keys=True), "",
              "P3 exclusions and adduct mapping (records; adduct_unmapped records stay in the universe):", "",
              "```json", json.dumps(prov, indent=2, sort_keys=True), "```", "",
              "Synthetic chimeric set: " + json.dumps({k: sc(v) for k, v in run.get("chimeric", {}).items()}, sort_keys=True), "",
              "Instrument class, P2 bin -> console class, base pass (record counts): "
              + json.dumps(sc(run.get("instrument_classes", {})), sort_keys=True), "",
              "Console skips and errors per source, set and pass (records, summed over shards; `fire` stops on any error):", "",
              "```json", json.dumps({k: sc(v) for k, v in skip_counts(run.get("fire", {})).items()}, indent=1, sort_keys=True),
              "```", "",
              "## Definitions", "",
              "- Universe: every representative of the ion mode; conditional population: representatives where the rule "
              "has at least one hypothesis (the console lists the outcome). Precision is the same in both.",
              "- A cluster fires when any member fires; truth is uniform within a cluster. Cluster-bootstrap intervals "
              f"resample the conditional clusters ({BOOTSTRAP} replicates by default, seeded) and give record- and cluster-level precision.",
              "- Decoy FDR: base-pass conditional population; mean over +5, +9 and -7 Da of decoy fires / real fires, "
              "in clusters for the gate and in records (report only).",
              "- Occupancy: (axis) the fraction of the nominal-mass span of the rule's hypothesis precursors within the rule's "
              "precursor tolerance of a hypothesis; (records) the fraction of motif-negative representatives of the mode "
              "whose precondition holds. Either above 5 % makes the rule `non_conditional`.",
              "- Prior stress: P2's statistics() on the cluster universe table (library-wide prevalence of the mode, halved, "
              "LR+ lower bound). The conditional stress grid is report-only.",
              "- Per-source precision uses library sources (VS20, VS21); provenance buckets are reported only.",
              "", "## Not done in this phase", ""] + [f"- {x}" for x in NOT_DONE]
    (Path(out) / f"report{sfx}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------------------------- measure stage
def load_rules(ruleset_dir):
    res = rf.load_files(rf.read_directory(ruleset_dir))
    if not res.ok:
        raise RuntimeError("Rule set does not load: " + str(res.errors[:3]))
    if res.ruleset_sha256 != EXPECTED_RULESET_SHA256 or res.rule_set_version != EXPECTED_RULESET_VERSION:
        raise RuntimeError("Rule set is not the pinned 0.1.0 set")
    raw = ruleset_file_hashes(ruleset_dir)
    rules = []
    for rid, info in sorted(res.rules.items()):
        d = info.doc
        pre = d["preconditions"]
        tol = round(pre["precursor_tolerance_da"] * rf.PDA_PER_DA)
        hyps = [(h.precursor_pda, h.adduct_source) for h, _ in rf.iter_hypotheses(info, res.tables)]
        rules.append(dict(rule_id=rid, version=d["version"], ion_mode=pre["ion_mode"], motif_id=d["claim"]["motif_id"],
                          competitors=[c["motif_id"] for c in pre["competitors"]], measure=d["measure"],
                          source_checked=d["source_checked"], primary_mode=d["primary_precondition_mode"],
                          origin=d["parameters_origin"], product_tolerance=d["product_tolerance"],
                          precursor_tolerance_da=pre["precursor_tolerance_da"], convention=d["convention"],
                          loss_reference=d["loss_reference"], semantic_sha256=res.semantic_sha256(rid),
                          rule_file_sha256=raw.get(info.file),
                          axis={"adduct_ambiguous": axis_occupancy([p for p, _ in hyps], tol),
                                "mz_only": axis_occupancy([p for p, s in hyps if s == "rule"], tol)}))
    return rules, res, raw


def load_truth(out, source):
    import pyarrow.parquet as pq
    cols = ["record_id", "source", "mode", "motifs", "cluster_id", "instrument", "energy", "bucket"]
    t = pq.read_table(p2.row_path(Path(out), "prepared", source), columns=cols).to_pylist()
    return {(source, r["record_id"]): dict(mode=r["mode"], motifs=frozenset(r["motifs"]), cluster=r["cluster_id"],
                                           instrument=r["instrument"], energy=r["energy"], bucket=r["bucket"], source=source)
            for r in t}


def load_set(out, source, kind):
    """Evidence of one set (representatives or chimeric) joined to its index."""
    msp_dir = Path(out) / "msp" / source
    index = read_rows(msp_dir / f"{kind}-index{suffix(source)}.tsv")
    by_file = defaultdict(list)
    for r in index:
        uid = (source, int(r["record_id"])) + ((float(r["level"]),) if kind == "chimeric" else ())
        by_file[r["file_name"]].append((int(r["record_index"]), uid))
    evidence = []
    for name, ids in sorted(by_file.items()):
        ids.sort()
        if [i for i, _ in ids] != list(range(len(ids))):
            raise RuntimeError("Index rows are not contiguous for " + name)
        passes = PASSES if kind == "representatives" else ("base", "mz_only_curated")
        evidence.append(read_evidence(Path(out) / "fire" / source / kind / Path(name).stem, name, [u for _, u in ids], passes))
    return evidence, [u for ids in by_file.values() for _, u in ids]


def write_fires(sink, by_rule, truth, rules, kind, level=None):
    motif = {r["rule_id"]: r["motif_id"] for r in rules}
    for p in sorted(by_rule):
        for rid in sorted(by_rule[p]):
            for uid, fired, nh, veto, adduct in by_rule[p][rid]:
                sink.add(dict(set=kind, level=level, source=uid[0], record_id=uid[1], pass_id=p, rule_id=rid,
                              n_hypotheses=nh, fired=fired, veto_blocked=veto, adduct=adduct,
                              motif_match=motif.get(rid) in truth[uid[:2]]["motifs"]))


def measure(args, run):
    import pyarrow as pa
    rules, res, raw = load_rules(args.ruleset_dir)
    run["ruleset"] = dict(path=str(args.ruleset_dir), sha256=res.ruleset_sha256, version=res.rule_set_version,
                          files_raw_sha256=raw, files_canonical_sha256=dict(sorted(res.file_hashes.items())))
    versions = {r["rule_id"]: r["version"] for r in rules}
    truth, main_ev, chim_ev, chim_universe = {}, [], [], defaultdict(list)
    for source in args.sources:
        if source == "vs21":
            vs21_guard(args)
        truth.update(load_truth(args.out, source))
        ev, _ = load_set(args.out, source, "representatives")
        main_ev += ev
        ev, uids = load_set(args.out, source, "chimeric")
        chim_ev += ev
        for u in uids:
            chim_universe[u[2]].append(u)
    by_rule, notev, inst = collect(main_ev, versions)
    chim_by_level, chim_notev = {}, {}
    for level in CHIMERIC_LEVELS:
        sel = [{p: [rec for rec in recs if rec[0][2] == level] for p, recs in ev.items()} for ev in chim_ev]
        chim_by_level[level], chim_notev[level], _ = collect(sel, versions)
    refused = sorted({r["rule_id"] for shard in run.get("fire", {}).values() for r in shard.get("refusals", {}).get("base", [])})
    pairs = Counter(f"{truth[u[:2]]['instrument']} -> {c}" for u, c in inst.items())
    run["instrument_classes"] = dict(sorted(pairs.items()))
    result = measure_core(rules, truth, by_rule, notev, args.sources, args.seed, args.bootstrap,
                          dict(by_rule=chim_by_level, universe=chim_universe, notev=chim_notev), refused)
    local = any_vs21(args.sources)
    schema = pa.schema([("set", pa.string()), ("level", pa.float64()), ("source", pa.string()), ("record_id", pa.int64()),
                        ("pass_id", pa.string()), ("rule_id", pa.string()), ("n_hypotheses", pa.int64()),
                        ("fired", pa.bool_()), ("veto_blocked", pa.bool_()), ("adduct", pa.string()), ("motif_match", pa.bool_())])
    sink = p2.ParquetRows(args.out / ("fires" + (".local-vs21" if local else "") + ".parquet"), schema)
    write_fires(sink, by_rule, truth, rules, "representatives")
    for level in CHIMERIC_LEVELS:
        write_fires(sink, chim_by_level[level], truth, rules, "chimeric", level)
    sink.close()
    doc = write_aggregates(args.out, result, args.out.name, rules, local,
                           extra=dict(refused_rule_ids=refused, sources=args.sources))
    write_report(args.out, doc, run, local)
    run["acceptance_artifact_audit"] = audit_aggregates(args.out)
    run["status_counts"] = doc["status_counts"]


# ------------------------------------------------------------------------------------------- VS21 and runs
FREEZE_KEYS = ("ruleset_sha256_expected", "ruleset_files", "console_sha256", "console_files", "code", "parameters")


def freeze_hashes(ruleset_dir, console, seed, bootstrap):
    return dict(ruleset_sha256_expected=EXPECTED_RULESET_SHA256,
                ruleset_files=ruleset_file_hashes(ruleset_dir), console_sha256=p2.sha256(console),
                console_files=console_files(console),
                code={p.relative_to(ROOT).as_posix(): p2.sha256(p) for p in CODE_FILES},
                parameters=dict(seed=seed, bootstrap=bootstrap, console_ini=CONSOLE_INI))


def vs21_guard(args, out_root=OUT_ROOT):
    """The VS21 path runs only with --vs21 and a freeze file (named and placed as `freeze` writes it) whose hashes
    match the rule set, the console's bin folder, the measurement code and the parameters (seed, bootstrap, console
    INI) now. It never opens a VS21 file."""
    if not getattr(args, "vs21", False):
        raise PermissionError("VS21 needs the explicit --vs21 flag")
    freeze = getattr(args, "freeze", None)
    if not freeze or not Path(freeze).is_file():
        raise PermissionError("VS21 needs a freeze file (--freeze)")
    if not getattr(args, "ruleset_dir", None) or not getattr(args, "console", None):
        raise PermissionError("VS21 needs --ruleset-dir and --console to check the freeze")
    try:
        check_freeze_path(freeze, out_root)
    except ValueError as exc:
        raise PermissionError(str(exc)) from None
    frozen = json.loads(Path(freeze).read_text(encoding="utf-8"))
    now = freeze_hashes(args.ruleset_dir, args.console, getattr(args, "seed", None), getattr(args, "bootstrap", None))
    diffs = [k for k in FREEZE_KEYS if frozen.get(k) != now[k]]
    if not frozen.get("console_identity"):
        diffs.append("console_identity")
    if diffs:
        raise PermissionError("Freeze mismatch, VS21 refused: " + ", ".join(diffs))
    return True


def check_out(out, limit, commit, out_root=OUT_ROOT):
    """Run directory: p3cs-<yyyymmdd>-<sha> (or p3cssmoke- with --limit) directly under the ws8-p3 root, the sha
    a prefix of the code commit."""
    out = Path(out).resolve()
    if out.parent != Path(out_root).resolve():
        raise ValueError("--out must be directly under " + Path(out_root).as_posix())
    m = re.fullmatch(RUN_NAME, out.name)
    if not m:
        raise ValueError("--out must be named p3cs-<yyyymmdd>-<sha> or p3cssmoke-<yyyymmdd>-<sha>")
    if bool(m.group(1)) != (limit is not None):
        raise ValueError("A smoke run (p3cssmoke-) needs --limit and a full run (p3cs-) must not have it")
    if not commit.startswith(m.group(3)):
        raise ValueError("The run id's sha must be a prefix of the code commit")
    return out


def check_freeze_path(path, out_root=OUT_ROOT):
    path = Path(path).resolve()
    if path.parent != Path(out_root).resolve() or not re.fullmatch(FREEZE_NAME, path.name):
        raise ValueError("--freeze must be freeze-<yyyymmdd>-<sha>.json directly under " + Path(out_root).as_posix())
    return path


def check_run_state(run, args):
    if run.get("status") == "complete":
        raise ValueError("Completed run is immutable; use a new code commit/run ID")
    if "parameters" in run:
        for field in ("sources", "workers", "limit", "seed", "bootstrap"):
            if run["parameters"].get(field) != getattr(args, field):
                raise ValueError("Stage arguments must match the run manifest: " + field)


def detach(argv, out):
    """P2's WMI launch (p2_diagnostics.main): the child is created by the WMI service, outside the app's job."""
    command = [sys.executable, "-B", str(Path(__file__).resolve())] + [x for x in argv if x != "--detach"] + ["--detached-child"]
    env = dict(os.environ, WS8_P3_COMMAND=subprocess.list2cmdline(command), WS8_P3_CWD=str(ROOT))
    script = ("$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$env:WS8_P3_COMMAND; "
              "CurrentDirectory=$env:WS8_P3_CWD}; if ($r.ReturnValue -ne 0) { throw 'WMI detached launch failed' }; Write-Output $r.ProcessId")
    return subprocess.check_output(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", script],
                                   env=env, text=True, creationflags=NO_WINDOW).strip()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("stage", choices=("prepare", "fire", "measure", "all", "freeze"))
    parser.add_argument("--out", type=Path, help="run directory p3cs[smoke]-<yyyymmdd>-<sha> under " + OUT_ROOT.as_posix())
    parser.add_argument("--sources", default="vs20")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=float, help="deterministic per-file fraction (0,1]; 0.01 for the smoke")
    parser.add_argument("--ruleset-dir", type=Path, help="the RuleSet folder of the console's source (0.1.0)")
    parser.add_argument("--console", type=Path, help="MsfinderConsoleApp.exe")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--bootstrap", type=int, default=BOOTSTRAP)
    parser.add_argument("--vs21", action="store_true", help="allow the VS21 source (needs a matching --freeze)")
    parser.add_argument("--freeze", type=Path, help="freeze file to write (stage freeze) or to check (VS21)")
    parser.add_argument("--detach", action="store_true", help="launch outside the caller's Windows job (WMI), log to run.log")
    parser.add_argument("--detached-child", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    args.sources = args.sources.split(",")
    if args.sources not in (["vs20"], ["vs20", "vs21"]) or not 1 <= args.workers <= 8:
        # VS21 needs VS20 in the same run: its records are deduplicated against the VS20 first blocks.
        parser.error("Use --sources vs20 or vs20,vs21 and 1..8 workers")
    if args.limit is not None and not 0 < args.limit <= 1:
        parser.error("--limit must be a fraction in (0,1]")
    if args.stage in ("prepare", "measure", "all", "freeze") and not args.ruleset_dir:
        parser.error("--ruleset-dir is required for " + args.stage)
    if args.stage in ("fire", "all", "freeze") and not args.console:
        parser.error("--console is required for " + args.stage)
    commit = git("rev-parse", "HEAD")
    if args.stage == "freeze":
        path = check_freeze_path(args.freeze)
        if path.exists():
            parser.error("Freeze file exists; a freeze is never rewritten")
        p2.write_json(path, dict(label="ws8-p3-freeze", code_commit=commit, console_identity=console_identity(args.console),
                                 **freeze_hashes(args.ruleset_dir, args.console, args.seed, args.bootstrap)))
        print("freeze written", path.name, flush=True)
        return
    if any_vs21(args.sources):
        try:
            vs21_guard(args)
        except PermissionError as exc:
            parser.error(str(exc))
    if not args.out:
        parser.error("--out is required")
    try:
        args.out = check_out(args.out, args.limit, commit)
    except ValueError as exc:
        parser.error(str(exc))
    if git("status", "--porcelain", "--untracked-files=no"):
        parser.error("Run from a clean, committed pin")
    if not git("branch", "-r", "--contains", "HEAD"):
        parser.error("Run from a pushed commit")
    args.out.mkdir(parents=True, exist_ok=True)
    if args.detach:
        if sys.platform != "win32":
            parser.error("Detached launch requires Windows")
        print("Detached P3 process:", detach(argv, args.out), flush=True)
        return
    if args.detached_child:
        sys.stdout = open(args.out / "run.log", "a", encoding="utf-8", buffering=1)
        sys.stderr = sys.stdout
    run_path = args.out / "run.json"
    run = json.loads(run_path.read_text(encoding="utf-8")) if run_path.exists() else {}
    try:
        check_run_state(run, args)
    except ValueError as exc:
        parser.error(str(exc))
    priority_name = p2.priority()
    context = p2.process_context()
    if args.detached_child and context["launcher"]["in_windows_job"]:
        raise RuntimeError("Detached launcher still belongs to a Windows job")
    import pyarrow
    versions = dict(python=sys.version, pyarrow=pyarrow.__version__)
    if args.stage in ("prepare", "all"):
        import rdkit
        versions["rdkit"] = rdkit.__version__
    run.update(label="ws8-p3-class-rules", run_id=args.out.name, code_commit=commit, pin_path=str(ROOT),
               interpreter=sys.executable, versions=dict(run.get("versions", {}), **versions),
               code_sha256={p.relative_to(ROOT).as_posix(): p2.sha256(p) for p in CODE_FILES},
               parameters=dict(sources=args.sources, workers=args.workers, limit=args.limit, seed=args.seed,
                               bootstrap=args.bootstrap, priority=priority_name, criteria=CRITERIA,
                               threads={v: os.environ[v] for v in p2.THREAD_VARS}, passes=list(PASSES),
                               primary_pass=PRIMARY_PASS, chimeric=dict(levels=list(CHIMERIC_LEVELS), window_da=CHIMERIC_WINDOW_DA,
                                                                        partner="same mode, other first block, seeded"),
                               console_ini=CONSOLE_INI, record_unit="one representative per first block, ion mode and "
                               "curated precursor type (alias-mapped when the adduct table maps it)",
                               cluster="(truth Hill formula, curated precursor type, truth motif set) within a source",
                               expected_ruleset_sha256=EXPECTED_RULESET_SHA256),
               process=dict(detached_child=args.detached_child, launch_method="WMI" if args.detached_child else "foreground", **context),
               status="running", not_done=NOT_DONE)
    run.setdefault("timings_seconds", {})
    p2.write_json(run_path, run)
    stages = ("prepare", "fire", "measure") if args.stage == "all" else (args.stage,)
    try:
        for stage in stages:
            start = time.perf_counter()
            globals()[stage](args, run)
            run["timings_seconds"][stage] = time.perf_counter() - start
            p2.write_json(run_path, run)
        if args.stage in ("measure", "all"):
            run["status"] = "complete"
        run["files"] = {str(p.relative_to(args.out)).replace("\\", "/"): p2.sha256(p) for p in sorted(args.out.rglob("*"))
                        if p.is_file() and p.name not in ("run.json", "run.log")}
        p2.write_json(run_path, run)
        print("completed", args.stage, run["timings_seconds"], flush=True)
    except Exception as exc:
        run.update(status="failed", error_type=type(exc).__name__)
        p2.write_json(run_path, run)
        # Never expose an exception text that may carry library data or paths in a shared log.
        print("FAILED", type(exc).__name__, flush=True)
        import traceback
        for frame in traceback.extract_tb(exc.__traceback__):
            print(Path(frame.filename).name, frame.lineno, frame.name, flush=True)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
