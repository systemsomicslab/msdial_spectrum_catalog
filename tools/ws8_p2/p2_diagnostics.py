"""WS8 P2 library-only diagnostic calibration; optional chemistry stays outside src/."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import re
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from tools.ws8_p1.mcs_core import standardise, first_block, hill_formula, RejectedStructure

THREAD_VARS = ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
for _var in THREAD_VARS:
    os.environ[_var] = "1"
OUT_ROOT = Path("F:/MSDIAL_Public_Reanalysis/annotation/ws8-p2")
LIB_ROOT = Path("D:/13_MSDIAL_Public_Reanalysis/libraries")
P1_MOTIFS = ROOT / "tools/ws8_p1/motifs-p1.tsv"
P1_RESULT = Path("F:/MSDIAL_Public_Reanalysis/annotation/ws8-p1/p1-20261113-8153661/motif_precision.csv")
PUBLIC = [("21901200/MSMS-Public_all-pos-VS20.msp", "16c0e3513f9328215e8dcde5c44e2c4ffda93b7b8ed76b501797a2df5e8751e7"),
          ("21904103/MSMS-Public_all-neg-VS20.msp", "556091ff7e902907b3df055c04ec588fa544d5c2a56adf34f365444ccd13a1eb")]
RULE = dict(n_min=100, wilson_target=0.95, source_n_min=30, source_precision=0.95,
            prior_fraction=0.5, prior_target=0.90, interval="wilson95", lr_interval="log95",
            zero_correction="add 0.5 to all four cells if any cell is zero")
VARIANTS = {"default": (0.01, 0.01), "tol_0.005": (0.005, 0.01),
            "tol_0.02": (0.02, 0.01), "intensity_0.05": (0.01, 0.05),
            "intensity_0.10": (0.01, 0.10)}
Z = 1.959963984540054


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_csv(path, rows, fields=None):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def formula_counts(formula):
    formula = re.sub(r"[+-]\d*$", "", formula)  # RDKit Hill formula of permanent ions.
    parts = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    if "".join(e + n for e, n in parts) != formula:
        raise ValueError("Invalid elemental formula")
    return {e: int(n or 1) for e, n in parts}


def formula_mass(formula, kind, mode):
    from rdkit import Chem
    table = Chem.GetPeriodicTable()
    mass = sum(table.GetMostCommonIsotopeMass(e) * n for e, n in formula_counts(formula).items())
    return mass + (0 if kind == "neutral_loss" else -0.000549 if mode == "pos" else 0.000549)


def diagnostics():
    rows = read_tsv(HERE / "diagnostics-p2.tsv")
    for r in rows:
        for field in ("mz_or_loss", "tolerance_da", "min_rel_intensity"):
            r[field] = float(r[field])
        r["precursor_types"] = r["precursor_types"].split(";")
        if not r["source"] or abs(r["mz_or_loss"] - formula_mass(r["formula"], r["kind"], r["ion_mode"])) > 1e-5:
            raise ValueError("Diagnostic table fails source or formula validation")
    if len({r["diagnostic_id"] for r in rows}) != len(rows):
        raise ValueError("Diagnostic IDs must be unique")
    return rows


def motif_rows():
    rows = read_tsv(P1_MOTIFS) + read_tsv(HERE / "motifs-p2-extra.tsv")
    if len({r["motif_id"] for r in rows}) != len(rows):
        raise ValueError("Motif IDs collide")
    return rows


def wilson(tp, n):
    if n == 0:
        return 0.0, 1.0
    p, z2 = tp / n, Z * Z
    mid = (p + z2 / (2 * n)) / (1 + z2 / n)
    half = Z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    return max(0, mid - half), min(1, mid + half)


def statistics(tp, fp, fn, tn):
    """Uncorrected point rates; corrected LR+ and log interval when a cell is zero."""
    n, obs, positive, negative = tp + fp + fn + tn, tp + fp, tp + fn, fp + tn
    lo, hi = wilson(tp, obs)
    recall = tp / positive if positive else None
    fpr = fp / negative if negative else None
    lr = lrlo = lrhi = stress = None
    if positive and negative:
        a, b, c, d = (x + 0.5 for x in (tp, fp, fn, tn)) if min(tp, fp, fn, tn) == 0 else (tp, fp, fn, tn)
        lr = (a / (a + c)) / (b / (b + d))
        se = math.sqrt(max(0, 1 / a - 1 / (a + c) + 1 / b - 1 / (b + d)))
        lrlo, lrhi = lr * math.exp(-Z * se), lr * math.exp(Z * se)
        prior = 0.5 * positive / n
        stress = prior * lrlo / (1 - prior + prior * lrlo)
    return dict(tp=tp, fp=fp, fn=fn, tn=tn, n=n, n_obs=obs,
                precision=tp / obs if obs else None, wilson_lo=lo, wilson_hi=hi,
                recall=recall, fpr=fpr, lr_plus=lr, lr_lo=lrlo, lr_hi=lrhi,
                prevalence=positive / n if n else None, prior_stress_precision=stress)


def status(row, sources=(), n_min=100):
    if row["n_obs"] < n_min:
        return "not_estimable"
    return "pass" if (row["wilson_lo"] >= 0.95 and row["prior_stress_precision"] is not None
                      and row["prior_stress_precision"] >= 0.90
                      and all(s["n_obs"] < 30 or s["precision"] >= 0.95 for s in sources)) else "fail"


def nominal_mass(peaks, decimals):
    return bool(peaks) and (all(abs(m - round(m)) <= 0.001 + 1e-10 for m, _ in peaks)
                            or all(d < 2 for d in decimals))


def observe_one(peaks, precursor, adduct, diagnostic, tolerance=None, floor=None):
    if adduct not in diagnostic["precursor_types"]:
        return False, 0.0
    base = max((i for m, i in peaks if abs(m - precursor) > 1.0), default=0)
    if base <= 0:
        return False, 0.0
    loss = diagnostic["kind"] == "neutral_loss"
    target = precursor - diagnostic["mz_or_loss"] if loss else diagnostic["mz_or_loss"]
    tol = diagnostic["tolerance_da"] if tolerance is None else tolerance
    minimum = diagnostic["min_rel_intensity"] if floor is None else floor
    rel = max((i / base for m, i in peaks if abs(m - target) <= tol + 1e-10
               and (not loss or precursor - m >= 1.0)), default=0.0)
    return rel >= minimum and rel > 0, rel


def formula_gate(truth_formula, diagnostic_formula):
    """Require residue heavy-element counts; hydrogen depends on charge/adduct."""
    truth = formula_counts(truth_formula)
    return all(truth.get(e, 0) >= n for e, n in formula_counts(diagnostic_formula).items() if e != "H")


def instrument_bin(value):
    s = value.upper()
    return "QTOF-type" if "TOF" in s else "Orbitrap/FT" if any(t in s for t in ("ORBITRAP", "FT", "Q EXACTIVE", "Q-EXACTIVE")) else "other"


def ion_mode(value, fallback):
    """MSP libraries use both long labels and the P/N abbreviations."""
    s = value.strip().lower()
    return {"p": "pos", "positive": "pos", "pos": "pos", "n": "neg", "negative": "neg", "neg": "neg", "": fallback}.get(s)


def energy_bin(value):
    """Single explicit energies only. Keep normalized energy apart from eV."""
    s = value.strip().lower()
    m = re.fullmatch(r"(?:nce\s*[=:]?\s*)?(\d+(?:\.\d+)?)\s*(ev|v|%|nce)?", s)
    if not m:
        return "unparseable"
    unit = "NCE" if "%" in s or "nce" in s else "eV" if m[2] in ("ev", "v") else "unspecified"
    n = float(m[1])
    return unit + (":<20" if n < 20 else ":20-<40" if n < 40 else ":>=40")


def priority():
    import ctypes
    if sys.platform != "win32":
        raise RuntimeError("Real runs require Windows BelowNormal priority")
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.GetCurrentProcess.restype = ctypes.c_void_p
    k.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
    k.GetPriorityClass.argtypes = (ctypes.c_void_p,)
    handle = k.GetCurrentProcess()
    if not k.SetPriorityClass(handle, 0x4000) or k.GetPriorityClass(handle) != 0x4000:
        raise RuntimeError("BelowNormal verification failed")
    return "BelowNormal"


def process_context():
    """Inspect our own launcher, not other jobs or campaign processes.

    Windows' venv launcher puts the real interpreter in its own job. Verify the
    launcher is outside a job rather than incorrectly rejecting that venv job.
    """
    import ctypes
    from ctypes import wintypes
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    nt = ctypes.WinDLL("ntdll")
    class BasicInfo(ctypes.Structure):
        _fields_ = [("exit", ctypes.c_void_p), ("peb", ctypes.c_void_p),
                    ("affinity", ctypes.c_void_p), ("base_priority", ctypes.c_void_p),
                    ("pid", ctypes.c_size_t), ("parent", ctypes.c_size_t)]
    k.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    k.OpenProcess.restype = wintypes.HANDLE
    k.CloseHandle.argtypes = (wintypes.HANDLE,)
    k.IsProcessInJob.argtypes = (wintypes.HANDLE, wintypes.HANDLE, ctypes.POINTER(wintypes.BOOL))
    k.QueryFullProcessImageNameW.argtypes = (wintypes.HANDLE, wintypes.DWORD, wintypes.LPWSTR, ctypes.POINTER(wintypes.DWORD))
    nt.NtQueryInformationProcess.argtypes = (wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.ULONG, ctypes.c_void_p)
    def info(pid):
        handle = k.OpenProcess(0x400 | 0x1000, False, pid)
        if not handle:
            raise RuntimeError("Cannot inspect own process ancestry")
        try:
            basic, inside = BasicInfo(), wintypes.BOOL()
            buf, size = ctypes.create_unicode_buffer(32768), wintypes.DWORD(32768)
            if nt.NtQueryInformationProcess(handle, 0, ctypes.byref(basic), ctypes.sizeof(basic), None) != 0:
                raise RuntimeError("Cannot inspect own parent process")
            if not k.IsProcessInJob(handle, None, ctypes.byref(inside)) or not k.QueryFullProcessImageNameW(handle, 0, buf, ctypes.byref(size)):
                raise RuntimeError("Cannot verify own process context")
            return dict(pid=pid, parent_pid=basic.parent, image=Path(buf.value).name, in_windows_job=bool(inside.value))
        finally:
            k.CloseHandle(handle)
    current = info(os.getpid())
    launcher = info(current["parent_pid"])
    return dict(current=current, launcher=launcher)


QUERIES = None


def worker_init():
    global QUERIES
    priority()
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    QUERIES = [(r["motif_id"], Chem.MolFromSmarts(r["smarts"])) for r in motif_rows()]
    if any(q is None for _, q in QUERIES):
        raise ValueError("Invalid motif SMARTS")


def chemistry(smiles):
    try:
        mol = standardise(smiles)
        return dict(key=first_block(mol), formula=hill_formula(mol),
                    motifs=[mid for mid, q in QUERIES if mol.HasSubstructMatch(q, useChirality=False)])
    except Exception as exc:
        return dict(rejected=str(exc) if isinstance(exc, RejectedStructure) else "standardisation")


def iter_msp(path):
    """Allowlist metadata and retain lexical m/z precision for the resolution gate."""
    allowed = {"precursormz", "precursortype", "ionmode", "instrumenttype", "collisionenergy", "smiles"}
    meta, peaks, decimals, in_peaks = {}, [], [], False
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                if meta:
                    yield meta, peaks, decimals
                meta, peaks, decimals, in_peaks = {}, [], [], False
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
                if k in ("num peaks", "numpeaks"):
                    in_peaks = True
        if meta:
            yield meta, peaks, decimals


class ParquetRows:
    def __init__(self, path, schema):
        import pyarrow.parquet as pq
        self.writer = pq.ParquetWriter(path, schema, compression="zstd")
        self.schema, self.rows = schema, []

    def add(self, row):
        self.rows.append(row)
        if len(self.rows) >= 4000:
            self.flush()

    def flush(self):
        import pyarrow as pa
        if self.rows:
            self.writer.write_table(pa.Table.from_pylist(self.rows, schema=self.schema))
            self.rows.clear()

    def close(self):
        self.flush()
        self.writer.close()


def parquet_rows(path):
    import pyarrow.parquet as pq
    for batch in pq.ParquetFile(path).iter_batches(batch_size=4000):
        yield from batch.to_pylist()


def row_path(out, stem, source):
    return out / (stem + (".local-vs21" if source == "vs21" else "") + ".parquet")


def prepare(args, run):
    import pyarrow as pa
    schema = pa.schema([("record_id", pa.int64()), ("raw_smiles", pa.string()), ("mode", pa.string()),
                        ("adduct", pa.string()), ("precursor", pa.float64()), ("instrument", pa.string()),
                        ("energy", pa.string()), ("mz", pa.list_(pa.float64())),
                        ("intensity", pa.list_(pa.float64()))])
    prepared_schema = pa.schema([f for f in schema if f.name != "raw_smiles"] + [
        ("key", pa.string()), ("formula", pa.string()), ("motifs", pa.list_(pa.string())),
        ("representative", pa.bool_()), ("spectrum_count", pa.int64())])
    public_keys = set()
    counts = Counter()
    run["inputs"] = [dict(name=p.name, sha256=sha256(p)) for p in (P1_MOTIFS, HERE / "motifs-p2-extra.tsv", HERE / "diagnostics-p2.tsv")]
    run["prepared"] = {}
    for source in ("vs20", "vs21"):
        if source == "vs21" and source not in args.sources:
            continue
        paths = [LIB_ROOT / p for p, _ in PUBLIC] if source == "vs20" else [
            next(LIB_ROOT.rglob(f"MSMS-Private-{mode}-VS21.msp")) for mode in ("pos", "neg")]
        raw_path = row_path(args.out, "raw", source)
        sink = ParquetRows(raw_path, schema)
        raw_smiles = set()
        for mode, path in zip(("pos", "neg"), paths):
            digest = sha256(path)
            if source == "vs20" and digest != dict(PUBLIC)[str(path.relative_to(LIB_ROOT)).replace("\\", "/")]:
                raise RuntimeError("Public input checksum mismatch")
            run["inputs"].append(dict(name=path.name, sha256=digest))
            for ordinal, (meta, peaks, decimals) in enumerate(iter_msp(path)):
                prefix = source + ":" + mode + ":"
                counts[prefix + "input"] += 1
                raw_smiles.add(meta.get("smiles", ""))
                # Sampling never limits the VS20 identity universe used to remove VS21 overlaps.
                if args.limit and ordinal % max(1, round(1 / args.limit)) != 0:
                    counts[prefix + "smoke_unsampled"] += 1
                    continue
                if not peaks:
                    counts[prefix + "no_peaks"] += 1
                    continue
                if nominal_mass(peaks, decimals):
                    counts[prefix + "nominal_mass"] += 1
                    continue
                try:
                    precursor = float(meta.get("precursormz", ""))
                    if not math.isfinite(precursor) or precursor <= 0:
                        raise ValueError
                except ValueError:
                    counts[prefix + "precursor"] += 1
                    continue
                actual_mode = ion_mode(meta.get("ionmode", ""), mode)
                if actual_mode != mode:
                    counts[prefix + "mode_mismatch"] += 1
                    continue
                adduct = meta.get("precursortype", "").replace(" ", "")
                if not adduct:
                    counts[prefix + "missing_precursor_type"] += 1
                    continue
                sink.add(dict(record_id=ordinal * 2 + (mode == "neg"), raw_smiles=meta.get("smiles", ""),
                              mode=mode, adduct=adduct, precursor=precursor,
                              instrument=instrument_bin(meta.get("instrumenttype", "")),
                              energy=energy_bin(meta.get("collisionenergy", "")),
                              mz=[p[0] for p in peaks], intensity=[p[1] for p in peaks]))
        sink.close()
        smiles_order = sorted(raw_smiles)
        with mp.get_context("spawn").Pool(args.workers, initializer=worker_init) as pool:
            chem = dict(zip(smiles_order, pool.imap(chemistry, smiles_order, chunksize=64)))
        for result in chem.values():
            if "rejected" in result:
                counts[source + ":unique_smiles:" + result["rejected"]] += 1
        if source == "vs20":
            public_keys = {r["key"] for r in chem.values() if "key" in r}
            run["vs20_dedup_universe"] = len(public_keys)
        reps, ns = {}, Counter()
        for row in parquet_rows(raw_path):
            c = chem[row["raw_smiles"]]
            if "rejected" in c:
                counts[source + ":structure_rejected_records"] += 1
                continue
            if source == "vs21" and c["key"] in public_keys:
                counts[source + ":vs20_overlap_records"] += 1
                continue
            unit = (c["key"], row["mode"], row["adduct"])
            ns[unit] += 1
            candidate = (len(row["mz"]), -row["record_id"])
            if unit not in reps or candidate > reps[unit]:
                reps[unit] = candidate
        sink = ParquetRows(row_path(args.out, "prepared", source), prepared_schema)
        kept = 0
        for row in parquet_rows(raw_path):
            c = chem[row.pop("raw_smiles")]
            if "key" not in c or (source == "vs21" and c["key"] in public_keys):
                continue
            unit = (c["key"], row["mode"], row["adduct"])
            sink.add(dict(row, **c, representative=-reps[unit][1] == row["record_id"], spectrum_count=ns[unit]))
            kept += 1
        sink.close()
        raw_path.unlink()
        run["prepared"][source] = dict(spectra=kept, representatives=len(reps), compounds=len({u[0] for u in reps}))
        print("prepare", source, run["prepared"][source], flush=True)
    run["exclusions"] = dict(sorted(counts.items()))


def add_confusion(counter, key, truth, observed, weight=1.0):
    cells = counter.setdefault(key, [0.0, 0.0, 0.0, 0.0])
    cells[0 if truth and observed else 1 if observed else 2 if truth else 3] += weight


def observe(args, run):
    import pyarrow as pa
    ds = diagnostics()
    counters = {}
    schema = pa.schema([("key", pa.string()), ("record_id", pa.int64()), ("diagnostic_id", pa.string()),
                        ("ion_mode", pa.string()), ("precursor_type", pa.string()), ("observed", pa.bool_()),
                        ("intensity", pa.float64()), ("motif_match", pa.bool_())])
    for source in args.sources:
        sink = ParquetRows(row_path(args.out, "observations", source), schema)
        n = 0
        for r in parquet_rows(row_path(args.out, "prepared", source)):
            peaks = list(zip(r["mz"], r["intensity"]))
            for d in ds:
                if d["ion_mode"] != r["mode"] or r["adduct"] not in d["precursor_types"]:
                    continue
                truth = d["motif_id"] in r["motifs"]
                observed, intensity = observe_one(peaks, r["precursor"], r["adduct"], d)
                common = (d["diagnostic_id"], d["ion_mode"], source)
                # All spectra weighted to total one per compound/mode/precursor type.
                for adduct in ("all", r["adduct"]):
                    add_confusion(counters, common + (adduct, "all_spectra_weighted", "all"), truth, observed, 1 / r["spectrum_count"])
                if not r["representative"]:
                    continue
                sink.add(dict(key=r["key"], record_id=r["record_id"], diagnostic_id=d["diagnostic_id"],
                              ion_mode=r["mode"], precursor_type=r["adduct"], observed=observed,
                              intensity=intensity, motif_match=truth))
                n += 1
                for variant, (tol, floor) in VARIANTS.items():
                    obs = observed if variant == "default" else observe_one(peaks, r["precursor"], r["adduct"], d, tol, floor)[0]
                    for adduct in ("all", r["adduct"]):
                        add_confusion(counters, common + (adduct, variant, "all"), truth, obs)
                if formula_gate(r["formula"], d["formula"]):
                    for adduct in ("all", r["adduct"]):
                        add_confusion(counters, common + (adduct, "formula_gated", "all"), truth, observed)
                add_confusion(counters, common + ("all", "instrument", r["instrument"]), truth, observed)
                add_confusion(counters, common + ("all", "collision_energy", r["energy"]), truth, observed)
        sink.close()
        print("observe", source, n, "representative diagnostic rows", flush=True)
    rows = [dict(diagnostic_id=k[0], ion_mode=k[1], source=k[2], precursor_type=k[3], variant=k[4],
                 stratum=k[5], tp=v[0], fp=v[1], fn=v[2], tn=v[3]) for k, v in sorted(counters.items())]
    write_csv(args.out / "confusion_counts.csv", rows)


def precision(args, run):
    with open(args.out / "confusion_counts.csv", encoding="utf-8") as f:
        counts = list(csv.DictReader(f))
    counters = {}
    for r in counts:
        key = tuple(r[k] for k in ("diagnostic_id", "ion_mode", "source", "precursor_type", "variant", "stratum"))
        cells = [float(r[k]) for k in ("tp", "fp", "fn", "tn")]
        counters[key] = cells
        pooled = key[:2] + ("pooled",) + key[3:]
        counters.setdefault(pooled, [0., 0., 0., 0.])
        counters[pooled] = [a + b for a, b in zip(counters[pooled], cells)]
    # Include empty cells, so an unobserved diagnostic is explicit, never silently absent.
    for d in diagnostics():
        for source in args.sources + ["pooled"]:
            for adduct in ["all"] + d["precursor_types"]:
                for variant in list(VARIANTS) + ["all_spectra_weighted", "formula_gated"]:
                    counters.setdefault((d["diagnostic_id"], d["ion_mode"], source, adduct, variant, "all"), [0., 0., 0., 0.])
    rows = []
    motif_ids = {d["diagnostic_id"]: d["motif_id"] for d in diagnostics()}
    for key, cells in sorted(counters.items()):
        row = dict(zip(("diagnostic_id", "ion_mode", "source", "precursor_type", "variant", "stratum"), key))
        row["motif_id"] = motif_ids[key[0]]
        row.update(statistics(*cells))
        sources = [statistics(*counters.get(key[:2] + (s,) + key[3:], [0., 0., 0., 0.])) for s in args.sources]
        row["status"] = status(row, sources if key[2] == "pooled" else ())
        rows.append(row)
    write_csv(args.out / "diagnostics_precision.csv", rows)
    # P1 measured multiple families; retain every family beside each diagnostic result.
    if P1_RESULT.exists():
        run["inputs"] = [r for r in run["inputs"] if r["name"] != "P1 motif_precision.csv"] + [dict(name="P1 motif_precision.csv", sha256=sha256(P1_RESULT))]
        with open(P1_RESULT, encoding="utf-8") as f:
            p1 = list(csv.DictReader(f))
        cross = []
        for r in rows:
            if (r["source"], r["precursor_type"], r["variant"], r["stratum"]) != ("pooled", "all", "default", "all"):
                continue
            for p in p1:
                if p["motif_id"] == r["motif_id"]:
                    cross.append(dict(diagnostic_id=r["diagnostic_id"], ion_mode=r["ion_mode"], motif_id=r["motif_id"],
                                      p1_family=p["family"], p1_n=p["n"], p1_precision=p["precision"],
                                      p1_wilson_lo=p["wilson_lo"], p1_status=p["status"],
                                      p2_n=r["n_obs"], p2_precision=r["precision"], p2_wilson_lo=r["wilson_lo"], p2_status=r["status"]))
        write_csv(args.out / "p1_crosscheck.csv", cross)


def candidate_bins(record):
    base = max((i for m, i in zip(record["mz"], record["intensity"]) if abs(m - record["precursor"]) > 1), default=0)
    if not base:
        return set()
    result = set()
    for mz, intensity in zip(record["mz"], record["intensity"]):
        if intensity / base < 0.01:
            continue
        result.add(("product_ion", int(math.floor(mz * 100 + 0.5))))
        loss = record["precursor"] - mz
        if loss >= 12:
            result.add(("neutral_loss", int(math.floor(loss * 100 + 0.5))))
    return result


def mining_compounds(records):
    """Collapse adduct representatives to one compound/mode, OR-ing 0.01-Da bins."""
    result = {}
    for r in records:
        if not r["representative"]:
            continue
        key = (r["key"], r["mode"])
        if key not in result:
            result[key] = dict(mode=r["mode"], motifs=set(r["motifs"]), bins=set())
        result[key]["bins"].update(candidate_bins(r))
    return list(result.values())


def select_candidates(vs20, motif_ids, min_positive=200, top=20):
    """No held-out argument: selection cannot depend on VS21."""
    selected = []
    for mode in ("pos", "neg"):
        records = [r for r in vs20 if r["mode"] == mode]
        total = Counter(b for r in records for b in r["bins"])
        for motif in motif_ids:
            positives = [r for r in records if motif in r["motifs"]]
            if len(positives) < min_positive:
                continue
            present = Counter(b for r in positives for b in r["bins"])
            candidates = [(b, n) for b, n in present.items() if n >= math.ceil(0.05 * len(positives))]
            candidates.sort(key=lambda bn: (-bn[1] / total[bn[0]], -bn[1], bn[0]))
            for (kind, mass_bin), tp in candidates[:top]:
                fp = total[(kind, mass_bin)] - tp
                s = statistics(tp, fp, len(positives) - tp, len(records) - len(positives) - fp)
                selected.append(dict(motif_id=motif, ion_mode=mode, kind=kind, mass_bin=mass_bin,
                                     mz_or_loss=mass_bin / 100, vs20_n=s["n_obs"], vs20_precision=s["precision"],
                                     vs20_positive_compounds=len(positives)))
    return selected


def score_candidates(selected, vs21):
    scored = []
    for c in selected:
        cells = [0, 0, 0, 0]
        for r in vs21:
            if r["mode"] != c["ion_mode"]:
                continue
            truth, obs = c["motif_id"] in r["motifs"], (c["kind"], c["mass_bin"]) in r["bins"]
            cells[0 if truth and obs else 1 if obs else 2 if truth else 3] += 1
        s = statistics(*cells)
        scored.append(dict(c, **{"vs21_" + k: v for k, v in s.items()}, status=status(s, n_min=30)))
    return scored


def assign_formulas(candidates):
    """Enumerate all nonnegative CHNOPS compositions up to the largest selected mass.

    Hydrogen is bounded only by remaining mass (no valence heuristic). Formula is
    unique only if exactly one composition is within 3 mDa of the bin centre.
    """
    from rdkit import Chem
    table = Chem.GetPeriodicTable()
    masses = {e: table.GetMostCommonIsotopeMass(e) for e in ("C", "H", "N", "O", "P", "S")}
    targets = sorted(set(c["mz_or_loss"] + (0 if c["kind"] == "neutral_loss" else 0.000549 if c["ion_mode"] == "pos" else -0.000549) for c in candidates))
    found = {t: [] for t in targets}
    if not targets:
        return
    # Stop after two formulas: ambiguity is proven and no further enumeration is needed.
    def visit(t, index, mass, counts):
        if len(found[t]) >= 2:
            return
        elements = ("S", "P", "N", "O", "C")
        if index == len(elements):
            h = round((t - mass) / masses["H"])
            if h >= 0 and abs(mass + h * masses["H"] - t) <= 0.003 + 1e-10:
                comp = dict(counts, H=h)
                formula = "".join(e + (str(comp[e]) if comp[e] != 1 else "") for e in ("C", "H", "N", "O", "P", "S") if comp.get(e, 0))
                if formula:
                    found[t].append(formula)
            return
        e = elements[index]
        for n in range(int((t + 0.003 - mass) / masses[e]) + 1):
            visit(t, index + 1, mass + n * masses[e], dict(counts, **{e: n}))
            if len(found[t]) >= 2:
                break
    for t in targets:
        visit(t, 0, 0, {})
    for c in candidates:
        t = c["mz_or_loss"] + (0 if c["kind"] == "neutral_loss" else 0.000549 if c["ion_mode"] == "pos" else -0.000549)
        c["formula"] = found[t][0] if len(found[t]) == 1 else ""
        c["formula_status"] = "unique" if len(found[t]) == 1 else "ambiguous" if found[t] else "unassigned"


def mine(args, run):
    public = mining_compounds(parquet_rows(row_path(args.out, "prepared", "vs20"))) if "vs20" in args.sources else []
    selected = select_candidates(public, [r["motif_id"] for r in motif_rows()])
    # Freeze discovery before reading the held-out spectra.
    write_json(args.out / "mined_selection-vs20.json", selected)
    private = mining_compounds(parquet_rows(row_path(args.out, "prepared", "vs21"))) if "vs21" in args.sources else []
    scored = score_candidates(selected, private)
    assign_formulas(scored)
    prevalence_rows = []
    for source, compounds in (("vs20", public), ("vs21", private)):
        if source not in args.sources:
            continue
        for mode in ("pos", "neg"):
            records = [r for r in compounds if r["mode"] == mode]
            for motif in motif_rows():
                n = sum(motif["motif_id"] in r["motifs"] for r in records)
                prevalence_rows.append(dict(source=source, ion_mode=mode, motif_id=motif["motif_id"],
                                            n_compounds=len(records), n_positive=n,
                                            prevalence=n / len(records) if records else None))
    write_csv(args.out / "motif_prevalence.csv", prevalence_rows)
    fields = list(scored[0]) if scored else ["motif_id", "ion_mode", "kind", "mz_or_loss", "formula", "vs20_n", "vs20_precision", "vs21_n_obs", "vs21_precision", "vs21_wilson_lo", "status"]
    write_csv(args.out / "mined_candidates.csv", scored, fields)
    run["mining"] = dict(selected=len(selected), status_counts=dict(Counter(r["status"] for r in scored)),
                         selection_source="vs20", evaluation_source="vs21", n_min=30,
                         positive_compound_min=200, positive_frequency_min=0.05, top_per_motif_mode=20,
                         observation="same 0.01-Da bin, >=1% relative intensity; OR across all precursor-type representatives")
    print("mine", run["mining"]["status_counts"], flush=True)


def choose(args, run):
    with open(args.out / "diagnostics_precision.csv", encoding="utf-8") as f:
        raw = list(csv.DictReader(f))
    numeric = set(statistics(0, 0, 0, 0))
    rows = [{k: (float(v) if v else None) if k in numeric else v for k, v in r.items()} for r in raw]
    by_pair = []
    for r in rows:
        if (r["source"], r["precursor_type"], r["variant"], r["stratum"]) != ("pooled", "all", "default", "all"):
            continue
        sources = [s for s in rows if s["diagnostic_id"] == r["diagnostic_id"] and s["source"] in args.sources
                   and (s["precursor_type"], s["variant"], s["stratum"]) == ("all", "default", "all")]
        s = statistics(*(r[k] for k in ("tp", "fp", "fn", "tn")))
        by_pair.append(dict(diagnostic_id=r["diagnostic_id"], ion_mode=r["ion_mode"], motif_id=r["motif_id"],
                            **s, status=status(s, sources), by_source={x["source"]: {k: x[k] for k in numeric} for x in sources}))
    result = dict(label="ws8-p2-diagnostics", run_id=args.out.name, rule=RULE,
                  table_hashes={p.name: sha256(p) for p in (HERE / "diagnostics-p2.tsv", P1_MOTIFS, HERE / "motifs-p2-extra.tsv")},
                  by_pair=by_pair, si={"diagnostics_allowed": [dict(diagnostic_id=r["diagnostic_id"], ion_mode=r["ion_mode"])
                                                              for r in by_pair if r["status"] == "pass"]})
    write_json(args.out / "diagnostics_precision.json", result)
    audit_outputs(args.out)
    report(args, run, result, rows)


def audit_outputs(out):
    """P2-A4/A5: recompute choices from sufficient statistics and scan aggregates."""
    result = json.loads((out / "diagnostics_precision.json").read_text(encoding="utf-8"))
    with open(out / "diagnostics_precision.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for r in result["by_pair"]:
        raw = next(x for x in rows if x["diagnostic_id"] == r["diagnostic_id"] and
                   (x["source"], x["precursor_type"], x["variant"], x["stratum"]) == ("pooled", "all", "default", "all"))
        calc = statistics(*(float(raw[k]) for k in ("tp", "fp", "fn", "tn")))
        for k, v in calc.items():
            assert (r[k] is None and v is None) or (r[k] is not None and v is not None and math.isclose(r[k], v, abs_tol=1e-12)), k
        assert r["status"] == raw["status"] == status(calc, list(r["by_source"].values()))
    for p in out.iterdir():
        if p.suffix in (".csv", ".json", ".md"):
            content = p.read_text(encoding="utf-8")
            assert not re.search(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", content), "Aggregate contains an identifier"
            if p.suffix == ".csv":
                columns = next(csv.reader(content.splitlines()), [])
                assert not any(re.search(r"(^|_)(key|first.?block|smiles|name)($|_)", c, re.I) for c in columns), "Forbidden aggregate column"
    return dict(P2_A4="pass", P2_A5="pass")


def report(args, run, result, rows):
    lines = ["# WS8 P2 diagnostic evidence", "", f"Run: `{args.out.name}`. Library-only measurements; mined candidates are report only.", "",
             "| Diagnostic | Mode | n_obs | Precision | Wilson lower | Recall | LR+ | Prior stress | Status |",
             "|---|---|---:|---:|---:|---:|---:|---:|---|"]
    def fmt(v):
        return "NA" if v is None else f"{v:.4g}"
    for r in result["by_pair"]:
        lines.append("| " + " | ".join([r["diagnostic_id"], r["ion_mode"]] + [fmt(r[k]) for k in ("n_obs", "precision", "wilson_lo", "recall", "lr_plus", "prior_stress_precision")] + [r["status"]]) + " |")
    lines += ["", "Counts use one representative per constitution, mode and precursor type. Mode-level pools sum those units across allowed precursor types.",
              "Weighted all-spectra Wilson/log intervals are descriptive sensitivity intervals using fractional counts, not independent replicate confidence intervals.",
              "Formula gating requires at least the residue's heavy-element counts; H is ignored because charge/adduct changes H stoichiometry.",
              "Collision energies are stratified by explicit units; stepped/ranged/mixed fields remain unparseable. Ties in representative peak counts keep the earliest record.",
              "Mined bins are 0.01 Da with half-up rounding and a >=1% floor. A compound is counted once per mode, OR across all precursor types; NL is PRECURSORMZ minus peak m/z and >=12 Da.",
              "Discovery keeps the top 20 combined product-ion/loss bins per motif/mode by VS20 precision, then positive support, kind and bin. Formula uniqueness searches all nonnegative CHNOPS compositions within 3 mDa of the bin centre.",
              "", "Mining: " + json.dumps(run.get("mining", {}), sort_keys=True), "",
              "Exclusions (record units except explicitly unique-SMILES counters):", "", "```json", json.dumps(run.get("exclusions", {}), indent=2), "```", "",
              "P1 cross-check is in p1_crosscheck.csv; all specified tolerance, intensity, formula, instrument and energy variants are in diagnostics_precision.csv.",
              "Both the broad P1 pentoside and strict P2 pentopyranoside are reported in motif_prevalence.csv (compounds per mode).",
              "The formula-derived negative taurine ion is 124.007384 Da, correcting the neutral-mass rounding shown in the starting examples; the formula/electron rule is authoritative."]
    (args.out / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pairs = result["by_pair"]
    fig, axes = plt.subplots(1, 2, figsize=(13, max(6, len(pairs) * 0.25)), sharey=True)
    for i, r in enumerate(pairs):
        if r["precision"] is not None:
            axes[0].errorbar(r["precision"], i, xerr=[[max(0, r["precision"] - r["wilson_lo"] )], [max(0, r["wilson_hi"] - r["precision"])]], fmt="o", color="tab:blue")
        if r["recall"] is not None:
            lo, hi = wilson(r["tp"], r["tp"] + r["fn"])
            axes[1].errorbar(r["recall"], i, xerr=[[max(0, r["recall"] - lo)], [max(0, hi - r["recall"])]], fmt="o", color="tab:orange")
    axes[0].set_yticks(range(len(pairs)), [r["diagnostic_id"] for r in pairs])
    for ax, title in zip(axes, ("Precision (Wilson 95%)", "Recall (Wilson 95%)")):
        ax.set(xlim=(-0.02, 1.02), title=title)
        ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    (args.out / "figures").mkdir(exist_ok=True)
    fig.savefig(args.out / "figures/precision_recall.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "observe", "precision", "mine", "choose", "all"))
    parser.add_argument("--sources", default="vs20,vs21")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=float, help="Deterministic per-file fraction (0,1]; 0.01 for smoke")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--detach", action="store_true", help="Launch outside the caller's Windows job and log to run.log")
    parser.add_argument("--detached-child", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.sources = args.sources.split(",")
    if args.sources not in (["vs20"], ["vs21"], ["vs20", "vs21"]) or not 1 <= args.workers <= 8:
        parser.error("Use --sources vs20[,vs21] or vs21 and 1..8 workers")
    if args.limit is not None and not 0 < args.limit <= 1:
        parser.error("--limit must be a fraction in (0,1]")
    args.out = args.out.resolve()
    if args.out.parent != OUT_ROOT.resolve() or not re.fullmatch(r"p2(?:smoke)?-\d{8}-[a-f0-9]{7,40}", args.out.name):
        parser.error("--out must be the assigned ws8-p2 run directory")
    args.out.mkdir(parents=True, exist_ok=True)
    if args.detach:
        if sys.platform != "win32":
            parser.error("Detached launch requires Windows")
        command = [sys.executable, "-B", str(Path(__file__).resolve())] + [x for x in sys.argv[1:] if x != "--detach"] + ["--detached-child"]
        # WMI's service creates the child independently of nested application jobs.
        # CREATE_BREAKAWAY_FROM_JOB alone did not escape all host jobs in the smoke.
        env = dict(os.environ, WS8_P2_COMMAND=subprocess.list2cmdline(command), WS8_P2_CWD=str(ROOT))
        script = "$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$env:WS8_P2_COMMAND; CurrentDirectory=$env:WS8_P2_CWD}; if ($r.ReturnValue -ne 0) { throw 'WMI detached launch failed' }; Write-Output $r.ProcessId"
        pid = subprocess.check_output(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", script],
                                      env=env, text=True, creationflags=0x08000000).strip()
        print("Detached P2 process:", pid, flush=True)
        return
    if args.detached_child:
        sys.stdout = open(args.out / "run.log", "a", encoding="utf-8", buffering=1)
        sys.stderr = sys.stdout
    run_path = args.out / "run.json"
    run = json.loads(run_path.read_text(encoding="utf-8")) if run_path.exists() else {}
    if run.get("status") == "complete":
        parser.error("Completed run is immutable; use a new code commit/run ID")
    if "parameters" in run:
        for field, value in (("sources", args.sources), ("workers", args.workers), ("limit", args.limit)):
            if run["parameters"].get(field) != value:
                parser.error("Stage arguments must match the preparation manifest")
    priority_name = priority()
    context = process_context()
    if args.detached_child and context["launcher"]["in_windows_job"]:
        raise RuntimeError("Detached venv launcher still belongs to a Windows job")
    import rdkit, numpy, matchms, pandas, pyarrow
    run.update(label="ws8-p2", run_id=args.out.name, code_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
               interpreter=sys.executable, versions=dict(python=sys.version, rdkit=rdkit.__version__, numpy=numpy.__version__,
                                                        matchms=matchms.__version__, pandas=pandas.__version__, pyarrow=pyarrow.__version__),
               parameters=dict(sources=args.sources, workers=args.workers, limit=args.limit, priority=priority_name,
                               threads={v: os.environ[v] for v in THREAD_VARS}, variants=VARIANTS, rule=RULE,
                               resolution="all peaks within 0.001 Da of integer OR all lexical m/z <2 decimals",
                               ion_mode_aliases="P/Positive/pos => pos; N/Negative/neg => neg; blank uses file mode",
                               representative="most positive finite peaks after resolution gate, earliest record tie",
                               pooling="one unit per first block/mode/precursor type, summed across types",
                               formula_gate="heavy-element residue counts; H ignored", electron_mass=0.000549,
                               precursor_base_exclusion_da=1.0, neutral_loss_min_below_precursor_da=1.0,
                               standardisation=dict(implementation="tools.ws8_p1.mcs_core.standardise",
                                                    min_heavy_atoms=4, max_heavy_atoms=80, remove_stereo=True,
                                                    remove_isotopes=True, reject_radicals=True,
                                                    largest_fragment=True, uncharger="RDKit defaults",
                                                    allowed_elements="C H N O P S F Cl Br I Si B Se"),
                               multiprocessing="spawn", worker_priority="Win32 Set/GetPriorityClass verified by initializer",
                               mining=dict(bin_da=0.01, bin_rounding="half up", min_rel_intensity=0.01,
                                           min_neutral_loss_da=12, min_positive_compounds=200,
                                           min_positive_frequency=0.05, top_per_motif_mode=20,
                                           ranking="precision, positive support, kind, bin",
                                           formula_tolerance_da=0.003, formula_elements="C H N O P S",
                                           formula_compositions="nonnegative; no valence heuristic",
                                           heldout_n_min=30, discovery="vs20 only", evaluation="vs21 only")),
               process=dict(detached_child=args.detached_child, launch_method="WMI" if args.detached_child else "foreground", **context),
               status="running")
    run.setdefault("timings_seconds", {})
    write_json(run_path, run)
    stages = ("prepare", "observe", "precision", "mine", "choose") if args.stage == "all" else (args.stage,)
    try:
        for stage in stages:
            start = time.perf_counter()
            globals()[stage](args, run)
            run["timings_seconds"][stage] = time.perf_counter() - start
            write_json(run_path, run)
        if args.stage in ("choose", "all"):
            run["acceptance_artifact_audit"] = audit_outputs(args.out)
            run["status"] = "complete"
        run["files"] = {str(p.relative_to(args.out)).replace("\\", "/"): sha256(p) for p in sorted(args.out.rglob("*"))
                        if p.is_file() and p.name not in ("run.json", "run.log")}
        write_json(run_path, run)
        print("completed", args.stage, run["timings_seconds"], flush=True)
    except Exception as exc:
        run.update(status="failed", error_type=type(exc).__name__)
        write_json(run_path, run)
        # Never expose an exception containing real library data or paths in a shared log.
        print("FAILED", type(exc).__name__, flush=True)
        import traceback
        for frame in traceback.extract_tb(exc.__traceback__):
            print(Path(frame.filename).name, frame.lineno, frame.name, flush=True)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
