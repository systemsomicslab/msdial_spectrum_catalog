"""WS8 P1: local-only, reproducible SI size calibration. No package dependencies on import."""
import argparse
import collections
import csv
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import re
import subprocess
import sys
import time

try:
    from .mcs_core import mcs_params, mcs_params_sha256
except ImportError:
    from mcs_core import mcs_params, mcs_params_sha256

HERE = Path(__file__).resolve().parent
OUTPUT_ROOT = Path("F:/MSDIAL_Public_Reanalysis/annotation/ws8-p1")
LIB_ROOT = Path("D:/13_MSDIAL_Public_Reanalysis/libraries")
RESOURCE = Path("D:/0_SourceCode/msemblator/external/msfinder/MSFINDER ver 3.61/Resources")
DB_SHA = "2255e022b5b34937d361b23f722b09f7a5bc600654521c2e253d780484c68e7e"
FAMILIES = ("spectral_nn", "spectral_nn_any", "spectral_nn_scaffold", "same_formula_spectral",
            "same_formula_db_similar", "same_formula_db_random", "truth_in_set")
BAR_FAMILIES = ("spectral_nn", "same_formula_db_random")
BINS = ("4-15", "16-25", "26-40", "41-80")
THREAD_VARS = ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
for _var in THREAD_VARS:
    os.environ[_var] = "1"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def wilson(success, n):
    if not n:
        return (0.0, 1.0)
    z = 1.959963984540054
    p, z2 = success / n, z * z
    center = (p + z2 / (2 * n)) / (1 + z2 / n)
    half = z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    return max(0.0, center - half), min(1.0, center + half)


def quantile(values, q):
    """Linear sample quantile, identical to numpy's default, without numpy."""
    if not values:
        return None
    values = sorted(values)
    index = (len(values) - 1) * q
    lo = int(index)
    return values[lo] + (values[min(lo + 1, len(values) - 1)] - values[lo]) * (index - lo)


def cumulative(values, threshold):
    selected = [ok for size, ok in values if size >= threshold]
    success = sum(selected)
    return {"n": len(selected), "success": success, "precision": success / len(selected) if selected else None,
            "wilson_lo": wilson(success, len(selected))[0], "wilson_hi": wilson(success, len(selected))[1]}


def choose_axis(curves, chance95, target=0.95, n_min=100):
    """Both required families must support the returned threshold on written curves."""
    if chance95 is None:
        return None
    by_family = {f: {float(r["threshold"]): r for r in curves if r["family"] == f} for f in BAR_FAMILIES}
    first = {}
    for family, rows in by_family.items():
        passing = [x for x, r in rows.items() if x > chance95 and int(r["n"]) >= n_min
                   and float(r["wilson_lo"]) >= target]
        if not passing:
            return None
        first[family] = min(passing)
    threshold = max(first.values())
    common = sorted(set.intersection(*(set(rows) for rows in by_family.values())))
    # Wilson confidence need not be monotone: recheck both families at the final bar.
    for x in common:
        if x < threshold:
            continue
        if all(int(rows[x]["n"]) >= n_min and float(rows[x]["wilson_lo"]) >= target
               for rows in by_family.values()):
            deciding = max(BAR_FAMILIES, key=lambda f: first[f])
            row = by_family[deciding][x]
            return {"value": x, "deciding_family": deciding, "n": int(row["n"]),
                    "wilson_lo": float(row["wilson_lo"]), "family_bars": first,
                    "evidence": {f: {"n": int(rows[x]["n"]), "wilson_lo": float(rows[x]["wilson_lo"])}
                                 for f, rows in by_family.items()}}
    return None


def priority():
    if sys.platform != "win32":
        raise RuntimeError("P1 real runs require Windows BelowNormal priority")
    import ctypes
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    kernel.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
    kernel.GetPriorityClass.argtypes = (ctypes.c_void_p,)
    handle = kernel.GetCurrentProcess()
    if not kernel.SetPriorityClass(handle, 0x4000) or kernel.GetPriorityClass(handle) != 0x4000:
        raise RuntimeError("Cannot verify BelowNormal priority")
    return "BelowNormal"


def worker_init():
    priority()
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")  # Never print real-data identifiers or SMILES.


def chem_record(smiles):
    try:
        try:
            from .mcs_core import standardise, first_block, hill_formula, RejectedStructure
        except ImportError:
            from mcs_core import standardise, first_block, hill_formula, RejectedStructure
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = standardise(smiles)
        return {"key": first_block(mol), "smiles": Chem.MolToSmiles(mol), "formula": hill_formula(mol),
                "atoms": mol.GetNumHeavyAtoms(), "mass": int(round(Descriptors.ExactMolWt(mol)))}
    except Exception as exc:
        return {"rejected": str(exc) if isinstance(exc, RejectedStructure) else "standardisation"}


def atom_bin(n):
    return BINS[0 if n <= 15 else 1 if n <= 25 else 2 if n <= 40 else 3]


def sample_count(available, cap, pilot):
    return math.ceil(min(available, cap) * (0.05 if pilot else 1))


def seed_generators(seed, sources):
    import numpy as np
    # Fixed registry, independent of missing cells and scheduling.
    labels = [(s, stratum, f, k) for s in ("vs20", "vs21", "db")
              for stratum in [f"{lipid}:{b}" for lipid in (0, 1) for b in BINS] + ["all"]
              for f in ("truths", "motifs", "unconditional") + FAMILIES for k in (0, 2, 3, 4, 5)]
    children = np.random.SeedSequence(seed).spawn(len(labels))
    return {label: np.random.Generator(np.random.PCG64(child)) for label, child in zip(labels, children)}


def dump(path, data):
    with open(path, "wb") as stream:
        pickle.dump(data, stream, protocol=5)


def load(path):
    with open(path, "rb") as stream:
        return pickle.load(stream)


def iter_msp(path):
    """Streaming allowlist parser; record titles are used only for a boolean lipid flag."""
    record, peaks, in_peaks = {}, [], False
    allowed = {"name", "precursormz", "precursortype", "ionmode", "smiles"}
    with open(path, encoding="utf-8-sig", errors="replace") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                if record:
                    yield record, peaks
                record, peaks, in_peaks = {}, [], False
                continue
            if in_peaks:
                try:
                    parts = line.split()
                    mz, intensity = float(parts[0]), float(parts[1])
                    if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
                        peaks.append((mz, intensity))
                except (ValueError, IndexError):
                    pass
            elif ":" in line:
                key, value = line.split(":", 1)
                key, value = key.lower().strip(), value.strip()
                if key in allowed:
                    record[key] = value
                if key in ("num peaks", "numpeaks"):
                    in_peaks = True
        if record:
            yield record, peaks


def motifs():
    from rdkit import Chem
    with open(HERE / "motifs-p1.tsv", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    result = {r["motif_id"]: Chem.MolFromSmarts(r["smarts"]) for r in rows}
    if any(q is None for q in result.values()):
        raise ValueError("Invalid motif SMARTS")
    return result


def prepare(args, run):
    import pandas as pd
    db_path = RESOURCE / "MsfinderStructureDB-VS15.esd"
    if sha256(db_path) != DB_SHA:
        raise RuntimeError("Structure DB checksum mismatch; stopped")
    inputs = [{"name": db_path.name, "sha256": DB_SHA},
              {"name": "InchikeyClassyfireDB-VS5.icd", "sha256": sha256(RESOURCE / "InchikeyClassyfireDB-VS5.icd")},
              {"name": "motifs-p1.tsv", "sha256": sha256(HERE / "motifs-p1.tsv")}]
    rejection = collections.Counter()
    raw_db = []
    with open(db_path, encoding="utf-8-sig") as stream:
        next(stream)
        for line in stream:
            row = line.rstrip("\n").split("\t")
            raw_db.append((row[6], row[18] not in ("N/A", "", "NA")))
    db, lipid_keys = {}, set()
    with mp.get_context("spawn").Pool(args.workers, initializer=worker_init) as pool:
        for raw, record in zip(raw_db, pool.imap(chem_record, (r[0] for r in raw_db), chunksize=128)):
            if "rejected" in record:
                rejection["db:" + record["rejected"]] += 1
                continue
            db.setdefault(record["key"], record)
            if raw[1]:
                lipid_keys.add(record["key"])
    del raw_db
    # Reporting join only; no parent label is emitted to aggregate files.
    parents = {}
    with open(RESOURCE / "InchikeyClassyfireDB-VS5.icd", encoding="utf-8-sig") as stream:
        next(stream)
        for line in stream:
            row = line.rstrip("\n").split("\t")
            if len(row) >= 2:
                parents[row[0].split("-")[0]] = row[1]
    libraries, all_truths, public_keys = {}, {}, set()
    for source in ("vs20", "vs21"):
        # VS20 is always parsed for the VS21 exclusion, even when only VS21 is requested.
        if source == "vs21" and source not in args.sources:
            continue
        paths = ([LIB_ROOT / "21901200/MSMS-Public_all-pos-VS20.msp",
                  LIB_ROOT / "21904103/MSMS-Public_all-neg-VS20.msp"] if source == "vs20" else
                 [next(LIB_ROOT.rglob(f"MSMS-Private-{mode}-VS21.msp")) for mode in ("pos", "neg")])
        # First reduce repeated spectra by raw SMILES and ion mode, then standardise once.
        raw_reps = {}
        for path in paths:
            inputs.append({"name": path.name, "sha256": sha256(path)})
            for metadata, peaks in iter_msp(path):
                adduct = metadata.get("precursortype", "").replace(" ", "")
                if adduct not in ("[M+H]+", "[M-H]-") or not peaks:
                    rejection[source + ":unusable_spectrum"] += 1
                    continue
                try:
                    precursor = float(metadata.get("precursormz", ""))
                    if not math.isfinite(precursor) or precursor <= 0:
                        raise ValueError
                except ValueError:
                    rejection[source + ":precursor"] += 1
                    continue
                identity = (metadata.get("smiles", ""), adduct)
                lipid = bool(re.match(r"^[A-Za-z][A-Za-z0-9-]*\s+\d+:\d+", metadata.get("name", "")))
                lipid = lipid or (identity in raw_reps and raw_reps[identity]["lipid_name"])
                if identity not in raw_reps or len(peaks) > len(raw_reps[identity]["peaks"]):
                    raw_reps[identity] = {"adduct": adduct, "mode": "pos" if adduct == "[M+H]+" else "neg",
                                          "precursor": precursor, "peaks": sorted(peaks), "lipid_name": lipid}
                elif lipid:
                    raw_reps[identity]["lipid_name"] = True
        reps = {}
        with mp.get_context("spawn").Pool(args.workers, initializer=worker_init) as pool:
            for (identity, spectrum), record in zip(raw_reps.items(), pool.imap(
                    chem_record, (x[0] for x in raw_reps), chunksize=64)):
                if "rejected" in record:
                    rejection[source + ":" + record["rejected"]] += 1
                    continue
                record.update(spectrum)
                key = (record["key"], record["adduct"])
                if key in reps:
                    record["lipid_name"] = record["lipid_name"] or reps[key]["lipid_name"]
                if key not in reps or len(record["peaks"]) > len(reps[key]["peaks"]):
                    reps[key] = record
                elif record["lipid_name"]:
                    reps[key]["lipid_name"] = True
        del raw_reps
        truths = {}
        for record in reps.values():
            key = record["key"]
            if source == "vs21" and key in public_keys:
                continue
            if key not in truths or len(record["peaks"]) > len(truths[key]["peaks"]):
                truths[key] = dict(record, lipid_name=record["lipid_name"] or (key in truths and truths[key]["lipid_name"]))
            elif record["lipid_name"]:
                truths[key]["lipid_name"] = True
        if source == "vs20":
            public_keys = set(truths)
        for key, record in truths.items():
            record["lipid"] = int(key in lipid_keys or record["lipid_name"])
            record["lipid_source"] = "lipidmaps" if key in lipid_keys else "shorthand" if record["lipid_name"] else "none"
            record["heavy_bin"] = atom_bin(record["atoms"])
            record["stratum"] = f'{record["lipid"]}:{record["heavy_bin"]}'
            record["direct_parent"] = parents.get(key, "")
            record["source"] = source
        if source in args.sources:
            libraries[source] = list(reps.values())
            all_truths[source] = truths
    generators = seed_generators(args.seed, args.sources)
    from rdkit import Chem
    motif_queries = motifs()
    selected = []
    for source, truths in all_truths.items():
        chosen = set()
        for stratum in sorted({r["stratum"] for r in truths.values()}):
            candidates = sorted(k for k, r in truths.items() if r["stratum"] == stratum)
            count = sample_count(len(candidates), args.cap, args.pilot)
            rng = generators[source, stratum, "truths", 0]
            for index in rng.choice(len(candidates), count, replace=False):
                key = candidates[index]
                chosen.add(key)
                selected.append(dict(truths[key], motif_oversample=False))
        # Extra truths affect motif precision alone; each truth is sampled at most once.
        for motif_id, query in motif_queries.items():
            candidates = [k for k in sorted(truths) if k not in chosen
                          and Chem.MolFromSmiles(truths[k]["smiles"]).HasSubstructMatch(query)]
            rng = generators[source, "all", "motifs", 0]
            count = sample_count(len(candidates), 150, args.pilot)
            for index in rng.choice(len(candidates), count, replace=False):
                key = candidates[index]
                chosen.add(key)
                selected.append(dict(truths[key], motif_oversample=True))
        rows = [{k: v for k, v in r.items() if k not in ("peaks",)} for r in selected if r["source"] == source]
        pd.DataFrame(rows).to_parquet(args.out / row_filename("truths", source), index=False)
    dump(args.out / "prepared.local-vs21.pkl", {"db": db, "libraries": libraries, "truths": selected})
    dump(args.out / "libraries.local-vs21.pkl", libraries)
    run.update(inputs=inputs, rejection_counts=dict(rejection), truth_counts=dict(collections.Counter(
        r["source"] + ":" + r["stratum"] + (":motif" if r["motif_oversample"] else "") for r in selected)),
               representative_counts={s: len(rows) for s, rows in libraries.items()}, db_usable=len(db))


def row_filename(stem, source):
    return stem + (".local-vs21" if source == "vs21" else "") + ".parquet"


def neighbour_init(path):
    global SPECTRA, KERNEL
    worker_init()
    import logging
    import numpy as np
    from matchms import Spectrum
    from matchms.similarity import ModifiedCosineGreedy
    logging.getLogger("matchms").setLevel(logging.ERROR)
    SPECTRA = {}
    for source, rows in load(path).items():
        groups = collections.defaultdict(list)
        for r in rows:
            peaks = np.array(r["peaks"], dtype=float)
            spectrum = Spectrum(mz=peaks[:, 0], intensities=peaks[:, 1],
                                metadata={"precursor_mz": r["precursor"]}, metadata_harmonization=False)
            groups[r["adduct"]].append((r, spectrum))
        SPECTRA[source] = groups
    KERNEL = ModifiedCosineGreedy(tolerance=0.02, mz_power=0.0, intensity_power=1.0)


def neighbour_query(truth):
    import numpy as np
    from matchms import Spectrum
    peaks = np.array(truth["peaks"], dtype=float)
    query = Spectrum(mz=peaks[:, 0], intensities=peaks[:, 1],
                     metadata={"precursor_mz": truth["precursor"]}, metadata_harmonization=False)
    ranked, pairs = [], 0
    for candidate, spectrum in SPECTRA[truth["source"]][truth["adduct"]]:
        if candidate["key"] == truth["key"] or abs(candidate["precursor"] - truth["precursor"]) > 200:
            continue
        pairs += 1
        score = KERNEL.pair(spectrum, query)
        if int(score["matches"]) >= 6:
            ranked.append((float(score["score"]), candidate["key"], candidate["formula"], candidate["smiles"]))
    ranked.sort(key=lambda x: (-x[0], x[1]))
    return {"source": truth["source"], "key": truth["key"], "pairs": pairs,
            "spectral_nn": [r for r in ranked if r[0] >= 0.6][:5], "spectral_nn_any": ranked[:5],
            "same_formula_spectral": [r for r in ranked if r[2] == truth["formula"]][:5]}


def neighbours(args, run):
    truths = load(args.out / "prepared.local-vs21.pkl")["truths"]
    results = []
    with mp.get_context("spawn").Pool(args.workers, initializer=neighbour_init,
                                     initargs=(args.out / "libraries.local-vs21.pkl",)) as pool:
        for index, row in enumerate(pool.imap(neighbour_query, truths, chunksize=1), 1):
            results.append(row)
            if index % 100 == 0:
                print(f"neighbours: {index}/{len(truths)} queries", flush=True)
    dump(args.out / "neighbours.local-vs21.pkl", results)
    run["spectral_pairs"] = sum(r["pairs"] for r in results)


def sets(args, run):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    try:
        from .mcs_core import first_block
    except ImportError:
        from mcs_core import first_block
    data = load(args.out / "prepared.local-vs21.pkl")
    db = data["db"]
    by_formula, by_mass = collections.defaultdict(list), collections.defaultdict(list)
    for key in sorted(db):
        r = db[key]
        by_formula[r["formula"]].append(key)
        by_mass[r["mass"]].append(key)
    nn = {(r["source"], r["key"]): r for r in load(args.out / "neighbours.local-vs21.pkl")}
    rngs = seed_generators(args.seed, args.sources)
    fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = {}
    pending, shortfalls = [], collections.Counter()
    def fingerprint(smiles):
        if smiles not in fps:
            fps[smiles] = fp_generator.GetFingerprint(Chem.MolFromSmiles(smiles))
        return fps[smiles]
    for truth in data["truths"]:
        source, stratum, key = truth["source"], truth["stratum"], truth["key"]
        neighbours_row = nn[source, key]
        members = {f: [(r[1], r[3]) for r in neighbours_row[f]]
                   for f in ("spectral_nn", "spectral_nn_any", "same_formula_spectral")}
        members["spectral_nn_scaffold"] = []
        # Do not replace missing scaffolds with later neighbours: preserve the first-k definition.
        for member_key, smiles in members["spectral_nn"]:
            scaffold = MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smiles))
            scaffold_key = first_block(scaffold) if scaffold.GetNumAtoms() else ""
            members["spectral_nn_scaffold"].append((scaffold_key, Chem.MolToSmiles(scaffold) if scaffold_key != key else ""))
        exact = [x for x in by_formula[truth["formula"]] if x != key]
        truth_fp = fingerprint(truth["smiles"])
        similar = sorted(exact, key=lambda x: (-DataStructs.TanimotoSimilarity(truth_fp, fingerprint(db[x]["smiles"])), x))
        members["same_formula_db_similar"] = [(x, db[x]["smiles"]) for x in similar[:5]]
        random_order, relaxation = [], {}
        # Incremental draws use a distinct SeedSequence child for each k while preserving nesting.
        for k in (2, 3, 4, 5):
            pool_keys, level = exact, "none"
            if len(exact) < k:
                pool_keys = sorted(set(exact) | set(by_mass[truth["mass"]]))
                level = "nominal_mass_0"
                if len([x for x in pool_keys if x != key]) < k:
                    pool_keys = sorted(set(pool_keys) | set(by_mass[truth["mass"] - 1]) | set(by_mass[truth["mass"] + 1]))
                    level = "nominal_mass_1"
            available = [x for x in pool_keys if x != key and x not in random_order]
            rng = rngs[source, stratum, "same_formula_db_random", k]
            need = min(k - len(random_order), len(available))
            random_order.extend(available[int(i)] for i in rng.choice(len(available), need, replace=False))
            relaxation[k] = level
        members["same_formula_db_random"] = [(x, db[x]["smiles"]) for x in random_order]
        members["truth_in_set"] = [(key, truth["smiles"])] + members["same_formula_db_random"][:4]
        for family in FAMILIES:
            for k in (2, 3, 4, 5):
                selected = members[family][:k]
                if len(selected) < k or any(not smiles for _, smiles in selected):
                    shortfalls[f"{source}:{stratum}:{family}:{k}"] += 1
                    continue
                if family != "truth_in_set" and key in [x[0] for x in selected]:
                    raise AssertionError("Truth leaked into a candidate set")
                pending.append({"source": source, "family": family, "k": k, "truth_key": key,
                                "truth_smiles": truth["smiles"], "truth_atoms": truth["atoms"],
                                "member_keys": [x[0] for x in selected], "member_smiles": [x[1] for x in selected],
                                "lipid": truth["lipid"], "heavy_bin": truth["heavy_bin"], "stratum": stratum,
                                "formula_relaxed": relaxation[k] if family == "same_formula_db_random" else "none",
                                "motif_oversample": truth["motif_oversample"], "chance_kind": "matched" if family == "same_formula_db_random" else ""})
    for k in (2, 3, 4, 5):
        for b in BINS:
            formulas = sorted(f for f, keys in by_formula.items() if len(keys) >= k + 1 and atom_bin(db[keys[0]]["atoms"]) == b)
            rng = rngs["db", "0:" + b, "unconditional", k]
            count = math.ceil(1000 * (0.05 if args.pilot else 1)) if formulas else 0
            for _ in range(count):
                formula = formulas[int(rng.integers(len(formulas)))]
                keys = by_formula[formula]
                selected = [keys[int(i)] for i in rng.choice(len(keys), k, replace=False)]
                pending.append({"source": "db", "family": "same_formula_db_random", "k": k,
                                "truth_key": "", "truth_smiles": "", "truth_atoms": 0,
                                "member_keys": selected, "member_smiles": [db[x]["smiles"] for x in selected],
                                "lipid": -1, "heavy_bin": b, "stratum": "all:" + b, "formula_relaxed": "none",
                                "motif_oversample": False, "chance_kind": "unconditional"})
            if not formulas:
                shortfalls[f"unconditional:{b}:{k}"] += 1000
    dump(args.out / "pending.local-vs21.pkl", pending)
    run.update(shortfalls=dict(shortfalls), requested_mcs=len(pending), seed_registry="SeedSequence.spawn: source,stratum,family,k; fixed lexical registry in code")


def mcs_init():
    global MOTIF_QUERIES
    worker_init()
    MOTIF_QUERIES = motifs()


def mcs_one(row):
    from rdkit import Chem
    try:
        from .mcs_core import mcs_over_set, contained
    except ImportError:
        from mcs_core import mcs_over_set, contained
    members = [Chem.MolFromSmiles(s) for s in row["member_smiles"]]
    truth = Chem.MolFromSmiles(row["truth_smiles"]) if row["truth_smiles"] else None
    start = time.monotonic()
    smarts, atoms, timed_out = mcs_over_set(members)
    elapsed = time.monotonic() - start
    smallest = min(m.GetNumHeavyAtoms() for m in members)
    result = dict(row, mcs_smarts=smarts, mcs_atoms=atoms, smallest_atoms=smallest,
                  fraction_smallest=atoms / smallest if smallest else 0,
                  fraction_truth=atoms / row["truth_atoms"] if truth is not None else None,
                  timed_out=timed_out, contained=contained(truth, smarts) if truth is not None and not timed_out else None,
                  motif_matches=[motif for motif, q in MOTIF_QUERIES.items() if all(m.HasSubstructMatch(q) for m in members)],
                  truth_motifs=[motif for motif, q in MOTIF_QUERIES.items() if truth is not None and truth.HasSubstructMatch(q)])
    # Runtime is kept outside the deterministic sets artifact for P1-A4.
    return result, elapsed


def mcs(args, run):
    import pandas as pd
    pending = load(args.out / "pending.local-vs21.pkl")
    results, timings = [], []
    with mp.get_context("spawn").Pool(args.workers, initializer=mcs_init) as pool:
        for index, (result, elapsed) in enumerate(pool.imap(mcs_one, pending, chunksize=1), 1):
            results.append(result)
            timings.append(elapsed)
            if index % 1000 == 0:
                print(f"mcs: {index}/{len(pending)} sets", flush=True)
    for source in args.sources + ["db"]:
        rows = [r for r in results if r["source"] == source]
        pd.DataFrame(rows).to_parquet(args.out / row_filename("sets" if source != "db" else "chance_sets", source), index=False)
    dump(args.out / "results.local-vs21.pkl", results)
    run.update(mcs_median_seconds=quantile(timings, 0.5), mcs_cpu_seconds=sum(timings),
               mcs_timeouts=sum(r["timed_out"] for r in results), mcs_empty=sum(not r["timed_out"] and not r["mcs_atoms"] for r in results),
               truth_in_set_checks=sum(r["family"] == "truth_in_set" and not r["timed_out"] and bool(r["mcs_atoms"]) for r in results))
    checks = [r for r in results if r["family"] == "truth_in_set" and not r["timed_out"] and r["mcs_atoms"]]
    if not all(r["contained"] for r in checks):
        raise AssertionError("truth_in_set failed")


def scope_labels(row):
    for source in ("all", row["source"]):
        for lipid in ("all", str(row["lipid"])):
            for heavy_bin in ("all", row["heavy_bin"]):
                yield source, lipid, heavy_bin


def write_csv(path, rows):
    if not rows:
        raise ValueError("No aggregate rows")
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def curves(args, run):
    results = load(args.out / "results.local-vs21.pkl")
    grouped, chance_groups = collections.defaultdict(list), collections.defaultdict(list)
    for row in results:
        if row["motif_oversample"]:
            continue
        if row["chance_kind"]:
            for source, lipid, b in scope_labels(row):
                chance_groups[row["chance_kind"], source, lipid, b, row["k"]].append(row)
        if row["source"] != "db":
            for source, lipid, b in scope_labels(row):
                grouped[source, lipid, b, row["family"], row["k"]].append(row)
    output = []
    for (source, lipid, b, family, k), rows in sorted(grouped.items()):
        valid = [r for r in rows if not r["timed_out"] and r["mcs_atoms"]]
        timeout_count = sum(r["timed_out"] for r in rows)
        base = {"source": source, "lipid": lipid, "heavy_bin": b, "family": family, "k": k,
                "total_sets": len(rows), "timeouts": timeout_count,
                "empty": sum(not r["timed_out"] and not r["mcs_atoms"] for r in rows),
                "timeout_flag": timeout_count / len(rows) > 0.02}
        for axis, field in (("heavy_atoms", "mcs_atoms"), ("fraction_smallest", "fraction_smallest"), ("fraction_truth", "fraction_truth")):
            values = [(r[field], bool(r["contained"])) for r in valid]
            thresholds = list(range(1, 81)) if axis == "heavy_atoms" else [i / 20 for i in range(21)]
            for threshold in thresholds:
                output.append(dict(base, axis=axis, curve="cumulative", threshold=threshold, upper="", **cumulative(values, threshold)))
            intervals = ([(i, i + 1) for i in range(1, 31)] + [(31, 41), (41, 61), (61, 81)] if axis == "heavy_atoms" else
                         [(i / 20, (i + 1) / 20 if i < 20 else math.inf) for i in range(21)])
            for lower, upper in intervals:
                cell = [(v, ok) for v, ok in values if lower <= v < upper]
                output.append(dict(base, axis=axis, curve="point", threshold=lower,
                                   upper=upper if math.isfinite(upper) else "inf", **cumulative(cell, lower)))
    write_csv(args.out / "curves_correctness.csv", output)
    chance_output = []
    for (kind, source, lipid, b, k), rows in sorted(chance_groups.items()):
        valid = [r for r in rows if not r["timed_out"] and r["mcs_atoms"]]
        for axis, field in (("heavy_atoms", "mcs_atoms"), ("fraction_smallest", "fraction_smallest")):
            values = [r[field] for r in valid]
            chance_output.append({"kind": kind, "source": source, "lipid": lipid, "heavy_bin": b, "k": k,
                                  "axis": axis, "n": len(values), "total_sets": len(rows),
                                  "timeouts": sum(r["timed_out"] for r in rows),
                                  "empty": sum(not r["timed_out"] and not r["mcs_atoms"] for r in rows),
                                  "mean": sum(values) / len(values) if values else None,
                                  **{f"p{q}": quantile(values, q / 100) for q in (50, 90, 95, 99)},
                                  "shortfall_n200": max(0, 200 - len(values)),
                                  "timeout_flag": sum(r["timed_out"] for r in rows) / len(rows) > 0.02})
    write_csv(args.out / "curves_chance.csv", chance_output)
    motif_rows = []
    for motif_id in motifs():
        for family in FAMILIES:
            selected = [r for r in results if r["source"] != "db" and r["family"] == family and motif_id in r["motif_matches"]]
            n = len(selected)
            success = sum(motif_id in r["truth_motifs"] for r in selected)
            lo, hi = wilson(success, n)
            motif_rows.append({"motif_id": motif_id, "family": family, "n": n, "success": success,
                               "precision": success / n if n else None, "wilson_lo": lo, "wilson_hi": hi,
                               "motif_oversample_sets": sum(r["motif_oversample"] for r in selected),
                               "status": "not_estimable" if n < args.n_min else "pass" if lo >= args.target else "fail"})
    write_csv(args.out / "motif_precision.csv", motif_rows)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure_dir = args.out / "figures"
    figure_dir.mkdir(exist_ok=True)
    for axis in ("heavy_atoms", "fraction_smallest", "fraction_truth"):
        fig, panels = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
        for k, panel in zip((2, 3, 4, 5), panels.flat):
            for family in ("spectral_nn", "same_formula_db_random", "same_formula_db_similar", "spectral_nn_any"):
                rows = [r for r in output if r["source"] == r["lipid"] == r["heavy_bin"] == "all"
                        and r["family"] == family and r["k"] == k and r["axis"] == axis and r["curve"] == "cumulative" and r["n"]]
                panel.plot([r["threshold"] for r in rows], [r["precision"] for r in rows], label=family)
            chance = [r for r in chance_output if r["kind"] == "matched" and r["source"] == r["lipid"] == r["heavy_bin"] == "all"
                      and r["k"] == k and r["axis"] == axis]
            if chance and chance[0]["p95"] is not None:
                panel.axvline(chance[0]["p95"], linestyle=":", color="black", label="chance95")
            panel.axhline(args.target, linestyle="--", color="gray")
            panel.set(title=f"k = {k}", xlabel=axis, ylabel="P(contained | size >= x)", ylim=(0, 1.02))
        panels[0, 0].legend(fontsize=7)
        fig.savefig(figure_dir / f"correctness_{axis}.png", dpi=160)
        plt.close(fig)


def read_csv(path):
    with open(path, encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def choose(args, run):
    correct = read_csv(args.out / "curves_correctness.csv")
    chance = read_csv(args.out / "curves_chance.csv")
    def bars(scope):
        source, lipid, b = scope
        by_k = {}
        for k in (2, 3, 4, 5):
            c = [r for r in correct if (r["source"], r["lipid"], r["heavy_bin"]) == scope
                 and int(r["k"]) == k and r["curve"] == "cumulative"]
            ch = [r for r in chance if (r["source"], r["lipid"], r["heavy_bin"]) == scope
                  and int(r["k"]) == k and r["kind"] == "matched"]
            chances = {r["axis"]: float(r["p95"]) if r["p95"] else None for r in ch}
            a = choose_axis([r for r in c if r["axis"] == "heavy_atoms"], chances.get("heavy_atoms"), args.target, args.n_min)
            f = choose_axis([r for r in c if r["axis"] == "fraction_smallest"], chances.get("fraction_smallest"), args.target, args.n_min)
            by_k[str(k)] = {"heavy_atoms": int(a["value"]) if a else None, "fraction": f["value"] if f else None,
                            "status": "estimable" if a and f else "not_estimable",
                            "axis_status": {"heavy_atoms": "estimable" if a else "not_estimable", "fraction": "estimable" if f else "not_estimable"},
                            "n_at_bar": {"heavy_atoms": a["n"] if a else None, "fraction": f["n"] if f else None},
                            "wilson_lo_at_bar": {"heavy_atoms": a["wilson_lo"] if a else None, "fraction": f["wilson_lo"] if f else None},
                            "chance95_atoms": chances.get("heavy_atoms"), "chance95_fraction": chances.get("fraction_smallest"),
                            "deciding_family": {"heavy_atoms": a["deciding_family"] if a else None, "fraction": f["deciding_family"] if f else None},
                            "evidence": {"heavy_atoms": a, "fraction": f}}
        return by_k
    by_k = bars(("all", "all", "all"))
    global_bar = {"heavy_atoms": None, "fraction": None, "n_at_bar": {}, "wilson_lo_at_bar": {},
                  "chance95_atoms": None, "chance95_fraction": None, "deciding_family": {}, "deciding_k": {}}
    for axis in ("heavy_atoms", "fraction"):
        estimable = [(int(k), r) for k, r in by_k.items() if r[axis] is not None]
        if estimable:
            k, row = max(estimable, key=lambda item: item[1][axis])
            global_bar[axis] = row[axis]
            global_bar["n_at_bar"][axis] = row["n_at_bar"][axis]
            global_bar["wilson_lo_at_bar"][axis] = row["wilson_lo_at_bar"][axis]
            global_bar["deciding_family"][axis] = row["deciding_family"][axis]
            global_bar["deciding_k"][axis] = k
            chance_field = "chance95_atoms" if axis == "heavy_atoms" else "chance95_fraction"
            global_bar[chance_field] = row[chance_field]
    global_bar["status"] = "estimable" if all(global_bar[a] is not None for a in ("heavy_atoms", "fraction")) else "not_estimable"
    global_bar["not_estimable_k"] = [k for k, r in by_k.items() if r["status"] != "estimable"]
    scopes = sorted({(r["source"], r["lipid"], r["heavy_bin"]) for r in correct})
    # Marginal strata requested by the rule, plus fully crossed reporting cells.
    strata = {f"source={s};lipid={l};heavy_bin={b}": bars((s, l, b)) for s, l, b in scopes if (s, l, b) != ("all", "all", "all")}
    above = []
    for label, rows in strata.items():
        for k, row in rows.items():
            for axis in ("heavy_atoms", "fraction"):
                if row[axis] is not None and by_k[k][axis] is not None and row[axis] > by_k[k][axis]:
                    above.append({"stratum": label, "k": int(k), "axis": axis, "bar": row[axis], "pooled": by_k[k][axis]})
    motif_output = []
    motif_rows = read_csv(args.out / "motif_precision.csv")
    for motif_id in sorted({r["motif_id"] for r in motif_rows}):
        rows = [r for r in motif_rows if r["motif_id"] == motif_id and r["family"] in BAR_FAMILIES]
        status = "pass" if all(r["status"] == "pass" for r in rows) else "fail" if any(r["status"] == "fail" for r in rows) else "not_estimable"
        motif_output.append({"motif_id": motif_id, "status": status,
                             "n": {r["family"]: int(r["n"]) for r in rows},
                             "precision": {r["family"]: float(r["precision"]) if r["precision"] else None for r in rows},
                             "wilson_lo": {r["family"]: float(r["wilson_lo"]) for r in rows}})
    result = {"label": "ws8-p1-si-bar", "run_id": args.out.name, "target": args.target, "interval": "wilson95",
              "n_min": args.n_min, "chance_quantile": 0.95, "families_used_for_bar": list(BAR_FAMILIES),
              "mcs_params": mcs_params, "mcs_params_sha256": mcs_params_sha256,
              "motifs_file_sha256": sha256(HERE / "motifs-p1.tsv"), "by_k": by_k, "global": global_bar,
              "by_stratum": strata, "strata_above_pooled": above, "motifs": motif_output, "inputs": run["inputs"]}
    write_json(args.out / "si_bar.json", result)
    report = ["# WS8 P1 SI size experiment", "", f"Run: {args.out.name}", "",
              f"Global atoms: {global_bar['heavy_atoms']}; fraction: {global_bar['fraction']}; status: {global_bar['status']}.",
              "", "| k | atoms | fraction | status |", "|---|---|---|---|"]
    report.extend(f"| {k} | {r['heavy_atoms']} | {r['fraction']} | {r['status']} |" for k, r in by_k.items())
    report += ["", f"MCS sets: {run.get('requested_mcs')}; timeouts: {run.get('mcs_timeouts')}; empty: {run.get('mcs_empty')}.",
               f"Truth-in-set checks: {run.get('truth_in_set_checks')} (all contained).", "",
               "Named motifs: " + ", ".join(r["motif_id"] + ": " + r["status"] for r in motif_output), "",
               f"Strata above pooled: {len(above)}.", "", "Shortfalls (source, stratum, family, k):", "",
               "```json", json.dumps(run.get("shortfalls", {}), sort_keys=True), "```", "",
               "Chance cells with fewer than 200 valid sets are explicitly marked in curves_chance.csv.",
               "Timeouts and empty MCS are excluded from size curves. Motif oversamples affect motif precision alone.",
               "Global values are maxima over estimable k; unavailable k are explicitly listed."]
    (args.out / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def privacy_audit(out):
    audited = [out / name for name in ("curves_correctness.csv", "curves_chance.csv", "motif_precision.csv", "si_bar.json", "report.md")]
    for path in audited:
        text = path.read_text(encoding="utf-8")
        if re.search(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", text) or re.search(r"\b[A-Z]{14}\b", text):
            raise AssertionError("Identifier in aggregate artifact")
        if path.suffix == ".csv":
            columns = next(csv.reader(text.splitlines()))
            if any(re.search(r"key|smiles|name", c, re.I) for c in columns):
                raise AssertionError("Identifying aggregate column")
    return {"aggregate_files_checked": len(audited), "passed": True}


def finish(args, run):
    run["files"] = {str(p.relative_to(args.out)).replace("\\", "/"): sha256(p)
                    for p in sorted(args.out.rglob("*")) if p.is_file() and p.name not in ("run.json", "run.log")}
    write_json(args.out / "run.json", run)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "neighbours", "sets", "mcs", "curves", "choose", "all"))
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--sources", default="vs20,vs21")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20261113)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cap", type=int, choices=(300, 600), default=600)
    parser.add_argument("--target", type=float, default=0.95)
    parser.add_argument("--n-min", type=int, default=100)
    parser.add_argument("--pilot-run", type=Path)
    parser.add_argument("--reuse-prepare", type=Path)
    parser.add_argument("--launch", action="store_true", help="Detach from the app job with Windows breakaway")
    args = parser.parse_args(argv)
    args.sources = args.sources.split(",")
    if not 1 <= args.workers <= 16 or not set(args.sources) <= {"vs20", "vs21"}:
        parser.error("Use 1..16 workers and sources vs20[,vs21]")
    args.out = args.out.resolve()
    if args.out.parent != OUTPUT_ROOT or not re.fullmatch(r"p1(?:pilot)?-\d{8}-[0-9a-f]{7,40}", args.out.name):
        parser.error("Output must be the assigned root and P1 run-id format")
    args.out.mkdir(parents=True, exist_ok=True)
    for variable, directory in (("NUMBA_CACHE_DIR", "cache/numba"), ("MPLCONFIGDIR", "cache/matplotlib"),
                                ("TEMP", "cache/tmp"), ("TMP", "cache/tmp")):
        target = args.out / directory
        target.mkdir(parents=True, exist_ok=True)
        os.environ[variable] = str(target)
    if args.launch:
        command = [sys.executable, "-B", str(Path(__file__).resolve())] + [a for a in sys.argv[1:] if a != "--launch"]
        with open(args.out / "run.log", "ab", buffering=0) as log:
            process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                                       creationflags=0x01000000 | 0x00000008 | 0x08000000 | 0x00004000)
        print(json.dumps({"launched_pid": process.pid, "run_id": args.out.name, "job_breakaway": True}))
        return
    worker_init()
    import numpy as np
    import rdkit
    import matchms
    run_path = args.out / "run.json"
    run = json.loads(run_path.read_text(encoding="utf-8")) if run_path.exists() else {
        "label": "ws8-p1", "run_id": args.out.name, "code_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip(),
        "environment": {"interpreter": sys.executable, "python": sys.version, "rdkit": rdkit.__version__,
                        "matchms": matchms.__version__, "numpy": np.__version__},
        "parameters": {"pilot": args.pilot, "sources": args.sources, "workers": args.workers, "seed": args.seed,
                       "per_stratum_cap": args.cap, "target": args.target, "n_min": args.n_min,
                       "unconditional_per_cell": 50 if args.pilot else 1000, "motif_extra_cap": 8 if args.pilot else 150,
                       "spectral": {"class": "matchms.similarity.ModifiedCosineGreedy", "tolerance": 0.02,
                                    "mz_power": 0, "intensity_power": 1, "min_matches": 6, "score_floor": 0.6,
                                    "precursor_window": 200, "prep": "positive finite peaks; sorted m/z; raw intensity; precursor kept"}},
        "mcs_params": mcs_params, "mcs_params_sha256": mcs_params_sha256,
        "resource_control": {"workers": args.workers, "parent_priority": priority(), "worker_priority": "BelowNormal; verified by initializer",
                             "thread_environment": {v: os.environ[v] for v in THREAD_VARS}, "spawn": True, "job_breakaway": True},
        "timings": {}, "status": "running"}
    if not args.pilot and args.stage == "all":
        if args.pilot_run is None:
            parser.error("A full run requires --pilot-run")
        pilot = json.loads((args.pilot_run / "run.json").read_text(encoding="utf-8"))
        projected = pilot["projected_full_seconds_at_16_workers"]
        args.cap = 600 if projected <= 21600 else 300
        run["parameters"]["per_stratum_cap"] = args.cap
        run["pilot_projection"] = {"run_id": args.pilot_run.name, "seconds_at_16": projected,
                                   "cap_decision": args.cap, "reason": "within_6h" if args.cap == 600 else "projection_above_6h"}
    stages = ("prepare", "neighbours", "sets", "mcs", "curves", "choose") if args.stage == "all" else (args.stage,)
    finish(args, run)
    try:
        for stage in stages:
            print(f"stage {stage} started", flush=True)
            start = time.monotonic()
            if stage == "prepare" and args.reuse_prepare:
                import shutil
                previous = json.loads((args.reuse_prepare / "run.json").read_text(encoding="utf-8"))
                for path in args.reuse_prepare.iterdir():
                    if path.name.startswith(("prepared.", "libraries.", "truths.")):
                        shutil.copyfile(path, args.out / path.name)
                for key in ("inputs", "rejection_counts", "truth_counts", "representative_counts", "db_usable"):
                    run[key] = previous[key]
                run["prepare_reused_from"] = args.reuse_prepare.name
            else:
                globals()[stage](args, run)
            run["timings"][stage] = time.monotonic() - start
            finish(args, run)
            print(f"stage {stage} finished: {run['timings'][stage]:.1f}s", flush=True)
        if args.pilot and "mcs" in run["timings"]:
            run["projected_full_seconds_at_16_workers"] = sum(
                seconds * (20 if stage in ("neighbours", "sets", "mcs") else 1) *
                (args.workers / 16 if stage in ("neighbours", "mcs") else 1)
                for stage, seconds in run["timings"].items())
            run["projection_method"] = "Measured stages: 20x neighbour/set/MCS pilot work, prepare/report once; parallel stages scaled to 16 workers"
        if (args.out / "si_bar.json").exists():
            run["privacy_audit"] = privacy_audit(args.out)
        run["status"] = "complete"
        finish(args, run)
    except BaseException as exc:
        run["status"] = "failed"
        run["error_type"] = type(exc).__name__
        finish(args, run)
        # Never log exception text from chemistry or real-data parsing.
        print("Run failed: " + type(exc).__name__, flush=True)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
