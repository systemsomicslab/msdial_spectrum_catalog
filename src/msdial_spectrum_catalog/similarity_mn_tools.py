"""Deterministic parity fixtures and aggregate-only MSP smoke checks."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys
import time

from .parsers import read_msp
from . import similarity_mn as mn

FORMAT = "msdial-similarity-fixtures"
VS20_NAME = "MSMS-Public_all-neg-VS20.msp"
VS20_MD5 = "2cc3f3843d7c139e807611553b966202"
FIXTURE_ROOT = Path(__file__).resolve().parents[2] / "similarity-fixtures"


def prep_params(prep=mn.MSDIAL_MN_V1):
    return {"prep": "mn", **{"prep." + k: v for k, v in asdict(prep).items()}}


def kernel_ids(prep="none", tol=0.025):
    params = {"prep": "none"} if prep == "none" else prep_params(prep)
    return [mn.format_kernel_id(base, {**params, "tol": tol}) for base in
            ("msdial.moddot.legacy", "msdial.bonanza.v1", "msdial.cosine_matched.legacy")] + [
        mn.format_kernel_id("gnps.modified_cosine", {**params, "tol": tol,
                            "assign": assign, "mzpow": 0, "ipow": 1})
        for assign in ("greedy", "optimal")] + [
        mn.format_kernel_id("entropy.li2021", {**params, "bin": 0.05, "weighting": "none"})]


def smoke_ids():
    default = kernel_ids(mn.MSDIAL_MN_V1)
    gnps = kernel_ids(replace(mn.MSDIAL_MN_V1, frame=0, precwin=17, topk=6))
    return [default[0], default[1], default[3], gnps[3], kernel_ids()[5]]


def json_bytes(value):
    return (json.dumps(value, indent=1, ensure_ascii=True, allow_nan=False) + "\n").encode("utf-8")


def _envelope(kind, cases, **extra):
    return {"format": FORMAT, "format_version": 1, "kind": kind, **extra, "cases": cases}


def _spectrum(value):
    return mn.Spectrum(value["precursor_mz"], tuple(map(tuple, value["peaks"])))


def _input(precursor, peaks):
    return {"precursor_mz": precursor, "peaks": [list(p) for p in peaks]}


def _encoded_score(value):
    return "NaN" if value is not None and math.isnan(value) else value


def _pair_case(name, purpose, a, b, tol=0.025, swap=True, tags=()):
    return {"id": name, "purpose": purpose, "tags": list(tags),
            "implementations": ["python", "csharp"], "a": a, "b": b,
            "expected": _expectations(a, b, tol, swap)}


def _expectations(a, b, tol=0.025, swap=True):
    expected = []
    for identity in kernel_ids(tol=tol) + kernel_ids(mn.MSDIAL_MN_V1, tol):
        result = mn.score(identity, _spectrum(a), _spectrum(b))
        entry = {"kernel_id": identity, "score": _encoded_score(result.score),
                 "abs_tol": 1e-9, "swap_invariant": swap,
                 "implementations": ["python"] if "assign=optimal" in identity else ["python", "csharp"]}
        if "assign=optimal" not in identity:
            entry["matched_peak_count"] = result.matched_peak_count
        if identity.startswith("msdial."):
            entry["msdial_pairs"] = {"e": result.details["e"], "pairs": result.details["pairs"]}
        expected.append(entry)
    return expected


def synthetic_pairs():
    a = _input(250, [(100, 3), (150, 4)])
    cases = [
        _pair_case("identity", "Identical positive spectra", a, a),
        _pair_case("disjoint", "No product or neutral-loss overlap", a, _input(270, [(80, 6), (190, 7)])),
        _pair_case("moddot-single-shared-fragment", "Matched-only single-match defect", a,
                   _input(264, [(100, 3), (180, 12)]), tags=["defect:moddot-matched-only"]),
        _pair_case("precursor-only", "Retained precursor pair through shift", _input(200, [(200, 7)]),
                   _input(214, [(214, 9)])),
        _pair_case("heavy-unmatched", "Two matches with dominant unmatched intensity",
                   _input(250, [(80, 90), (100, 3), (150, 4)]),
                   _input(264, [(100, 3), (164, 4), (220, 99)])),
        _pair_case("shift-identity", "Product and neutral loss make full identity", a,
                   _input(264, [(100, 3), (164, 4)])),
        _pair_case("neutral-loss-beats-product", "Smallest intensity difference wins over product",
                   _input(200, [(90, 9), (100, 3)]), _input(210, [(100, 10)])),
        _pair_case("equal-intensity-difference", "First e-index wins a tied intensity difference",
                   _input(200, [(99.99, 8), (100.01, 12)]), _input(210, [(100, 10)])),
        _pair_case("equal-precursor-asymmetry", "Equal precursor makes e the second argument",
                   _input(200, [(100, 10), (100.02, 25)]), _input(200, [(100.01, 20)]), swap=False),
        _pair_case("small-shift-plain-cosine", "No shifted candidates when shift <= tol",
                   _input(200, [(100, 10)]), _input(200.02, [(100.04, 10)])),
        _pair_case("inclusive-boundary", "Exactly representable tolerance boundary",
                   _input(200, [(100, 10)]), _input(200, [(100.25, 10)]), tol=0.25),
        _pair_case("greedy-below-optimal", "One large edge blocks two useful edges",
                   _input(200, [(100, 9), (100.04, 10)]),
                   _input(200, [(100.02, 10), (100.06, 9)]), swap=False),
        _pair_case("moddot-zero-nan", "Legacy zero denominator is NaN; All has a guard",
                   _input(200, [(100, 0)]), _input(214, [(100, 0)])),
        _pair_case("entropy-empty", "Empty spectrum is not comparable", _input(200, []), a),
        _pair_case("entropy-disjoint", "Distinct fixed bins", _input(200, [(90, 1)]),
                   _input(200, [(110, 1)])),
        _pair_case("entropy-mixture", "Half shared distribution has intermediate entropy",
                   _input(200, [(90, 1), (110, 1)]), _input(200, [(90, 1), (130, 1)])),
        _pair_case("empty-both", "All kernels handle empty inputs", _input(200, []), _input(200, [])),
        _pair_case("stable-equal-mass", "Equal masses retain input order",
                   _input(200, [(100, 8), (100, 12)]), _input(210, [(100, 10)])),
    ]
    return _envelope("kernel_pairs", cases)


def prep_cases():
    default = mn.MSDIAL_MN_V1
    inputs = [
        ("base-before-filter", "Above-precursor maximum controls relative cutoff", 200,
         [(100, 5), (150, 10), (250, 1000)], default),
        ("above-precursor-boundary", "Just below, equal and above P+0.05", 200,
         [(200.04999, 4), (200.05, 4), (200.05001, 4)], replace(default, frame=0)),
        ("absolute-cutoff", "Absolute cutoff is inclusive", 200,
         [(90, 4.99999), (100, 5), (110, 6)], replace(default, abs=5)),
        ("relative-equality", "Keep intensity exactly at one percent", 200,
         [(90, 0.99999), (100, 1), (150, 100)], default),
        ("frame-summing", "Sum in m/z order and retain intensity apex", 200,
         [(100.1, 3), (100.2, 7), (101.1, 5)], default),
        ("frame-apex-tie", "First strict greater apex means lower m/z", 200,
         [(100.1, 5), (100.9, 5)], default),
        ("unsorted", "Stable m/z ordering before framing", 200,
         [(101.3, 7), (100.8, 5), (100.1, 5), (100.1, 5)], default),
        ("empty", "Empty stays empty", 200, [], default),
        ("precursor-retained", "Default MN keeps the precursor", 200, [(200, 1)], default),
        ("power-scale", "Explicit intensity power and scale", 200,
         [(100, 2), (110, 4)], replace(default, pow=2, scale=10)),
        ("precursor-window", "Inclusive Python precursor-window removal", 200,
         [(182.99999, 1), (183, 1), (200, 1), (217, 1)],
         replace(default, above=None, precwin=17)),
        ("topk", "Sliding top-k uses strictly greater intensities and inclusive window", 200,
         [(50, 10), (75, 10), (100, 5), (100.00001, 4), (175, 3)],
         replace(default, topk=1, topkwin=50)),
        ("mass-range", "Keep min and max equality", 200,
         [(9.99999, 1), (10, 2), (150, 3), (150.00001, 4)],
         replace(default, minmz=10, maxmz=150)),
        ("frame-off", "No framing preserves order at the low-level API", 200,
         [(110, 3), (100, 4)], replace(default, frame=0)),
    ]
    cases = []
    for name, purpose, precursor, peaks, prep in inputs:
        identity = mn.format_kernel_id("msdial.bonanza.v1", {**prep_params(prep), "tol": 0.025})
        prep_id = identity.split("?", 1)[1].rsplit("&tol=", 1)[0]
        cases.append({"id": name, "purpose": purpose,
                      "implementations": ["python"] if prep.precwin or prep.topk else ["python", "csharp"],
                      "input": _input(precursor, peaks), "prep_id": prep_id,
                      "expected_peaks": [list(p) for p in mn.preprocess(mn.Spectrum(precursor, peaks), prep).peaks],
                      "abs_tol": 1e-9})
    return _envelope("prep", cases)


def id_cases():
    values = ["msdial.moddot.legacy", "msdial.bonanza.v1", "entropy.li2021",
              kernel_ids()[3], "gnps.modified_cosine",
              "msdial.bonanza.v1?prep=none", "msdial.bonanza.v1?prep=none&tol=0.025&unknown=1",
              "msdial.bonanza.v1?tol=0.0250&prep=none",
              "gnps.modified_cosine?prep=none&assign=greedy&ipow=1.0&mzpow=-0&tol=0.0250",
              "msdial.cosine_matched.legacy", "entropy.li2021?bin=0.05&prep=none&weighting=true",
              "msdial.bonanza.v1?prep=none&tol=0.025&tol=0.025",
              "msdial.bonanza.v1?prep=none&tol=NaN",
              "msdial.bonanza.v1?prep=none&tol=0",
              "unknown?prep=none&tol=1", "msdial.bonanza.v1?prep=none&tol=1e-2"]
    cases = []
    for index, value in enumerate(values):
        try:
            canonical = mn.canonical_kernel_id(value)
        except ValueError:
            canonical = "error"
        cases.append({"id": "id-%02d" % index, "input": value, "canonical": canonical})
    return _envelope("kernel_ids", cases)


README = """# Similarity parity fixtures, format version 1

These synthetic and public fixtures are shared unchanged with WS3 C# tests.
Each JSON case file has format, format_version, kind and cases. The manifest
pins the SHA256 and case count of all four case files (not itself or README).
UTF-8, ensure_ascii=True, indent=1, LF, trailing newline; no timestamps.

Run `python -m msdial_spectrum_catalog.similarity_mn_tools fixtures-check
similarity-fixtures`. `fixtures-write DIR` rebuilds synthetic definitions and
all expected values. Public inputs come from the committed kernels-vs20.json;
no library is needed to regenerate. For initial selection only, use
`fixtures-write DIR --vs20-msp PATH`, which validates the adjacent
zenodo-record.json and the full file MD5 before reading records.

Public source: MSMS-Public_all-neg-VS20.msp, MD5
2cc3f3843d7c139e807611553b966202. Attribution: MS-DIAL public mass spectral
library, Zenodo doi:10.5281/zenodo.21904103,
https://doi.org/10.5281/zenodo.21904103, CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/). The selection rule is the first
16 records in file order having a finite positive PRECURSORMZ and at least
three finite nonnegative-m/z, positive-intensity peaks; pair consecutive
selected records into eight pairs. Record indices are zero-based and stored
in each case. Only precursor and peaks are retained. Precursors and m/z are
rounded to five decimal places; intensities to six significant digits before
any computation. Source indices and normalized inputs make the selection
auditable. Source records remain otherwise unchanged in place.

Inputs are raw. Apply the prep in each expected kernel_id before scoring;
prep=none sorts stably by m/z only. msdial_pairs indices refer to these
processed lists; e says which argument has the e role. Compare scores with
abs_tol (default 1e-9), null as not comparable, and string NaN as IEEE NaN.
Compare matched counts only when present: optimal assignment can tie.
Implementations are attached to each expectation; optimal and prep cases
using precwin/topk are Python-only until C# supports them. swap_invariant
false permits the documented equal-precursor MS-DIAL role asymmetry.
prep.json uses a prep_id query without a base or tol. kernel-ids.json maps
input to canonical or error. All kernel identities are used verbatim across
languages. Cases with zero intensities also pin legacy NaN behavior.

The synthetic fixture definitions and generator are LGPL-3.0-or-later, as is
the repository; public VS20 inputs retain CC BY 4.0. No experimental record
is a fixture. See docs/similarity_kernels.md for formulas and source lines.
"""


def public_pairs(path=None):
    if path is None:
        existing = json.loads((FIXTURE_ROOT / "kernels-vs20.json").read_text(encoding="utf-8"))
        for case in existing["cases"]:
            case["expected"] = _expectations(case["a"], case["b"])
        return existing
    path = Path(path)
    metadata = json.loads(path.with_name("zenodo-record.json").read_text(encoding="utf-8"))
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    if path.name != VS20_NAME or metadata.get("md5") != VS20_MD5 or digest.hexdigest() != VS20_MD5:
        raise ValueError("VS20 provenance/MD5 mismatch; selection stopped")
    selected = []
    for record in read_msp(path):
        try:
            precursor = float(record.fields["PRECURSORMZ"])
        except (KeyError, ValueError):
            continue
        peaks = [(round(mz, 5), float(format(intensity, ".6g"))) for mz, intensity in record.peaks
                 if math.isfinite(mz) and math.isfinite(intensity) and mz >= 0 and intensity > 0]
        if math.isfinite(precursor) and precursor > 0 and len(peaks) >= 3:
            selected.append((record.index, _input(round(precursor, 5), peaks)))
        if len(selected) == 16:
            break
    if len(selected) != 16:
        raise ValueError("public source has fewer than 16 qualifying records")
    cases = []
    for index in range(0, 16, 2):
        (ai, a), (bi, b) = selected[index:index + 2]
        case = _pair_case("vs20-%02d" % (index // 2), "Consecutive qualifying public records", a, b)
        case["source_record_indices"] = [ai, bi]
        cases.append(case)
    return _envelope("kernel_pairs", cases, source={"filename": VS20_NAME, "md5": VS20_MD5,
                     "doi": "10.5281/zenodo.21904103", "license": "CC BY 4.0"})


def write_fixtures(directory, vs20_msp=None):
    files = {"kernel-ids.json": id_cases(), "prep.json": prep_cases(),
             "kernels.json": synthetic_pairs(), "kernels-vs20.json": public_pairs(vs20_msp)}
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {"format": FORMAT, "format_version": 1, "kernels_version": mn.KERNELS_VERSION, "files": []}
    for name, value in files.items():
        data = json_bytes(value)
        (directory / name).write_bytes(data)
        manifest["files"].append({"path": name, "sha256": hashlib.sha256(data).hexdigest(),
                                  "cases": len(value["cases"])})
    (directory / "manifest.json").write_bytes(json_bytes(manifest))
    (directory / "README.md").write_bytes(README.encode("utf-8"))
    return manifest


def _equal(actual, expected, tolerance):
    if expected == "NaN":
        return actual is not None and math.isnan(actual)
    if actual is None or expected is None:
        return actual is expected
    return math.isfinite(actual) and abs(actual - expected) <= tolerance


def check_fixtures(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("format"), manifest.get("format_version"), manifest.get("kernels_version")) != (FORMAT, 1, mn.KERNELS_VERSION):
        raise ValueError("invalid manifest version")
    names = [entry["path"] for entry in manifest["files"]]
    if set(names) != {"kernel-ids.json", "prep.json", "kernels.json", "kernels-vs20.json"} or len(names) != 4:
        raise ValueError("invalid manifest file set")
    total = expectations = 0
    for entry in manifest["files"]:
        data = (directory / entry["path"]).read_bytes()
        if hashlib.sha256(data).hexdigest() != entry["sha256"]:
            raise ValueError("manifest drift: " + entry["path"])
        value = json.loads(data)
        if value["format"] != FORMAT or value["format_version"] != 1 or len(value["cases"]) != entry["cases"]:
            raise ValueError("case format/count drift")
        for case in value["cases"]:
            total += 1
            label = entry["path"] + ":" + case["id"]
            if value["kind"] == "kernel_ids":
                try:
                    actual = mn.canonical_kernel_id(case["input"])
                except ValueError:
                    actual = "error"
                if actual != case["canonical"]:
                    raise ValueError(label + " identity drift")
            elif value["kind"] == "prep":
                _, params = mn.parse_kernel_id("msdial.bonanza.v1?" + case["prep_id"] + "&tol=0.025")
                prep = mn.MnPrep(**{k[5:]: v for k, v in params.items() if k.startswith("prep.")})
                actual = mn.preprocess(_spectrum(case["input"]), prep).peaks
                expected = case["expected_peaks"]
                if len(actual) != len(expected) or any(not _equal(x, y, case["abs_tol"])
                        for ap, ep in zip(actual, expected) for x, y in zip(ap, ep)):
                    raise ValueError(label + " preprocessing drift")
            elif value["kind"] == "kernel_pairs":
                for exp in case["expected"]:
                    expectations += 1
                    result = mn.score(exp["kernel_id"], _spectrum(case["a"]), _spectrum(case["b"]))
                    if result.kernel_id != exp["kernel_id"] or not _equal(result.score, exp["score"], exp.get("abs_tol", 1e-9)):
                        raise ValueError(label + " score drift: " + exp["kernel_id"])
                    if "matched_peak_count" in exp and result.matched_peak_count != exp["matched_peak_count"]:
                        raise ValueError(label + " count drift")
                    if "msdial_pairs" in exp and exp["msdial_pairs"] != {"e": result.details["e"], "pairs": result.details["pairs"]}:
                        raise ValueError(label + " pair drift")
                    if exp["swap_invariant"]:
                        reverse = mn.score(exp["kernel_id"], _spectrum(case["b"]), _spectrum(case["a"]))
                        if not _equal(reverse.score, exp["score"], exp.get("abs_tol", 1e-9)):
                            raise ValueError(label + " swap drift")
            else:
                raise ValueError("unknown fixture kind")
    return {"ok": True, "files": len(names), "cases": total, "kernel_expectations": expectations}


_PEAKID = re.compile(r"(?:^|\|)PEAKID=([^|]*)", re.IGNORECASE)


def read_smoke_msp(path):
    """Yield (Spectrum, PEAKID); never retain any other MSP field or comment."""
    precursor, peaks, peakid = 0.0, [], None
    active = reading_peaks = False
    with Path(path).open("r", encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                if active:
                    yield mn.Spectrum(precursor, peaks), peakid
                precursor, peaks, peakid = 0.0, [], None
                active = reading_peaks = False
                continue
            active = True
            if reading_peaks:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz, intensity = float(parts[0]), float(parts[1])
                    except ValueError:
                        continue
                    if math.isfinite(mz) and math.isfinite(intensity):
                        peaks.append((mz, intensity))
                continue
            key, separator, value = line.partition(":")
            if not separator:
                continue
            key = key.strip().upper()
            if key == "PRECURSORMZ":
                precursor = float(value.strip())
                if not math.isfinite(precursor):
                    raise ValueError("nonfinite precursor")
            elif key == "COMMENT":
                match = _PEAKID.search(value.strip())
                if match:
                    peakid = match[1]
            elif key.replace(" ", "") == "NUMPEAKS":
                reading_peaks = True
    if active:
        yield mn.Spectrum(precursor, peaks), peakid


def _quantile(values, fraction):
    if not values:
        return None
    position = (len(values) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    return values[lower] + (values[upper] - values[lower]) * (position - lower)


def smoke(path, max_dprec=400):
    start = time.perf_counter()
    records = list(read_smoke_msp(path))
    spectra = [spectrum for spectrum, _ in records if spectrum.peaks]
    ids = smoke_ids()
    mn_spectra = [mn.preprocess(s) for s in spectra]
    gnps_spectra = [mn.preprocess(s, replace(mn.MSDIAL_MN_V1, frame=0, precwin=17, topk=6)) for s in spectra]
    values, one = [[] for _ in ids], [0] * len(ids)
    count = defect = 0
    for i, j in itertools.combinations(range(len(spectra)), 2):
        if abs(spectra[i].precursor_mz - spectra[j].precursor_mz) > max_dprec:
            continue
        count += 1
        a, b = mn_spectra[i], mn_spectra[j]
        results = [mn.msdial_moddot_legacy(a, b), mn.msdial_bonanza_v1(a, b),
                   mn.modified_cosine(a, b, tol=0.025, assign="greedy"),
                   mn.modified_cosine(gnps_spectra[i], gnps_spectra[j], tol=0.025, assign="greedy"),
                   mn.entropy_li2021(spectra[i], spectra[j])]
        for k, result in enumerate(results):
            values[k].append(result.score)
            one[k] += result.matched_peak_count == 1
        defect += results[0].score >= 0.95 and results[2].score < 0.5
    summaries = []
    for identity, raw, matched_one in zip(ids, values, one):
        finite = sorted(x for x in raw if x is not None and math.isfinite(x))
        n = len(finite)
        summaries.append({"kernel_id": identity, "finite_scores": n,
            "not_comparable": sum(x is None for x in raw),
            "nonfinite_scores": sum(x is not None and not math.isfinite(x) for x in raw),
            **{key: _quantile(finite, q) for key, q in
               [("min", 0), ("p05", 0.05), ("p25", 0.25), ("p50", 0.5), ("p75", 0.75), ("p95", 0.95), ("max", 1)]},
            "fraction_ge_0_7": sum(x >= 0.7 for x in finite) / n if n else None,
            "fraction_ge_0_9": sum(x >= 0.9 for x in finite) / n if n else None,
            "fraction_matched_1": matched_one / count if count else None})
    return {"records": len(records), "records_with_peaks": len(spectra), "pair_count": count,
            "max_dprec": max_dprec, "kernels": summaries,
            "fraction_k1_ge_0_95_and_k3_lt_0_5": defect / count if count else None,
            "runtime_s": time.perf_counter() - start}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    writer = sub.add_parser("fixtures-write")
    writer.add_argument("directory", type=Path)
    writer.add_argument("--vs20-msp", type=Path)
    checker = sub.add_parser("fixtures-check")
    checker.add_argument("directory", type=Path)
    runner = sub.add_parser("smoke")
    runner.add_argument("--msp", type=Path, required=True)
    runner.add_argument("--max-dprec", type=float, default=400)
    args = parser.parse_args(argv)
    try:
        if args.command == "fixtures-write":
            result = write_fixtures(args.directory, args.vs20_msp)
        elif args.command == "fixtures-check":
            result = check_fixtures(args.directory)
        else:
            if not math.isfinite(args.max_dprec) or args.max_dprec < 0:
                raise ValueError("max-dprec must be finite and nonnegative")
            result = smoke(args.msp, args.max_dprec)
    except (ValueError, OSError, KeyError) as exc:
        # Never echo input record contents or paths from the smoke command.
        print("smoke input/check failed" if args.command == "smoke" else str(exc), file=sys.stderr)
        return 1
    print(json.dumps(result, ensure_ascii=True, allow_nan=False, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
