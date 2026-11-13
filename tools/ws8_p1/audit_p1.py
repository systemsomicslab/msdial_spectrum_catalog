"""Audit completed P1 artifacts without repeating spectral scoring or MCS work.

All reports contain aggregates only. Run with the D5 interpreter, from the repo root.
The measurement run.json and si_bar.json remain unchanged; finalization.json records
the audit code identity and hashes of the additional artifacts.
"""
import argparse
from bisect import bisect_left
from collections import defaultdict
import json
import math
from pathlib import Path
import subprocess

from tools.ws8_p1 import p1_si_bar as p1


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def close(actual, expected):
    if expected is None:
        return actual in (None, "")
    if isinstance(expected, bool):
        return str(actual) == str(expected)
    return math.isclose(float(actual), expected, rel_tol=1e-10, abs_tol=1e-12)


def distribution(rows, field):
    pairs = sorted((r[field], bool(r["contained"])) for r in rows
                   if not r["timed_out"] and r["mcs_atoms"])
    sizes, successes = [], [0]
    for value, success in pairs:
        sizes.append(value)
        successes.append(successes[-1] + success)
    return sizes, successes


def audit_curves(out, measured):
    groups, chances = defaultdict(list), defaultdict(list)
    for row in measured:
        if row["motif_oversample"]:
            continue
        for scope in p1.scope_labels(row):
            if row["source"] != "db":
                groups[(*scope, row["family"], row["k"])].append(row)
            if row["chance_kind"]:
                chances[(row["chance_kind"], *scope, row["k"])].append(row)
    fields = {"heavy_atoms": "mcs_atoms", "fraction_smallest": "fraction_smallest", "fraction_truth": "fraction_truth"}
    stats = {}
    correct = p1.read_csv(out / "curves_correctness.csv")
    for row in correct:
        key = (row["source"], row["lipid"], row["heavy_bin"], row["family"], int(row["k"]))
        cache_key = (*key, row["axis"])
        members = groups[key]
        if cache_key not in stats:
            sizes, successes = distribution(members, fields[row["axis"]])
            stats[cache_key] = (sizes, successes, sum(r["timed_out"] for r in members),
                                sum(not r["timed_out"] and not r["mcs_atoms"] for r in members))
        sizes, successes, timeouts, empty = stats[cache_key]
        start = bisect_left(sizes, float(row["threshold"]))
        stop = len(sizes) if row["curve"] == "cumulative" else bisect_left(sizes, float(row["upper"]))
        n, success = stop - start, successes[stop] - successes[start]
        lo, hi = p1.wilson(success, n)
        expected = dict(n=n, success=success, precision=success / n if n else None,
                        wilson_lo=lo, wilson_hi=hi, total_sets=len(members), timeouts=timeouts,
                        empty=empty, timeout_flag=timeouts / len(members) > 0.02)
        check(all(close(row[k], v) for k, v in expected.items()), "Correctness CSV disagrees with measured rows")
    chance = p1.read_csv(out / "curves_chance.csv")
    for row in chance:
        members = chances[(row["kind"], row["source"], row["lipid"], row["heavy_bin"], int(row["k"]))]
        values = sorted(r[fields[row["axis"]]] for r in members if not r["timed_out"] and r["mcs_atoms"])
        timeouts = sum(r["timed_out"] for r in members)
        expected = dict(n=len(values), total_sets=len(members), timeouts=timeouts,
                        empty=sum(not r["timed_out"] and not r["mcs_atoms"] for r in members),
                        mean=sum(values) / len(values) if values else None,
                        shortfall_n200=max(0, 200 - len(values)), timeout_flag=timeouts / len(members) > 0.02,
                        **{f"p{q}": p1.quantile(values, q / 100) for q in (50, 90, 95, 99)})
        check(all(close(row[k], v) for k, v in expected.items()), "Chance CSV disagrees with measured rows")
    return correct, chance


def audit_bars(out, correct, chance):
    result = json.loads((out / "si_bar.json").read_text(encoding="utf-8"))
    required = {"label", "run_id", "target", "interval", "n_min", "chance_quantile", "families_used_for_bar",
                "mcs_params", "mcs_params_sha256", "motifs_file_sha256", "by_k", "global", "by_stratum",
                "strata_above_pooled", "motifs", "inputs"}
    check(required <= result.keys(), "Missing SI output field")
    curves = defaultdict(list)
    for row in correct:
        if row["curve"] == "cumulative":
            curves[row["source"], row["lipid"], row["heavy_bin"], row["k"], row["axis"]].append(row)
    diagnostic = []
    scopes = {("all", "all", "all"): result["by_k"]}
    for label, bars in result["by_stratum"].items():
        parts = dict(item.split("=") for item in label.split(";"))
        scopes[parts["source"], parts["lipid"], parts["heavy_bin"]] = bars
    for scope, bars in scopes.items():
        for k, bar in bars.items():
            for axis, output_axis, chance_field in (("heavy_atoms", "heavy_atoms", "chance95_atoms"),
                                                   ("fraction_smallest", "fraction", "chance95_fraction")):
                rows = curves[(*scope, k, axis)]
                q = bar[chance_field]
                passing = []
                for family in p1.BAR_FAMILIES:
                    eligible = [r for r in rows if r["family"] == family and q is not None
                                and float(r["threshold"]) > q and int(r["n"]) >= result["n_min"]]
                    passing.append({float(r["threshold"]) for r in eligible if float(r["wilson_lo"]) >= result["target"]})
                    if scope == ("all", "all", "all"):
                        best = max(eligible, key=lambda r: float(r["wilson_lo"]), default=None)
                        diagnostic.append(dict(k=int(k), axis=axis, family=family, chance95=q,
                                               first_passing=min(passing[-1]) if passing[-1] else None,
                                               best_eligible={c: float(best[c]) for c in ("threshold", "n", "precision", "wilson_lo")} if best else None))
                common = set.intersection(*passing)
                expected = min(common) if common else None
                check(bar[output_axis] == expected, "SI bar disagrees with independently searched written curves")
                check(bar["axis_status"][output_axis] == ("estimable" if expected is not None else "not_estimable"), "Incorrect SI axis status")
                chance_row = [r for r in chance if (r["source"], r["lipid"], r["heavy_bin"], r["k"], r["axis"], r["kind"]) == (*scope, k, axis, "matched")]
                check(close(chance_row[0]["p95"] if chance_row else None, q), "SI chance95 mismatch")
    for axis in ("heavy_atoms", "fraction"):
        values = [r[axis] for r in result["by_k"].values() if r[axis] is not None]
        check(result["global"][axis] == (max(values) if values else None), "Global SI value mismatch")
    return result, diagnostic


def compare_pilots(pilot, rerun):
    metadata = [json.loads((p / "run.json").read_text(encoding="utf-8")) for p in (pilot, rerun)]
    check(metadata[0]["inputs"] == metadata[1]["inputs"] and metadata[0]["mcs_params"] == metadata[1]["mcs_params"]
          and metadata[0]["parameters"]["seed"] == metadata[1]["parameters"]["seed"], "Pilot input or seed mismatch")
    names = ("truths.parquet", "truths.local-vs21.parquet", "neighbours.local-vs21.pkl", "pending.local-vs21.pkl",
             "sets.parquet", "sets.local-vs21.parquet", "chance_sets.parquet", "motif_precision.csv")
    identical = {name: p1.sha256(pilot / name) == p1.sha256(rerun / name) for name in names}
    check(all(identical[n] for n in names[:4]), "Pilot sampling or construction is not reproducible")
    a, b = p1.load(pilot / "results.local-vs21.pkl"), p1.load(rerun / "results.local-vs21.pkl")
    check(len(a) == len(b), "Pilot row counts differ")
    changes = [(x, y) for x, y in zip(a, b) if x != y]
    check(all(x["timed_out"] != y["timed_out"] for x, y in changes), "Non-timeout pilot measurement difference")
    identity = ("source", "family", "k", "truth_key", "member_keys")
    check(all(all(x[f] == y[f] for f in identity) for x, y in zip(a, b)), "Pilot row identities differ")
    # Independently rebuild the complete aggregate difference listing.
    differences = []
    labels = ("kind", "source", "family", "k", "lipid", "heavy_bin", "axis", "curve", "threshold", "upper")
    for name in ("curves_correctness.csv", "curves_chance.csv", "motif_precision.csv"):
        left, right = p1.read_csv(pilot / name), p1.read_csv(rerun / name)
        check(len(left) == len(right), "Pilot aggregate row counts differ")
        for x, y in zip(left, right):
            check(all(x.get(f) == y.get(f) for f in labels), "Pilot aggregate coordinates differ")
            for field in x:
                if x[field] != y[field]:
                    differences.append(dict(artifact=name, **{f: x.get(f, "") for f in labels},
                                            metric=field, pilot=x[field], rerun=y[field]))
    return dict(pilot=pilot.name, rerun=rerun.name, artifact_identical=identical,
                set_count=len(a), changed_rows=len(changes), non_timeout_changes=0,
                timeout_counts=dict(pilot=sum(r["timed_out"] for r in a), rerun=sum(r["timed_out"] for r in b)),
                strict_measured_set_identity=all(identical[n] for n in ("sets.parquet", "sets.local-vs21.parquet")),
                status="timeout_exception" if changes else "pass", curve_difference_metrics=len(differences)), differences


def compare_full_to_pilot(measured, pilot):
    def identity(row):
        return row["source"], row["family"], row["k"], row["truth_key"], tuple(row["member_keys"])
    full = {identity(r): r for r in measured}
    common, timeout_changes, completed_changes = 0, 0, 0
    fields = ("mcs_smarts", "mcs_atoms", "contained", "fraction_smallest", "fraction_truth")
    for row in p1.load(pilot / "results.local-vs21.pkl"):
        other = full.get(identity(row))
        if other is None:
            continue
        common += 1
        timeout_changes += row["timed_out"] != other["timed_out"]
        if not row["timed_out"] and not other["timed_out"]:
            completed_changes += any(row[f] != other[f] for f in fields)
    check(not completed_changes, "Common completed pilot/full measurements differ")
    return dict(common_constructed_sets=common, timeout_flips=timeout_changes,
                jointly_completed_measurement_changes=completed_changes,
                sampling_note="Pilot and full use different sample sizes; only identical constructed sets are compared.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--rerun", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    check(out.parent == p1.OUTPUT_ROOT and out.is_dir(), "Audit requires an existing assigned run directory")
    p1.priority()
    run = json.loads((out / "run.json").read_text(encoding="utf-8"))
    check(run["status"] == "complete", "Measurement run is incomplete")
    check(p1.sha256(p1.RESOURCE / "MsfinderStructureDB-VS15.esd") == p1.DB_SHA, "Structure DB checksum mismatch")
    for name, digest in run["files"].items():
        check(p1.sha256(out / name) == digest, "Measurement artifact checksum mismatch")
    measured = p1.load(out / "results.local-vs21.pkl")
    check(len(measured) == run["requested_mcs"], "Measurement row count mismatch")
    import pandas as pd
    truth_keys = {}
    truth_counts = {}
    for source in ("vs20", "vs21"):
        truths = pd.read_parquet(out / p1.row_filename("truths", source))
        check({"key", "formula", "atoms", "lipid", "lipid_source", "direct_parent", "stratum"} <= set(truths.columns), "Truth schema incomplete")
        check(not truths["key"].duplicated().any() and set(truths["source"]) == {source}, "Truth source split or uniqueness failed")
        truth_keys[source] = set(truths["key"])
        truth_counts[source] = dict(base=int((~truths["motif_oversample"]).sum()), motif_extra=int(truths["motif_oversample"].sum()))
    check(not truth_keys["vs20"] & truth_keys["vs21"], "Truth sources overlap")
    for source in ("vs20", "vs21", "db"):
        name = p1.row_filename("sets" if source != "db" else "chance_sets", source)
        expected = pd.DataFrame([r for r in measured if r["source"] == source])
        try:
            pd.testing.assert_frame_equal(pd.read_parquet(out / name), expected, check_dtype=False)
        except AssertionError:
            raise AssertionError("Measured parquet disagrees with saved results") from None
    leaks = sum(r["truth_key"] in r["member_keys"] for r in measured if r["truth_key"] and r["family"] != "truth_in_set")
    check(leaks == 0, "Truth leaked into candidate sets")
    controls = [r for r in measured if r["family"] == "truth_in_set" and not r["timed_out"] and r["mcs_atoms"]]
    check(all(r["contained"] for r in controls), "Truth-in-set containment failed")
    check(all(not r["mcs_atoms"] and not r["mcs_smarts"] and r["contained"] is None
              for r in measured if r["timed_out"]), "A partial timeout result was used")
    correct, chance = audit_curves(out, measured)
    result, diagnostic = audit_bars(out, correct, chance)
    check(result["inputs"] == run["inputs"] and result["mcs_params_sha256"] == run["mcs_params_sha256"] == p1.mcs_params_sha256
          and result["motifs_file_sha256"] == p1.sha256(p1.HERE / "motifs-p1.tsv"), "SI provenance mismatch")
    for row in p1.read_csv(out / "motif_precision.csv"):
        selected = [r for r in measured if r["source"] != "db" and r["family"] == row["family"] and row["motif_id"] in r["motif_matches"]]
        n = len(selected)
        success = sum(row["motif_id"] in r["truth_motifs"] for r in selected)
        check(int(row["n"]) == n and int(row["success"]) == success and close(row["wilson_lo"], p1.wilson(success, n)[0]), "Motif CSV mismatch")
    privacy = p1.privacy_audit(out)
    check(all("path" not in r for r in run["inputs"]), "Input manifest contains a path")
    resource = run["resource_control"]
    check(resource["workers"] <= 16 and resource["parent_priority"] == "BelowNormal" and resource["job_breakaway"]
          and all(v == "1" for v in resource["thread_environment"].values()), "Resource contract failed")
    shortfalls = [{k: r[k] for k in ("kind", "source", "lipid", "heavy_bin", "k", "axis", "n", "shortfall_n200")}
                  for r in chance if int(r["shortfall_n200"]) > 0]
    pilot_comparison, differences = compare_pilots(args.pilot, args.rerun)
    full_comparison = compare_full_to_pilot(measured, args.pilot)
    p1.write_json(out / "pilot-reproducibility.json", pilot_comparison)
    if differences:
        p1.write_csv(out / "pilot-curve-differences.csv", differences)
    audit = dict(status="pass_with_A4_timeout_exception", sets=len(measured), truth_leaks=leaks,
                 truth_counts=truth_counts,
                 truth_in_set_valid=len(controls), timeouts=run["mcs_timeouts"], empty=run["mcs_empty"],
                 correctness_rows_checked=len(correct), chance_rows_checked=len(chance),
                 chance_shortfalls=shortfalls, privacy=privacy, bar_diagnostics=diagnostic,
                 measurement_artifacts_verified=len(run["files"]), pilot_reproducibility=pilot_comparison,
                 full_vs_pilot=full_comparison)
    p1.write_json(out / "acceptance-audit.json", audit)
    digest = p1.sha256(out / "si_bar.json")
    def entry(value, key):
        return dict(value=value, source=dict(pre_experiment="WS8 P1", run_id=out.name, sha256=digest, key=key))
    rules = dict(label="claim-rules-v1", status="provisional", si={
        "mcs_min_heavy_atoms": entry(result["global"]["heavy_atoms"], "global.heavy_atoms"),
        "mcs_min_fraction": entry(result["global"]["fraction"], "global.fraction"),
        "mcs_bar_by_k": entry(result["by_k"], "by_k"),
        "motifs_allowed": entry([r["motif_id"] for r in result["motifs"] if r["status"] == "pass"], "motifs[status=pass].motif_id")},
        consumption_note="Size route unavailable when both bars are null; do not coerce null to zero. Check status before use. Named motifs use the measured allowlist.")
    p1.write_json(out / "claim-rules-v1.source.json", rules)
    (out / "si_bar.json.sha256").write_text(digest + "  si_bar.json\n", encoding="utf-8")
    added = ["acceptance-audit.json", "pilot-reproducibility.json", "pilot-curve-differences.csv",
             "claim-rules-v1.source.json", "si_bar.json.sha256"]
    p1.write_json(out / "finalization.json", dict(
        measurement_commit=run["code_commit"], audit_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        audit_code_sha256=p1.sha256(Path(__file__)), run_json_sha256=p1.sha256(out / "run.json"),
        si_bar_json_sha256=digest, artifacts={n: p1.sha256(out / n) for n in added if (out / n).exists()}))
    print(json.dumps({k: audit[k] for k in ("status", "sets", "truth_leaks", "truth_in_set_valid", "correctness_rows_checked", "chance_rows_checked")}))
    print(json.dumps(dict(si_bar_sha256=digest, pilot_reproducibility=pilot_comparison)))


if __name__ == "__main__":
    main()
