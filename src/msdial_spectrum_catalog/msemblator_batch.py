"""Plan, run and summarise MS-Emblator over many units, one chunk at a time.

`plan` exports each unit's selections (benchmark, then unresolved) with msemblator_export and splits each
MSP into chunks of at most N records, in alignment-ID order. plan.json records the sha256 of every export
and every chunk, so an importer can prove that a chunk's results came from that chunk
(`msemblator_run.json.input_sha256 == chunk_sha256`) and that the chunk came from the export
(`source_msp_sha256 == <unit>_<selection>_selection.json:msp_sha256`).

`run` runs the driver once per chunk with MS-Emblator's interpreter, sequentially, under a lock that
allows one SIRIUS at a time on this machine. It resumes: a chunk whose record says exit 0 is skipped,
and a failed chunk is set aside and tried again, up to --max-attempts in all. It stops after a run of
login or web-service failures, when MS-Emblator's work folder is in use, or when <run>/STOP exists.

`status` and `throughput` read plan.json and the per-chunk records; they write nothing but
throughput.json.

The standard library only: this module is run with the catalog's interpreter from a pinned checkout.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from .msemblator_driver import EXIT_OK, EXIT_REFUSED, git_commit, git_dirty, sha256, subscription_limits
from .msemblator_export import export_unit

PLAN_SCHEMA = "msemblator-batch-plan.v1"
CAMPAIGN_MARKER = "_campaigns"  # the reanalysis (analysis) root holds the campaign ledgers here
SELECTION_ORDER = ("benchmark", "unresolved")
CHUNK_LIMITS = (20, 300)
SIRIUS_TIMEOUT_SECONDS = 1800  # MS-Emblator's single timeout for the whole SIRIUS structure workflow
LOCK_NAME = ".sirius.lock"
_ALIGNMENT_ID = re.compile(r"\|ALIGNMENTID=(\d+)\|")


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _refuse_analysis_root(path: Path) -> None:
    """Unit outputs are read-only inputs: refuse a run folder inside the reanalysis root, which is the folder
    that holds the campaign ledgers (_campaigns)."""
    resolved = path.resolve()
    for folder in (resolved, *resolved.parents):
        if (folder / CAMPAIGN_MARKER).is_dir():
            raise SystemExit(f"Refusing to write under the analysis root {folder}: {path}")


# ---------------------------------------------------------------------------------------------- plan

def load_units(devset: Path | None, pairs: list[str]) -> tuple[list[dict], dict]:
    """Units to plan, from a dev-set file (C6) or from <unit_key>=<output dir> pairs."""
    units, sources = [], {}
    if devset is not None:
        data = json.loads(devset.read_text(encoding="utf-8"))
        sources["devset"] = {"file": devset.name, "sha256": sha256(devset),
                             "content_sha256": data.get("content_sha256")}
        for unit in data["units"]:
            output = unit.get("output_dir") or str(Path(unit["workspace"]) / "output")
            units.append({"unit_key": unit["unit_key"], "output_dir": output,
                          **{key: unit.get(key) for key in ("campaign_id", "accession", "ion_mode", "instrument",
                                                            "dia", "ion_mode_both")}})
    for pair in pairs:
        key, sep, output = pair.partition("=")
        if not sep or not key or not output:
            raise SystemExit(f"--units takes <unit_key>=<output dir>, not {pair!r}")
        units.append({"unit_key": key, "output_dir": output})
    keys = [unit["unit_key"] for unit in units]
    if len(set(keys)) != len(keys):
        raise SystemExit("A unit_key appears twice")
    return units, sources


def msp_blocks(path: Path) -> list[str]:
    """The records of an exported MSP, each as its text with one trailing blank line, in alignment-ID order."""
    text = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    blocks = [block.strip("\n") + "\n\n" for block in re.split(r"\n\s*\n", text) if block.strip()]

    def alignment_id(block: str) -> int:
        match = _ALIGNMENT_ID.search(block)
        if not match:
            raise ValueError(f"{path.name}: a record has no ALIGNMENTID")
        return int(match.group(1))

    return sorted(blocks, key=alignment_id)


def chunk_size_from(seconds_per_structure_spectrum: float) -> int:
    """N = floor(0.5 x 1800 / t_struct), capped to [20, 300]: half of SIRIUS's whole-workflow timeout."""
    if seconds_per_structure_spectrum <= 0:
        return CHUNK_LIMITS[1]
    size = math.floor(0.5 * SIRIUS_TIMEOUT_SECONDS / seconds_per_structure_spectrum)
    return max(CHUNK_LIMITS[0], min(CHUNK_LIMITS[1], size))


def make_plan(units: list[dict], run_dir: Path, run_id: str, chunk_size: int, chunk_size_origin: str,
              selections: tuple[str, ...] = SELECTION_ORDER, max_spectra: int | None = None,
              sources: dict | None = None) -> dict:
    _refuse_analysis_root(run_dir)
    if (run_dir / "plan.json").exists():
        raise SystemExit(f"{run_dir / 'plan.json'} exists; a plan is written once")
    if chunk_size < 1:
        raise SystemExit("--chunk-size must be positive")
    exports, chunks = [], []
    for unit in units:
        key = unit["unit_key"]
        destination = run_dir / "input" / key
        for selection in selections:
            summary = export_unit(Path(unit["output_dir"]), key, destination, selection, max_spectra)
            source = destination / summary["msp"]
            blocks = msp_blocks(source)
            entry = {
                "unit_key": key,
                "selection": selection,
                "msp": f"input/{key}/{summary['msp']}",
                "source_msp_sha256": summary["msp_sha256"],
                "selection_json": f"input/{key}/{key}_{selection}_selection.json",
                "manifest": f"input/{key}/{summary['manifest']}",
                "records": len(blocks),
                "representative_spectra_by_tier": summary["representative_spectra_by_tier"],
                "counts": summary["counts"],
                "chunks": [],
            }
            for index in range(0, len(blocks), chunk_size):
                number = index // chunk_size
                name = f"{key}_{selection}_c{number:03d}.msp"
                path = destination / name
                with path.open("w", encoding="utf-8", newline="\n") as handle:
                    handle.write("".join(blocks[index:index + chunk_size]))
                chunk = {
                    "chunk_id": f"{key}/{selection}/c{number:03d}",
                    "unit_key": key,
                    "selection": selection,
                    "chunk": f"input/{key}/{name}",
                    "chunk_sha256": sha256(path),
                    "records": len(blocks[index:index + chunk_size]),
                    "manifest": entry["manifest"],
                    "output_dir": f"units/{key}/{selection}/c{number:03d}",
                }
                entry["chunks"].append(chunk["chunk_id"])
                chunks.append(chunk)
            exports.append(entry)
    # Benchmark first, so evaluation can start early; within a selection the largest unit runs last.
    size = {(e["unit_key"], e["selection"]): e["records"] for e in exports}
    rank = {selection: position for position, selection in enumerate(selections)}
    chunks.sort(key=lambda c: (rank[c["selection"]], size[(c["unit_key"], c["selection"])], c["unit_key"],
                               c["chunk_id"]))
    for order, chunk in enumerate(chunks):
        chunk["order"] = order
    totals = {selection: sum(e["records"] for e in exports if e["selection"] == selection)
              for selection in selections}
    by_tier: dict[str, int] = {}
    for entry in exports:
        if entry["selection"] == selections[0]:
            for tier, count in entry["representative_spectra_by_tier"].items():
                by_tier[tier] = by_tier.get(tier, 0) + count
    plan = {
        "schema": PLAN_SCHEMA,
        "run_id": run_id,
        "created_at": now(),
        "catalog_commit": git_commit(Path(__file__).resolve().parent),
        "chunk_size": chunk_size,
        "chunk_size_origin": chunk_size_origin,
        "selections": list(selections),
        "max_spectra": max_spectra,
        "sources": sources or {},
        "units": units,
        "exports": exports,
        "chunks": chunks,
        "totals": {"units": len(units), "chunks": len(chunks), "records_by_selection": totals,
                   "representative_spectra_by_tier": by_tier},
    }
    (run_dir / "plan.json").write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return plan


def export_entry(plan: dict, unit_key: str, selection: str) -> dict:
    """The export a chunk came from. The one place that knows the plan's key names."""
    for entry in plan["exports"]:
        if entry["unit_key"] == unit_key and entry["selection"] == selection:
            return entry
    raise KeyError(f"no export for {unit_key} {selection}")


# ---------------------------------------------------------------------------------------------- lock

def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.OpenProcess.restype = wintypes.HANDLE
        handle = kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not handle:
            return ctypes.get_last_error() == 5  # access denied: it exists
        try:
            code = wintypes.DWORD()
            if not kernel32.GetExitCodeProcess(handle, ctypes.byref(code)):
                return True
            return code.value == 259  # STILL_ACTIVE
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def acquire_lock(path: Path, run_id: str) -> tuple[bool, dict | None]:
    """Take the machine-wide SIRIUS lock; refuse while a live process holds it."""
    holder = {"pid": os.getpid(), "host": socket.gethostname(), "run_id": run_id, "started_at": now()}
    path.parent.mkdir(parents=True, exist_ok=True)
    for _ in range(2):
        try:
            descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            try:
                current = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                current = {}
            if current.get("host") == socket.gethostname() and not pid_alive(int(current.get("pid", 0))):
                path.unlink(missing_ok=True)  # the holder died without releasing it
                continue
            return False, current
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(holder) + "\n")
        return True, holder
    return False, None


def release_lock(path: Path) -> None:
    try:
        current = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return
    if current.get("pid") == os.getpid() and current.get("host") == socket.gethostname():
        path.unlink(missing_ok=True)


# ----------------------------------------------------------------------------------------------- run

def read_record(output: Path) -> dict | None:
    try:
        return json.loads((output / "msemblator_run.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def failed_attempts(output: Path) -> list[Path]:
    return sorted(output.parent.glob(f"{output.name}.failed-*"))


def set_aside(output: Path) -> Path | None:
    """Move a failed attempt out of the bundle path, keeping it for the record."""
    if not output.exists() or not any(output.iterdir()):
        if output.exists():
            output.rmdir()
        return None
    target = output.with_name(f"{output.name}.failed-{len(failed_attempts(output)) + 1}")
    output.rename(target)
    return target


def is_service_failure(record: dict | None) -> bool:
    if not record or record.get("exit_code") == EXIT_OK:
        return False
    return record.get("login_ok") is False or bool(record.get("service_error_lines"))


def append_event(run_dir: Path, event: dict) -> None:
    with (run_dir / "batch_events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=False) + "\n")


def kill_tree(process: subprocess.Popen) -> None:
    if sys.platform == "win32":
        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], capture_output=True)
    else:
        process.kill()
    process.wait()


def library_hashes(root: Path) -> dict:
    hashes = {}
    for name in ("MsfinderStructureDB_all.txt", "metfrag_StructureDB.txt", "sirius_structure_db.siriusdb"):
        path = root / "library" / name
        if path.is_file():
            hashes[name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    return hashes


def update_run_json(run_dir: Path, invocation: dict) -> None:
    path = run_dir / "run.json"
    record = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {
        "run_id": json.loads((run_dir / "plan.json").read_text(encoding="utf-8"))["run_id"],
        "plan_sha256": sha256(run_dir / "plan.json"),
        "invocations": []}
    invocations = [i for i in record["invocations"] if i.get("invocation_id") != invocation["invocation_id"]]
    record["invocations"] = invocations + [invocation]
    path.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def run_batch(plan_path: Path, msemblator_root: Path, python: Path | None = None, driver: Path | None = None,
              lock_path: Path | None = None, limit: int | None = None, max_attempts: int = 2,
              stop_after: int = 3, chunk_timeout: float | None = 6 * 3600, hash_libraries: bool = True,
              sirius_directory: Path | None = None) -> int:
    run_dir = plan_path.parent
    _refuse_analysis_root(run_dir)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    root = msemblator_root.resolve()
    python = python or root / ".venv" / "Scripts" / "python.exe"
    driver = (driver or Path(__file__).with_name("msemblator_driver.py")).resolve()
    lock_path = lock_path or run_dir.parent / LOCK_NAME
    acquired, holder = acquire_lock(lock_path, plan["run_id"])
    if not acquired:
        print(f"Refusing to start: {lock_path} is held by {holder}", file=sys.stderr)
        return EXIT_REFUSED
    catalog_root = Path(__file__).resolve().parents[2]
    invocation = {
        "invocation_id": f"{os.getpid()}-{now()}",
        "started_at": now(),
        "catalog_checkout": catalog_root.name,
        "catalog_commit": git_commit(catalog_root),
        "catalog_dirty": git_dirty(catalog_root),
        "msemblator_commit": git_commit(root),
        "msemblator_dirty": git_dirty(root),
        "interpreters": {"batch": platform.python_version(), "driver": python.name},
        "driver_sha256": sha256(driver),
        "parameter_file_sha256": sha256(root / "config" / "msemblator_parameter_file.yaml"),
        "library_sha256": library_hashes(root) if hash_libraries else "not computed",
        "host_logical_cpus": os.cpu_count(),
        "plan_sha256": sha256(plan_path),
    }
    try:
        from .msemblator_driver import affinity_mask
        invocation["affinity_mask"] = affinity_mask()
        if sirius_directory is not None:
            invocation["subscription_limits_before"] = subscription_limits(str(sirius_directory))
    except Exception as error:  # a limits query that fails must not stop the batch
        invocation["subscription_limits_before"] = f"not read: {type(error).__name__}"
    update_run_json(run_dir, invocation)
    outcome, ran, consecutive = "completed", 0, 0
    try:
        for chunk in plan["chunks"]:
            if (run_dir / "STOP").exists():
                outcome = "stopped: STOP file"
                break
            if limit is not None and ran >= limit:
                outcome = f"stopped: --limit {limit}"
                break
            output = run_dir / chunk["output_dir"]
            record = read_record(output)
            if record and record.get("exit_code") == EXIT_OK:
                continue
            attempts = len(failed_attempts(output)) + (1 if output.exists() and any(output.iterdir()) else 0)
            if attempts >= max_attempts:
                continue
            ran += 1
            while attempts < max_attempts:
                set_aside(output)
                output.parent.mkdir(parents=True, exist_ok=True)
                command = [str(python), str(driver), "--msemblator-root", str(root),
                           "--input", str(run_dir / chunk["chunk"]), "--output", str(output), "--mode", "2",
                           "--manifest", str(run_dir / chunk["manifest"])]
                started = now()
                environment = dict(os.environ, PYTHONIOENCODING="utf-8")
                process = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                           env=environment, cwd=str(root))
                reason = None
                try:
                    _, stderr = process.communicate(timeout=chunk_timeout)
                    exit_code = process.returncode
                except subprocess.TimeoutExpired:
                    kill_tree(process)
                    exit_code, stderr, reason = 2, b"", f"timeout after {chunk_timeout} s; process tree killed"
                    work = root / "work"
                    if work.exists():
                        output.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(work), str(output / "work"))
                attempts += 1
                record = read_record(output)
                if record is None and reason is None and exit_code != EXIT_REFUSED:
                    reason = "no msemblator_run.json"
                event = {"chunk_id": chunk["chunk_id"], "attempt": attempts, "started_at": started,
                         "finished_at": now(), "exit_code": exit_code,
                         "reason": reason or (record or {}).get("error"),
                         "failed_steps": (record or {}).get("failed_steps"),
                         "stderr_tail": stderr.decode("utf-8", "replace")[-2000:] if exit_code else ""}
                append_event(run_dir, event)
                if exit_code == EXIT_REFUSED:
                    return_reason = stderr.decode("utf-8", "replace").strip()[-500:]
                    outcome = f"stopped: driver refused {chunk['chunk_id']}: {return_reason}"
                    raise _Stop
                if exit_code == EXIT_OK:
                    break
            if is_service_failure(record):
                consecutive += 1
                if consecutive >= stop_after:
                    outcome = f"stopped: {consecutive} consecutive login or web-service failures"
                    break
            else:
                consecutive = 0
    except _Stop:
        pass
    finally:
        invocation["finished_at"] = now()
        invocation["outcome"] = outcome
        invocation["chunks_attempted"] = ran
        if sirius_directory is not None:
            try:
                invocation["subscription_limits_after"] = subscription_limits(str(sirius_directory))
            except Exception as error:
                invocation["subscription_limits_after"] = f"not read: {type(error).__name__}"
        update_run_json(run_dir, invocation)
        release_lock(lock_path)
    print(outcome)
    return EXIT_OK if outcome.startswith(("completed", "stopped: --limit")) else 1


class _Stop(Exception):
    pass


# -------------------------------------------------------------------------------- status, throughput

def directory_bytes(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file()) if path.is_dir() else 0


def collect(plan: dict, run_dir: Path) -> list[dict]:
    rows = []
    for chunk in plan["chunks"]:
        output = run_dir / chunk["output_dir"]
        record = read_record(output)
        if record is None:
            state = "failed" if failed_attempts(output) else "pending"
        else:
            state = "done" if record.get("exit_code") == EXIT_OK else "failed"
        rows.append({"chunk": chunk, "output": output, "record": record, "state": state,
                     "attempts": len(failed_attempts(output)) + (1 if record else 0)})
    return rows


def status(plan_path: Path) -> dict:
    run_dir = plan_path.parent
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    rows = collect(plan, run_dir)
    done = [r for r in rows if r["state"] == "done"]
    spectra_done = sum(r["chunk"]["records"] for r in done)
    seconds = sum(r["record"].get("duration_seconds", 0) for r in done)
    size = sum(directory_bytes(r["output"]) for r in done)
    intervals = sorted((r["record"]["started_at"], r["record"]["finished_at"], r["chunk"]["chunk_id"])
                       for r in rows if r["record"] and r["record"].get("finished_at"))
    overlaps = [(a[2], b[2]) for a, b in zip(intervals, intervals[1:]) if b[0] < a[1]]
    return {
        "run_id": plan["run_id"],
        "chunks": {"total": len(rows), "done": len(done), "failed": sum(r["state"] == "failed" for r in rows),
                   "pending": sum(r["state"] == "pending" for r in rows)},
        "spectra": {"planned": sum(r["chunk"]["records"] for r in rows), "done": spectra_done},
        "seconds_per_spectrum": round(seconds / spectra_done, 3) if spectra_done else None,
        "bytes_per_spectrum": round(size / spectra_done) if spectra_done else None,
        "failures": [{"chunk_id": r["chunk"]["chunk_id"], "attempts": r["attempts"],
                      "exit_code": (r["record"] or {}).get("exit_code"),
                      "reason": (r["record"] or {}).get("error"),
                      "failed_steps": (r["record"] or {}).get("failed_steps")}
                     for r in rows if r["state"] == "failed"],
        "overlapping_chunks": overlaps,
    }


def throughput(plan_path: Path, targets: dict[str, int], first: int | None = None) -> dict:
    """Per-stage seconds and bytes per spectrum from done chunks, and the extrapolation to the targets."""
    run_dir = plan_path.parent
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    rows = [r for r in collect(plan, run_dir) if r["state"] == "done"]
    rows.sort(key=lambda r: r["chunk"]["order"])
    if first is not None:
        rows = rows[:first]
    spectra = sum(r["chunk"]["records"] for r in rows)
    sent = sum((r["record"].get("name_map") or {}).get("sent_to_structure", 0) for r in rows)
    stages: dict[str, float] = {}
    for row in rows:
        for stage, seconds in (row["record"].get("stage_seconds") or {}).items():
            stages[stage] = stages.get(stage, 0.0) + seconds
    wall = sum(r["record"].get("duration_seconds", 0) for r in rows)
    work_bytes = sum(directory_bytes(r["output"] / "work") for r in rows)
    project_bytes = sum(v["bytes"] for r in rows for k, v in (r["record"].get("harvest") or {}).items()
                        if ".sirius" in k)
    total_bytes = sum(directory_bytes(r["output"]) for r in rows)
    per = (lambda value: round(value / spectra, 4) if spectra else None)
    t_struct = stages.get("structure.sirius", 0.0) / sent if sent else None
    result = {
        "run_id": plan["run_id"],
        "chunks": [r["chunk"]["chunk_id"] for r in rows],
        "spectra": spectra,
        "sent_to_structure": sent,
        "stage_seconds_total": {k: round(v, 3) for k, v in sorted(stages.items())},
        "seconds_per_spectrum": {k: per(v) for k, v in sorted(stages.items())},
        "wall_seconds_per_spectrum": per(wall),
        "structure_sirius_seconds_per_sent_spectrum": round(t_struct, 4) if t_struct else None,
        "chunk_size": chunk_size_from(t_struct) if t_struct else None,
        "bytes_per_spectrum": {"work": per(work_bytes), "sirius_projects": per(project_bytes),
                               "bundle": per(total_bytes)},
        "service_error_lines": sum(r["record"].get("service_error_lines", 0) for r in rows),
        "failed_steps": sorted({s for r in rows for s in r["record"].get("failed_steps") or []}),
        "extrapolation": {},
    }
    if spectra:
        for name, count in targets.items():
            result["extrapolation"][name] = {
                "spectra": count,
                "wall_days": round(wall / spectra * count / 86400, 2),
                "bundle_bytes": round(total_bytes / spectra * count),
            }
    (run_dir / "throughput.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


# ----------------------------------------------------------------------------------------------- CLI

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m msdial_spectrum_catalog.msemblator_batch",
                                     description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    plan_parser = commands.add_parser("plan", help="Export and chunk the selections; write plan.json")
    plan_parser.add_argument("--devset", type=Path)
    plan_parser.add_argument("--units", action="append", default=[], metavar="UNIT_KEY=OUTPUT_DIR")
    plan_parser.add_argument("--run-id", required=True)
    plan_parser.add_argument("--root", type=Path, required=True,
                             help="The msemblator annotation root; the plan goes to <root>/<run-id>")
    plan_parser.add_argument("--chunk-size", type=int, required=True)
    plan_parser.add_argument("--chunk-size-origin", default="given")
    plan_parser.add_argument("--selections", default=",".join(SELECTION_ORDER))
    plan_parser.add_argument("--max-spectra", type=int)

    run_parser = commands.add_parser("run", help="Run the driver on every pending chunk, one at a time")
    run_parser.add_argument("--plan", required=True, type=Path)
    run_parser.add_argument("--msemblator-root", required=True, type=Path)
    run_parser.add_argument("--python", type=Path, help="Default: <msemblator-root>/.venv/Scripts/python.exe")
    run_parser.add_argument("--lock", type=Path, help=f"Default: <plan dir>/../{LOCK_NAME}")
    run_parser.add_argument("--limit", type=int, help="Run at most this many chunks")
    run_parser.add_argument("--max-attempts", type=int, default=2)
    run_parser.add_argument("--stop-after", type=int, default=3)
    run_parser.add_argument("--chunk-timeout", type=float, default=6 * 3600)
    run_parser.add_argument("--no-library-hash", action="store_true")

    status_parser = commands.add_parser("status", help="Summarise plan.json against the chunk records")
    status_parser.add_argument("--plan", required=True, type=Path)

    throughput_parser = commands.add_parser("throughput", help="Seconds and bytes per spectrum; extrapolation")
    throughput_parser.add_argument("--plan", required=True, type=Path)
    throughput_parser.add_argument("--target", action="append", default=[], metavar="NAME=SPECTRA")
    throughput_parser.add_argument("--first", type=int)

    args = parser.parse_args(argv)
    if args.command == "plan":
        selections = tuple(s.strip() for s in args.selections.split(",") if s.strip())
        units, sources = load_units(args.devset, args.units)
        run_dir = args.root / args.run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        plan = make_plan(units, run_dir, args.run_id, args.chunk_size, args.chunk_size_origin, selections,
                         args.max_spectra, sources)
        print(json.dumps(plan["totals"], indent=2))
        return 0
    if args.command == "run":
        root = args.msemblator_root
        return run_batch(args.plan, root, args.python, lock_path=args.lock, limit=args.limit,
                         max_attempts=args.max_attempts, stop_after=args.stop_after,
                         chunk_timeout=args.chunk_timeout, hash_libraries=not args.no_library_hash,
                         sirius_directory=root / "external" / "sirius")
    if args.command == "status":
        print(json.dumps(status(args.plan), indent=2))
        return 0
    targets = {}
    for pair in args.target:
        name, _, count = pair.partition("=")
        targets[name] = int(count)
    print(json.dumps(throughput(args.plan, targets, args.first), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
