"""Plan, run and summarise MS-Emblator over many units, one chunk at a time.

`plan` exports each unit's selections (benchmark, then unresolved) with msemblator_export and splits each
MSP into chunks of at most N records, in alignment-ID order. plan.json records the sha256 of every export
and every chunk, so an importer can prove that a chunk's results came from that chunk
(`msemblator_run.json.input_sha256 == chunk_sha256`) and that the chunk came from the export
(`source_msp_sha256 == <unit>_<selection>_selection.json:msp_sha256`).

`run` runs the driver once per chunk with MS-Emblator's interpreter, sequentially and at below-normal
priority, under a lock that allows one SIRIUS at a time for every runner that shares the lock file. It
resumes: a chunk whose record says exit 0 for the planned input is skipped, and a failed chunk is set
aside and tried again, up to --max-attempts in all. A failure that may come from the login or the web
service is not retried in the same invocation. It uses up an attempt like any other, except the first
three failures of SIRIUS's stored login per chunk. The batch stops after a run of such failures, when
MS-Emblator's work folder is in use, or when <run>/STOP exists. It exits 0 only when no chunk is left failed.

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
import signal
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .msemblator_driver import (EXIT_OK, EXIT_REFUSED, git_commit, git_dirty, parse_log, rankers, scrub, sha256,
                                subscription_limits)
from .msemblator_export import export_unit

PLAN_SCHEMA = "msemblator-batch-plan.v1"
CAMPAIGN_MARKER = "_campaigns"  # the reanalysis (analysis) root holds the campaign ledgers here
SELECTION_ORDER = ("benchmark", "unresolved")
CHUNK_LIMITS = (20, 300)
SIRIUS_TIMEOUT_SECONDS = 1800  # MS-Emblator's single timeout for the whole SIRIUS structure workflow
LOCK_NAME = ".sirius.lock"
_ALIGNMENT_ID = re.compile(r"\|ALIGNMENTID=(\d+)\|")
_FAILED_INDEX = re.compile(r"\.failed-(\d+)$")

# `run` exit codes. 3 is the driver's own "refused" code, reused for a runner that refuses to start.
EXIT_STOPPED = 1  # STOP file, a refusing driver, consecutive service failures, an interrupt, a plan mismatch
EXIT_CHUNKS_FAILED = 5  # the batch went through, but at least one chunk is left failed
# Windows priority class for the driver; its Java and MS-FINDER children inherit it (charter C8).
BELOW_NORMAL_PRIORITY_CLASS = 0x00004000
DRIVER_PRIORITY = "below_normal"


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


REPLACE_RETRIES, REPLACE_PAUSE_SECONDS = 10, 0.2


def write_json(path: Path, data: dict) -> None:
    """Write JSON through a temporary file and a rename, so a kill mid-write leaves the old file whole.

    On Windows the rename fails while another process has the target open (a monitor reading run.json,
    a second runner reading the lock), so it is retried for about two seconds before giving up.
    """
    temporary = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for attempt in range(REPLACE_RETRIES):
        try:
            os.replace(temporary, path)
            return
        except PermissionError:
            if attempt == REPLACE_RETRIES - 1:
                temporary.unlink(missing_ok=True)
                raise
            time.sleep(REPLACE_PAUSE_SECONDS)


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
    write_json(run_dir / "plan.json", plan)
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
    """Take the SIRIUS lock; refuse while a live process holds it.

    The lock covers every runner that uses the same lock file (by default, one per annotation root). A
    holder is stale only when neither the runner nor the driver it last started is alive: a driver
    orphaned by a killed runner may still be running SIRIUS.
    """
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
            pids = [int(current.get(key) or 0) for key in ("pid", "driver_pid")]
            if current.get("host") == socket.gethostname() and not any(pid_alive(pid) for pid in pids):
                path.unlink(missing_ok=True)  # the holder and its driver died without releasing it
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


def note_driver(path: Path, holder: dict, driver_pid: int) -> None:
    """Record in the lock the driver this runner started, so a takeover can tell whether it still runs.
    Best effort: failing to note it must not end the batch while the driver runs, nor free the lock."""
    try:
        write_json(path, {**holder, "driver_pid": driver_pid})
    except OSError:
        pass


# ----------------------------------------------------------------------------------------------- run

def read_record(output: Path) -> dict | None:
    try:
        return json.loads((output / "msemblator_run.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def failed_attempts(output: Path) -> list[Path]:
    attempts = [path for path in output.parent.glob(f"{output.name}.failed-*") if _FAILED_INDEX.search(path.name)]
    return sorted(attempts, key=lambda path: int(_FAILED_INDEX.search(path.name).group(1)))


def set_aside(output: Path) -> Path | None:
    """Move a failed attempt out of the bundle path, keeping it for the record, as the next free
    .failed-<n> (an operator may have removed an earlier one)."""
    if not output.exists() or not any(output.iterdir()):
        if output.exists():
            output.rmdir()
        return None
    used = [int(_FAILED_INDEX.search(path.name).group(1)) for path in failed_attempts(output)]
    target = output.with_name(f"{output.name}.failed-{max(used, default=0) + 1}")
    output.rename(target)
    return target


def is_login_failure(record: dict | None) -> bool:
    """The driver could not use SIRIUS's stored login. Such an attempt says nothing about the chunk, so it
    is not retried at once and (up to FREE_LOGIN_FAILURES per chunk) does not use up an attempt.

    The driver has no exit code of its own for a web-service failure, and its service_error_lines counts
    matching lines anywhere in the log, incidental ones included, so neither exempts an attempt."""
    return bool(record) and record.get("exit_code") != EXIT_OK and record.get("login_ok") is False


def is_service_failure(record: dict | None) -> bool:
    """A failure that may come from the login or the web service rather than the chunk: it feeds only the
    --stop-after rule (a run of them stops the batch), never the attempt count."""
    if not record or record.get("exit_code") == EXIT_OK:
        return False
    return record.get("login_ok") is False or bool(record.get("service_error_lines"))


FREE_LOGIN_FAILURES = 3  # per chunk; later login failures count as attempts, so no chunk reruns for ever


def has_bundle(output: Path) -> bool:
    return output.is_dir() and any(output.iterdir())


def attempts_used(output: Path) -> int:
    """Attempts that count against --max-attempts: every set-aside or current bundle, interrupted ones
    included, except the first FREE_LOGIN_FAILURES that failed on the login."""
    bundles = failed_attempts(output) + ([output] if has_bundle(output) else [])
    counted, free = 0, FREE_LOGIN_FAILURES
    for bundle in bundles:
        if free and is_login_failure(read_record(bundle)):
            free -= 1
        else:
            counted += 1
    return counted


def append_event(run_dir: Path, event: dict) -> None:
    with (run_dir / "batch_events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=False) + "\n")


def kill_tree(process: subprocess.Popen) -> None:
    if sys.platform == "win32":
        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], capture_output=True)
    else:
        process.kill()
    process.wait()


def rescue_work(root: Path, output: Path) -> None:
    """Move MS-Emblator's shared work folder into the chunk's bundle after the driver was killed."""
    work = root / "work"
    if work.exists():
        output.mkdir(parents=True, exist_ok=True)
        shutil.move(str(work), str(output / "work"))


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
    write_json(path, record)


def recorded_plan_sha256(run_dir: Path) -> str | None:
    try:
        return json.loads((run_dir / "run.json").read_text(encoding="utf-8")).get("plan_sha256")
    except (OSError, ValueError):
        return None


def stderr_tail(stderr: bytes, characters: int) -> str:
    """The end of the driver's stderr, with account details removed, for the events and run.json. The
    whole text is scrubbed before it is cut, so a cut cannot split a detail away from its key."""
    return scrub(stderr.decode("utf-8", "replace"))[-characters:]


def start_driver(command: list[str], environment: dict, cwd: Path) -> subprocess.Popen:
    """Start the driver at below-normal priority where the platform has priority classes."""
    options: dict = {}
    if sys.platform == "win32":
        options["creationflags"] = BELOW_NORMAL_PRIORITY_CLASS
    else:
        options["preexec_fn"] = lambda: os.nice(5)
    return subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=environment,
                            cwd=str(cwd), **options)


def run_batch(plan_path: Path, msemblator_root: Path, python: Path | None = None, driver: Path | None = None,
              lock_path: Path | None = None, limit: int | None = None, max_attempts: int = 2,
              stop_after: int = 3, chunk_timeout: float | None = 6 * 3600, hash_libraries: bool = True,
              sirius_directory: Path | None = None, model_assets: str | None = None) -> int:
    run_dir = plan_path.parent
    _refuse_analysis_root(run_dir)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    plan_digest = sha256(plan_path)
    recorded = recorded_plan_sha256(run_dir)
    if recorded is not None and recorded != plan_digest:
        print(f"Refusing to start: plan.json has sha256 {plan_digest}, but run.json was started with {recorded}",
              file=sys.stderr)
        return EXIT_REFUSED
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
        "plan_sha256": plan_digest,
        "model_assets": model_assets,
        "driver_priority": DRIVER_PRIORITY,
        "lock": "<root>/" + LOCK_NAME if lock_path == run_dir.parent / LOCK_NAME else "given with --lock",
        "max_attempts": max_attempts,
    }
    try:
        from .msemblator_driver import affinity_mask
        invocation["affinity_mask"] = affinity_mask()
        if sirius_directory is not None:
            invocation["subscription_limits_before"] = subscription_limits(str(sirius_directory))
    except Exception as error:  # a limits query that fails must not stop the batch
        invocation["subscription_limits_before"] = f"not read: {type(error).__name__}"
    update_run_json(run_dir, invocation)
    # "crashed" until the loop says otherwise: an exception nobody expected must not be recorded as success.
    outcome, ran, consecutive, stopped = "crashed", 0, 0, False
    failed: list[dict] = []  # chunks this invocation ran and left failed
    exhausted: list[str] = []  # chunks skipped because earlier invocations used up their attempts
    try:
        for chunk in plan["chunks"]:
            if (run_dir / "STOP").exists():
                outcome, stopped = "stopped: STOP file", True
                break
            if limit is not None and ran >= limit:
                outcome = f"stopped: --limit {limit}"
                break
            output = run_dir / chunk["output_dir"]
            record = read_record(output)
            if record and record.get("exit_code") == EXIT_OK:
                if record.get("input_sha256") != chunk["chunk_sha256"]:
                    outcome = (f"stopped: {chunk['chunk_id']} has an exit-0 record for another input "
                               f"(input_sha256 {record.get('input_sha256')}, plan {chunk['chunk_sha256']})")
                    stopped = True
                    break
                continue
            attempts = attempts_used(output)
            if attempts >= max_attempts:
                exhausted.append(chunk["chunk_id"])
                continue
            ran += 1
            while attempts < max_attempts:
                set_aside(output)
                output.parent.mkdir(parents=True, exist_ok=True)
                command = [str(python), str(driver), "--msemblator-root", str(root),
                           "--input", str(run_dir / chunk["chunk"]), "--output", str(output), "--mode", "2",
                           "--manifest", str(run_dir / chunk["manifest"])]
                if model_assets:
                    command += ["--model-assets", model_assets]
                started = now()
                environment = dict(os.environ, PYTHONIOENCODING="utf-8")
                process = start_driver(command, environment, root)
                note_driver(lock_path, holder, process.pid)
                reason = None
                try:
                    _, stderr = process.communicate(timeout=chunk_timeout)
                    exit_code = process.returncode
                except subprocess.TimeoutExpired:
                    kill_tree(process)
                    exit_code, stderr, reason = 2, b"", f"timeout after {chunk_timeout} s; process tree killed"
                    rescue_work(root, output)
                except KeyboardInterrupt:
                    # Ctrl+C or Ctrl+Break: stop the chunk, keep what it wrote, and leave MS-Emblator's work
                    # folder free so a restart is not refused. The attempt counts as a failed one.
                    kill_tree(process)
                    rescue_work(root, output)
                    append_event(run_dir, {"chunk_id": chunk["chunk_id"], "attempt": attempts + 1,
                                           "started_at": started, "finished_at": now(), "exit_code": None,
                                           "reason": "interrupted (Ctrl+C or Ctrl+Break); process tree killed",
                                           "failed_steps": None, "stderr_tail": "",
                                           "counts_as_attempt": attempts_used(output) > attempts})
                    outcome = f"interrupted during {chunk['chunk_id']}"
                    raise _Stop
                record = read_record(output)
                if record is None and reason is None and exit_code != EXIT_REFUSED:
                    reason = "no msemblator_run.json"
                service, login = is_service_failure(record), is_login_failure(record)
                # The same rule attempts_used applies on resume: a bundle exists and is not a free login failure.
                used = attempts_used(output)
                counts = used > attempts
                append_event(run_dir, {"chunk_id": chunk["chunk_id"], "attempt": attempts + 1,
                                       "started_at": started, "finished_at": now(), "exit_code": exit_code,
                                       "reason": reason or (record or {}).get("error"),
                                       "failed_steps": (record or {}).get("failed_steps"),
                                       "stderr_tail": stderr_tail(stderr, 2000) if exit_code else "",
                                       "counts_as_attempt": counts})
                if exit_code == EXIT_REFUSED:
                    outcome = f"stopped: driver refused {chunk['chunk_id']}: {stderr_tail(stderr, 500).strip()}"
                    raise _Stop
                # A login or web-service failure (a quota or an outage) would fail the same way at once, and
                # an attempt that left nothing to count could repeat for ever: none is retried in this
                # invocation. Whether the attempt counts is decided above, not here.
                if exit_code == EXIT_OK or login or service or not counts:
                    break
                attempts = used
            if exit_code != EXIT_OK:
                failed.append({"chunk_id": chunk["chunk_id"], "exit_code": exit_code, "service_failure": service})
            if service:
                consecutive += 1
                if consecutive >= stop_after:
                    outcome = f"stopped: {consecutive} consecutive login or web-service failures"
                    stopped = True
                    break
            else:
                consecutive = 0
        else:
            left = len(failed) + len(exhausted)
            outcome = "completed" if not left else f"finished with {left} failed chunks"
    except _Stop:
        stopped = True
    except KeyboardInterrupt:  # between chunks
        outcome, stopped = "interrupted between chunks", True
    except Exception as error:
        outcome = f"crashed: {type(error).__name__}: {scrub(str(error))[:300]}"
        raise
    finally:
        # Record the outcome and free the lock first: whoever sent an interrupt may kill the process soon
        # after, and the limits query below starts SIRIUS, which takes seconds.
        invocation["finished_at"] = now()
        invocation["outcome"] = outcome
        invocation["chunks_attempted"] = ran
        invocation["chunks_failed"] = failed
        invocation["chunks_exhausted"] = exhausted
        try:
            update_run_json(run_dir, invocation)
        finally:
            release_lock(lock_path)  # even when run.json cannot be written
        if sirius_directory is not None and not outcome.startswith(("interrupted", "crashed")):
            try:
                invocation["subscription_limits_after"] = subscription_limits(str(sirius_directory))
            except Exception as error:
                invocation["subscription_limits_after"] = f"not read: {type(error).__name__}"
            update_run_json(run_dir, invocation)
    print(outcome)
    if stopped:
        return EXIT_STOPPED
    return EXIT_CHUNKS_FAILED if failed or exhausted else EXIT_OK


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
                     "partial": record is None and has_bundle(output),
                     "attempts": len(failed_attempts(output)) + (1 if has_bundle(output) else 0),
                     "attempts_counted": attempts_used(output) if state != "done" else None})
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
        "partial_bundles": [r["chunk"]["chunk_id"] for r in rows if r["partial"]],
        "failures": [{"chunk_id": r["chunk"]["chunk_id"], "attempts": r["attempts"],
                      "attempts_counted": r["attempts_counted"],
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


def read_sha256_list(path: Path) -> dict[str, str]:
    """`<sha256> *<relative path>` lines (sha256sum format) as {relative posix path: sha256}."""
    entries = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        digest, _, name = line.strip().partition(" ")
        if digest and name:
            entries[name.lstrip("*").strip().replace("\\", "/")] = digest.lower()
    return entries


def backfill_rankers(plan_path: Path, msemblator_root: Path, expected: Path) -> dict[str, list[str]]:
    """Add `rankers` to chunk records written before the driver recorded it, from each chunk's own log.

    The ranker file is the one the chunk's log names. Its sha256 must equal the entry in `expected`, a
    list written when the models were installed (before the run), or nothing is written. A record whose
    log names no ranking model is left alone: there is nothing to verify. The record says it was
    backfilled, from what, against which list and by which code; nothing else in it changes.
    Returns the chunk IDs backfilled and those whose log named no model.
    """
    run_dir = plan_path.parent
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    listed = read_sha256_list(expected)
    list_sha256 = sha256(expected)
    pending, no_model = [], []
    for chunk in plan["chunks"]:
        output = run_dir / chunk["output_dir"]
        record = read_record(output)
        log = output / "msemblator.log"
        if record is None or "rankers" in record or not log.is_file():
            continue
        loaded = parse_log(log.read_text(encoding="utf-8").splitlines())["ranker_models"]
        found = rankers(loaded, msemblator_root.resolve())
        if not found:
            no_model.append(chunk["chunk_id"])
            continue
        for step, entry in found.items():
            want = listed.get(f"models/{step}/{entry['file']}")
            if want is None or entry["sha256"] != want:
                raise SystemExit(f"Refusing to backfill {chunk['chunk_id']}: models/{step}/{entry['file']} "
                                 f"sha256 {entry['sha256']} does not match the install-time list ({want})")
        pending.append((chunk, output, record, found))
    changed = []
    for chunk, output, record, found in pending:
        record["rankers"] = found
        record["rankers_backfilled"] = {"from": "msemblator.log", "at": now(),
                                        "by_catalog_commit": git_commit(Path(__file__).resolve().parents[2]),
                                        "verified_against": expected.name, "verified_against_sha256": list_sha256,
                                        "match": True}
        write_json(output / "msemblator_run.json", record)
        changed.append(chunk["chunk_id"])
    return {"backfilled": changed, "no_model_in_log": no_model}


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
    run_parser.add_argument("--lock", type=Path,
                            help=f"Default: <plan dir>/../{LOCK_NAME}. Runners exclude each other only when they "
                                 "use the same lock file: give every runner on this machine the same one")
    run_parser.add_argument("--limit", type=int, help="Run at most this many chunks")
    run_parser.add_argument("--max-attempts", type=int, default=2)
    run_parser.add_argument("--stop-after", type=int, default=3)
    run_parser.add_argument("--chunk-timeout", type=float, default=6 * 3600)
    run_parser.add_argument("--no-library-hash", action="store_true")
    run_parser.add_argument("--model-assets", help="key=value,... identifying the installed model/library asset "
                                                   "set (e.g. zenodo_record=...,file=...,md5=...)")

    backfill_parser = commands.add_parser("backfill-rankers", help="Add the rankers field to older records")
    backfill_parser.add_argument("--plan", required=True, type=Path)
    backfill_parser.add_argument("--msemblator-root", required=True, type=Path)
    backfill_parser.add_argument("--expected-sha256", required=True, type=Path,
                                 help="sha256sum-format list of the installed models, written before the run")

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
        if hasattr(signal, "SIGBREAK"):  # Windows: Ctrl+Break would otherwise end Python without cleanup
            signal.signal(signal.SIGBREAK, signal.default_int_handler)
        return run_batch(args.plan, root, args.python, lock_path=args.lock, limit=args.limit,
                         max_attempts=args.max_attempts, stop_after=args.stop_after,
                         chunk_timeout=args.chunk_timeout, hash_libraries=not args.no_library_hash,
                         sirius_directory=root / "external" / "sirius", model_assets=args.model_assets)
    if args.command == "backfill-rankers":
        print(json.dumps(backfill_rankers(args.plan, args.msemblator_root, args.expected_sha256), indent=2))
        return 0
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
