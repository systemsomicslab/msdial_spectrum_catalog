"""Run MS-Emblator on one exported MSP with the SIRIUS login SIRIUS itself stores.

Run this file with MS-Emblator's own Python (its virtual environment), not the catalog's: it imports
MS-Emblator and nothing from this package, and it uses only the standard library itself.

MS-Emblator's command line requires --sirius_pass for modes 2 and 3 and logs in with it before each
SIRIUS step. A password on a command line is visible in the process list, and an agent must never hold
one. This driver runs MS-Emblator's own main() unchanged, except that its sirius_login is replaced by a
check that `sirius.exe login --show` reports a stored login: a person logs in once with
`sirius.exe login -u <account> -p`, and every run uses that login.

MS-Emblator keeps its intermediates in one shared work/ folder inside its checkout. The driver refuses to
start while that folder holds a file, and moves the folder into the run's output when the run ends, so
two runs never share intermediates and every run keeps its own.

Besides the four CSVs, a run keeps what an importer and a benchmark need to trust them:

- msemblator.log, every line MS-Emblator and its tools printed, time-stamped and scrubbed of account
  details; the stage timings and the failed steps are parsed from it;
- metfrag/, the MetFrag outputs, which MS-Emblator leaves in its checkout and clears on the next run;
- name_map.tsv, MS-Emblator's integer compound ID mapped back to the exported NAME and, through the
  export manifest, to unit and alignment ID;
- msemblator_run.json, with the SIRIUS version, the failed steps, the environment and a hash of every
  harvested file, including the CANOPUS summaries MS-Emblator writes but never reads.

Exit codes: 0 clean; 1 MS-Emblator reported failed steps; 2 exception (including a missing login or a
name map that is not 1:1); 3 refused to start, nothing written; 4 the CSVs exist but no CANOPUS summary
was found although the SIRIUS structure step did not fail.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path

PLACEHOLDER = "<stored-sirius-login>"
CSV_OUTPUTS = ("formula_summary.csv", "formula_score.csv", "structure_score.csv", "structure_summary.csv")
METFRAG_KEPT_IN_CHECKOUT = frozenset({"MetFragCommandLine-2.5.0.jar", "library_psv_v2.txt",
                                      "metfrag_StructureDB.txt"})
LIBRARY_DATABASES = ("MsfinderStructureDB_all.txt", "metfrag_StructureDB.txt", "sirius_structure_db.siriusdb")
MANIFEST_COLUMNS = ("name", "unit_id", "alignment_id")  # never reference_*: those are the benchmark truth

EXIT_OK, EXIT_FAILED_STEPS, EXIT_EXCEPTION, EXIT_REFUSED, EXIT_HARVEST = 0, 1, 2, 3, 4

_EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
_SID = re.compile(r'("?(?:sid|refresh_?token|access_?token|token)"?\s*[:=]\s*)("[^"]*"|\S+)', re.IGNORECASE)
_SERVICE_ERROR = re.compile(
    r"HTTP\D{0,12}\b(?:401|403)\b|\b(?:401|403) (?:unauthori[sz]ed|forbidden)|status\D{0,8}\b(?:401|403)\b"
    r"|unauthori[sz]ed|forbidden|terms of (?:service|use)|not accepted"
    r"|web-?service.{0,80}not returning|no stored login|login process|login failed",
    re.IGNORECASE)
_FAILED_STEPS = re.compile(r"(formula|structure) elucidation completed with failed steps:\s*(.+)$",
                           re.IGNORECASE)
_STAGE_MARK = re.compile(r"^(SIRIUS|MS-FINDER|MetFrag|msbuddy) processing (start|complete)\s*$")
_STAGE_TIME = re.compile(r"(SIRIUS|MetFrag|MSFinder|Summary generation) (?:processing )?time: ([0-9.]+) seconds")
_FORMULA_TOTAL = re.compile(r"Processing completed in ([0-9.]+) seconds")
_LOG_STAMP = re.compile(r"^(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d+\+00:00) (.*)$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scrub(text: str) -> str:
    """Remove account details (e-mail addresses, session IDs, tokens) from text that will be stored."""
    return _SID.sub(r"\1<redacted>", _EMAIL.sub("<account>", text))


def stored_login(sirius_directory: str) -> None:
    """Raise unless the SIRIUS in sirius_directory has a stored login. Prints nothing."""
    sirius = Path(sirius_directory) / "sirius.exe"
    result = subprocess.run([str(sirius), "login", "--show"], cwd=sirius.parent, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=600)
    text = result.stdout + result.stderr
    if "Not Logged in" in text or '"sid"' not in text:
        raise RuntimeError(
            f"SIRIUS at {sirius.parent} has no stored login; log in once with "
            "'sirius.exe login -u <account> -p' and run again.")


SUBSCRIPTION_LIMIT_KEYS = ("countQueries", "compoundLimit", "maxQueriesPerCompound", "compoundHashRecordingTime",
                           "expirationDate")


def subscription_limits(sirius_directory: str) -> dict:
    """The usage limits `login --show` reports, and nothing else of its output.

    Only the allowlisted keys are returned: the rest of the output names the account.
    """
    sirius = Path(sirius_directory) / "sirius.exe"
    result = subprocess.run([str(sirius), "login", "--show"], cwd=sirius.parent, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=600)
    text = result.stdout + result.stderr
    decoder, index = json.JSONDecoder(), 0
    while (start := text.find("{", index)) >= 0:
        try:
            value, index = decoder.raw_decode(text, start)
        except ValueError:
            index = start + 1
            continue
        if isinstance(value, dict) and any(key in value for key in SUBSCRIPTION_LIMIT_KEYS):
            return {key: value.get(key) for key in SUBSCRIPTION_LIMIT_KEYS
                    if isinstance(value.get(key), (bool, int, float, type(None)))}
    return {}


def git_commit(root: Path) -> str | None:
    try:
        return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def git_dirty(root: Path) -> bool | None:
    """True when the checkout differs from its commit, ignoring the untracked library/ the assets add."""
    try:
        lines = subprocess.run(["git", "-C", str(root), "status", "--porcelain"], capture_output=True,
                               text=True, check=True).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError):
        return None
    return any(line[3:].strip().strip('"').rstrip("/") != "library" for line in lines if line.strip())


def sirius_version(root: Path) -> str | None:
    """The SIRIUS version, from the name of the REST-service jar the bundled SIRIUS runs."""
    for jar in sorted((root / "external" / "sirius").rglob("sirius_rest_service-*-boot.jar")):
        match = re.match(r"sirius_rest_service-(.+)-boot\.jar$", jar.name)
        if match:
            return match.group(1)
    return None


def affinity_mask() -> str | None:
    """The process CPU affinity mask, as hex, where the platform tells it."""
    if hasattr(os, "sched_getaffinity"):
        return hex(sum(1 << cpu for cpu in os.sched_getaffinity(0)))
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        process, system = ctypes.c_size_t(), ctypes.c_size_t()
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        kernel32.GetProcessAffinityMask.argtypes = (wintypes.HANDLE, ctypes.POINTER(ctypes.c_size_t),
                                                    ctypes.POINTER(ctypes.c_size_t))
        if kernel32.GetProcessAffinityMask(kernel32.GetCurrentProcess(), ctypes.byref(process),
                                           ctypes.byref(system)):
            return hex(process.value)
    return None


def environment(root: Path) -> dict:
    """What the run depended on, apart from the input: code, models, databases, interpreter, CPUs."""
    models = {}
    names = {"formula": (), "structure": (), "ranker": ()}
    try:
        sys.path.insert(0, str(root / "src"))
        from msemblator import paths as msemblator_paths  # noqa: PLC0415 (MS-Emblator's own module)
        names = {"formula": msemblator_paths.FORMULA_MODEL_FILES,
                 "structure": msemblator_paths.STRUCTURE_MODEL_FILES,
                 "ranker": msemblator_paths.RANKER_MODEL_FILES}
    except ImportError:
        pass
    for kind, directory in (("formula", "formula"), ("structure", "structure")):
        for name in tuple(names[kind]) + tuple(names["ranker"]):
            path = root / "models" / directory / name
            if path.is_file():
                models[f"{directory}/{name}"] = sha256(path)
    databases = {}
    for name in LIBRARY_DATABASES:
        path = root / "library" / name
        if path.is_file():
            stat = path.stat()
            databases[name] = {"bytes": stat.st_size,
                               "mtime": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()}
    return {
        "msemblator_commit": git_commit(root),
        "msemblator_dirty": git_dirty(root),
        "model_sha256": models,
        "library_databases": databases,
        "python": platform.python_version(),
        "affinity_mask": affinity_mask(),
        "logical_cpus": os.cpu_count(),
    }


class Tee:
    """A text stream that passes writes on and appends each complete line, time-stamped and scrubbed, to a log."""

    def __init__(self, stream, log, lock: threading.Lock):
        self.stream, self.log, self.lock, self.pending = stream, log, lock, ""

    def write(self, text: str) -> int:
        if self.stream is not None:
            try:
                self.stream.write(text)
            except (OSError, ValueError, UnicodeError):
                pass
        with self.lock:
            self.pending += text
            *lines, self.pending = self.pending.split("\n")
            for line in lines:
                self._append(line)
        return len(text)

    def _append(self, line: str) -> None:
        stamp = datetime.now(timezone.utc).isoformat(timespec="microseconds")
        self.log.write(f"{stamp} {scrub(line.rstrip(chr(13)))}\n")
        self.log.flush()

    def flush(self) -> None:
        if self.stream is not None:
            try:
                self.stream.flush()
            except (OSError, ValueError):
                pass

    def close_pending(self) -> None:
        with self.lock:
            if self.pending:
                self._append(self.pending)
                self.pending = ""

    def isatty(self) -> bool:
        return False

    @property
    def encoding(self) -> str:
        return "utf-8"


def parse_log(lines: list[str]) -> dict:
    """Stage timings, failed steps and service-error lines from a driver log (time-stamped lines)."""
    section = None
    opened: dict[str, datetime] = {}
    stages: dict[str, float] = {}
    measured: dict[str, float] = {}  # the times MS-Emblator logs itself win over the start/complete marks
    failed: list[str] = []
    service_errors = 0
    formula_tool_end = None
    for raw in lines:
        match = _LOG_STAMP.match(raw.rstrip("\n"))
        if not match:
            continue
        stamp, text = datetime.fromisoformat(match.group(1)), match.group(2).strip()
        if text.startswith("Running formula elucidation") and not text.startswith("Running formula elucidation only"):
            section = "formula"
            continue
        if text.startswith("Running structure elucidation"):
            section = "structure"
            continue
        if _SERVICE_ERROR.search(text):
            service_errors += 1
        mark = _STAGE_MARK.match(text)
        if mark and section:
            tool = mark.group(1).lower().replace("-", "")
            key = f"{section}.{tool}"
            if mark.group(2) == "start":
                opened[key] = stamp
            elif key in opened:
                stages[key] = round((stamp - opened.pop(key)).total_seconds(), 3)
                if section == "formula":
                    formula_tool_end = stamp
            continue
        timed = _STAGE_TIME.search(text)
        if timed and section == "structure":
            tool = {"SIRIUS": "sirius", "MetFrag": "metfrag", "MSFinder": "msfinder",
                    "Summary generation": "summary"}[timed.group(1)]
            measured[f"structure.{tool}"] = float(timed.group(2))
            continue
        total = _FORMULA_TOTAL.search(text)
        if total and section == "formula":
            stages["formula.total"] = float(total.group(1))
            if formula_tool_end is not None:
                stages["formula.summary"] = round((stamp - formula_tool_end).total_seconds(), 3)
            continue
        failed_match = _FAILED_STEPS.search(text)
        if failed_match:
            step = failed_match.group(1).lower()
            failed += [f"{step}:{tool.strip()}" for tool in failed_match.group(2).split(",") if tool.strip()]
    stages.update(measured)
    return {"stage_seconds": stages, "failed_steps": failed, "service_error_lines": service_errors}


def read_msp_names(path: Path) -> list[tuple[str, dict[str, str]]]:
    """NAME and the COMMENT tokens of every record, in file order."""
    records: list[tuple[str, dict[str, str]]] = []
    name, tokens = None, {}
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            key, _, value = line.partition(":")
            key = key.strip().lower()
            if key == "name":
                if name is not None:
                    records.append((name, tokens))
                name, tokens = value.strip(), {}
            elif key == "comment" and name is not None:
                for token in value.strip().strip("|").split("|"):
                    token_key, sep, token_value = token.partition("=")
                    if sep:
                        tokens[token_key.strip()] = token_value.strip()
    if name is not None:
        records.append((name, tokens))
    return records


def build_name_map(input_msp: Path, work: Path, manifest: Path | None) -> tuple[list[dict], list[str]]:
    """Rows mapping MS-Emblator's integer ID to the exported NAME and its unit and alignment ID.

    Returns the rows and the problems that make the map unsound; an empty list means it is 1:1 and
    covers every input record.
    """
    problems: list[str] = []
    id_change = work / "save_folder" / "formula_fixed_msp" / "id_change.msp"
    formula_fixed = work / "save_folder" / "formula_fixed_msp" / "formula_fixed.msp"
    if not id_change.is_file():
        return [], [f"missing {id_change.name}"]
    input_names = [name for name, _ in read_msp_names(input_msp)]
    changed = read_msp_names(id_change)
    sent_to_structure = {name for name, _ in read_msp_names(formula_fixed)} if formula_fixed.is_file() else None
    keys: dict[str, dict[str, str]] = {}
    if manifest is not None:
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("name") in keys:
                    problems.append(f"manifest repeats NAME {row['name']}")
                keys[row.get("name", "")] = {column: row.get(column, "") for column in MANIFEST_COLUMNS}
    rows = []
    seen_ids, seen_names = set(), set()
    for msemblator_id, tokens in changed:
        original = tokens.get("ORIGNAME", "")
        if msemblator_id in seen_ids:
            problems.append(f"MS-Emblator ID {msemblator_id} repeats")
        if original in seen_names:
            problems.append(f"NAME {original} maps to two MS-Emblator IDs")
        seen_ids.add(msemblator_id)
        seen_names.add(original)
        key = keys.get(original)
        if manifest is not None and key is None:
            problems.append(f"NAME {original} is not in the manifest")
        unit_id = key["unit_id"] if key else tokens.get("UNIT", "")
        alignment_id = key["alignment_id"] if key else tokens.get("ALIGNMENTID", "")
        if tokens.get("UNIT") and unit_id != tokens["UNIT"] or \
                tokens.get("ALIGNMENTID") and str(alignment_id) != tokens["ALIGNMENTID"]:
            problems.append(f"NAME {original}: COMMENT keys differ from the manifest")
        if sent_to_structure is None:
            step = "not_run"
        else:
            step = "sent" if msemblator_id in sent_to_structure else "skipped_no_formula"
        rows.append({"msemblator_id": msemblator_id, "name": original, "unit_id": unit_id,
                     "alignment_id": alignment_id, "structure_step": step})
    if sorted(seen_names) != sorted(input_names):
        missing = sorted(set(input_names) - seen_names)
        extra = sorted(seen_names - set(input_names))
        problems.append(f"name map does not cover the input: {len(missing)} missing, {len(extra)} extra")
    if len(input_names) != len(set(input_names)):
        problems.append("input NAMEs repeat")
    return rows, problems


def write_name_map(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["msemblator_id", "name", "unit_id", "alignment_id",
                                                    "structure_step"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def copy_metfrag(root: Path, output: Path, since: float) -> list[str]:
    """Copy this run's MetFrag files out of the checkout, which the next run clears."""
    source = root / "external" / "metfrag"
    if not source.is_dir():
        return []
    copied = []
    for path in sorted(source.rglob("*")):
        if not path.is_file() or path.name in METFRAG_KEPT_IN_CHECKOUT or path.stat().st_mtime < since:
            continue
        destination = output / "metfrag" / path.relative_to(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        copied.append(destination.relative_to(output).as_posix())
    return copied


def harvest(output: Path) -> tuple[dict, list[str]]:
    """{relative path: {sha256, bytes}} for every product worth keeping, and the CANOPUS summary names."""
    work = output / "work"
    wanted: list[Path] = [output / name for name in CSV_OUTPUTS]
    for project in (work / "formula" / "sirius" / "output" / "sirius_project.sirius",
                    work / "structure" / "sirius_output" / "sirius_project.sirius"):
        if project.is_file():
            wanted.append(project)
        elif project.is_dir():
            wanted += sorted(p for p in project.rglob("*") if p.is_file())
    for folder in (work / "structure" / "sirius_output", work / "formula" / "sirius" / "output"):
        wanted += sorted(folder.rglob("*.tsv")) if folder.is_dir() else []
    for folder in (work / "formula" / "msfinder_output", work / "structure" / "msfinder_output",
                   output / "metfrag"):
        wanted += sorted(p for p in folder.rglob("*") if p.is_file()) if folder.is_dir() else []
    files = {}
    for path in wanted:
        if path.is_file():
            relative = path.relative_to(output).as_posix()
            files[relative] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    structure_out = work / "structure" / "sirius_output"
    canopus = sorted(p.relative_to(structure_out).as_posix() for p in structure_out.rglob("canopus*summary*.tsv")) \
        if structure_out.is_dir() else []
    return files, canopus


def run(args: argparse.Namespace) -> int:
    root = args.msemblator_root.resolve()
    work = root / "work"
    if work.exists() and any(path.is_file() for path in work.rglob("*")):
        print(f"Refusing to start: {work} holds files from another run.", file=sys.stderr)
        return EXIT_REFUSED
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        print(f"Refusing to start: {output} is not empty.", file=sys.stderr)
        return EXIT_REFUSED
    output.mkdir(parents=True, exist_ok=True)

    log_path = (args.log or output / "msemblator.log").resolve()
    log_handle = log_path.open("a", encoding="utf-8", newline="\n")
    lock = threading.Lock()
    stdout, stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = Tee(stdout, log_handle, lock), Tee(stderr, log_handle, lock)
    started = datetime.now(timezone.utc)
    record: dict = {
        "tool": "MS-Emblator",
        "msemblator_commit": git_commit(root),
        "sirius_version": sirius_version(root),
        "mode": args.mode,
        "input": args.input.name,
        "input_sha256": sha256(args.input),
        "input_records": len(read_msp_names(args.input)),
        "manifest": args.manifest.name if args.manifest else None,
        "parameter_file_sha256": sha256(root / "config" / "msemblator_parameter_file.yaml"),
        "sirius_login": "stored",
        "started_at": started.isoformat(),
    }
    exit_code = EXIT_EXCEPTION
    login_failures: list[str] = []
    try:
        try:
            record["environment"] = environment(root)
        except Exception as error:  # a missing fact about the environment must not cost the run
            record["environment"] = {"error": scrub(f"{type(error).__name__}: {error}")}
        jdk = root / "external" / "jdk21" / "bin"
        if jdk.is_dir():
            os.environ["PATH"] = str(jdk) + os.pathsep + os.environ.get("PATH", "")
        sys.path.insert(0, str(root / "src"))

        from msemblator.cli import msemblator as cli  # noqa: PLC0415 (imported after the tee is in place)
        from msemblator.runners import sirius_cmd, sirius_struc_cmd  # noqa: PLC0415
        from msemblator.workflows import struc_main  # noqa: PLC0415

        def checked_login(sirius_directory, username, password):
            try:
                stored_login(sirius_directory)
            except Exception as error:
                login_failures.append(scrub(str(error)))
                raise

        for module in (cli, sirius_cmd, sirius_struc_cmd, struc_main):
            module.sirius_login = checked_login

        sys.argv = ["msemblator", "--input", str(args.input), "--output", str(output), "--mode", str(args.mode)]
        if args.mode in (2, 3):
            sys.argv += ["--sirius_user", PLACEHOLDER, "--sirius_pass", PLACEHOLDER]
        exit_code = cli.main() or 0
    except Exception as error:  # the record must say how the run ended, whatever ended it
        record["error"] = scrub(f"{type(error).__name__}: {error}")
        exit_code = EXIT_EXCEPTION
    finally:
        record["finished_at"] = datetime.now(timezone.utc).isoformat()
        for stream in (sys.stdout, sys.stderr):
            if isinstance(stream, Tee):
                stream.close_pending()
        sys.stdout, sys.stderr = stdout, stderr
        log_handle.close()
        if work.exists():
            shutil.move(str(work), str(output / "work"))

    if login_failures:
        record["error"] = f"login: {login_failures[0]}"
        record["login_ok"] = False
        exit_code = EXIT_EXCEPTION
    else:
        record["login_ok"] = True
    parsed = parse_log(log_path.read_text(encoding="utf-8").splitlines())
    record["stage_seconds"] = parsed["stage_seconds"]
    record["failed_steps"] = parsed["failed_steps"]
    record["service_error_lines"] = parsed["service_error_lines"]
    record["metfrag_files"] = copy_metfrag(root, output, started.timestamp() - 1)

    rows, problems = build_name_map(args.input, output / "work", args.manifest)
    write_name_map(output / "name_map.tsv", rows)
    record["name_map_ok"] = not problems
    record["name_map"] = {
        "records": len(rows),
        "sent_to_structure": sum(row["structure_step"] == "sent" for row in rows),
        "skipped_no_formula": sum(row["structure_step"] == "skipped_no_formula" for row in rows),
        "problems": problems[:20],
    }
    files, canopus = harvest(output)
    record["harvest"] = files
    record["canopus_summary_files"] = canopus
    record["outputs"] = {name: sha256(output / name) for name in CSV_OUTPUTS if (output / name).is_file()}
    all_csvs = all(name in record["outputs"] for name in CSV_OUTPUTS)
    if exit_code in (EXIT_OK, EXIT_FAILED_STEPS) and args.mode in (2, 3) and all_csvs and not canopus \
            and "structure:SIRIUS" not in record["failed_steps"]:
        record["error"] = "harvest: no CANOPUS summary in work/structure/sirius_output"
        exit_code = EXIT_HARVEST
    if problems and exit_code in (EXIT_OK, EXIT_FAILED_STEPS, EXIT_HARVEST):
        record["error"] = f"name_map: {problems[0]}"
        exit_code = EXIT_EXCEPTION
    record["exit_code"] = exit_code
    record["duration_seconds"] = round((datetime.fromisoformat(record["finished_at"]) - started).total_seconds(), 3)
    text = scrub(json.dumps(record, indent=2, ensure_ascii=False))
    (output / "msemblator_run.json").write_text(text + "\n", encoding="utf-8")
    return exit_code


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--msemblator-root", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", type=int, choices=(1, 2, 3), default=2)
    parser.add_argument("--manifest", type=Path,
                        help="The export manifest; only its name, unit_id and alignment_id columns are read")
    parser.add_argument("--log", type=Path, help="Default: <output>/msemblator.log")
    return run(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
