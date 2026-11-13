"""Run MS-Emblator on one exported MSP with the SIRIUS login SIRIUS itself stores.

Run this file with MS-Emblator's own Python (its virtual environment), not the catalog's: it imports
MS-Emblator and nothing from this package.

MS-Emblator's command line requires --sirius_pass for modes 2 and 3 and logs in with it before each
SIRIUS step. A password on a command line is visible in the process list, and an agent must never hold
one. This driver runs MS-Emblator's own main() unchanged, except that its sirius_login is replaced by a
check that `sirius.exe login --show` reports a stored login: a person logs in once with
`sirius.exe login -u <account> -p`, and every run uses that login.

MS-Emblator keeps its intermediates in one shared work/ folder inside its checkout. The driver refuses to
start while that folder holds a file, and moves the folder into the run's output when the run ends, so
two runs never share intermediates and every run keeps its own.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

PLACEHOLDER = "<stored-sirius-login>"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stored_login(sirius_directory: str) -> None:
    """Raise unless the SIRIUS in sirius_directory has a stored login."""
    sirius = Path(sirius_directory) / "sirius.exe"
    result = subprocess.run([str(sirius), "login", "--show"], cwd=sirius.parent, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=600)
    text = result.stdout + result.stderr
    if "Not Logged in" in text or '"sid"' not in text:
        raise RuntimeError(
            f"SIRIUS at {sirius.parent} has no stored login; log in once with "
            "'sirius.exe login -u <account> -p' and run again.")


def git_commit(root: Path) -> str | None:
    try:
        return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--msemblator-root", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", type=int, choices=(1, 2, 3), default=2)
    args = parser.parse_args()
    root = args.msemblator_root.resolve()
    work = root / "work"
    if work.exists() and any(path.is_file() for path in work.rglob("*")):
        print(f"Refusing to start: {work} holds files from another run.", file=sys.stderr)
        return 3
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        print(f"Refusing to start: {output} is not empty.", file=sys.stderr)
        return 3
    output.mkdir(parents=True, exist_ok=True)

    jdk = root / "external" / "jdk21" / "bin"
    if jdk.is_dir():
        os.environ["PATH"] = str(jdk) + os.pathsep + os.environ.get("PATH", "")
    sys.path.insert(0, str(root / "src"))

    from msemblator.cli import msemblator as cli
    from msemblator.runners import sirius_cmd, sirius_struc_cmd
    from msemblator.workflows import struc_main

    def checked_login(sirius_directory, username, password):
        stored_login(sirius_directory)

    for module in (cli, sirius_cmd, sirius_struc_cmd, struc_main):
        module.sirius_login = checked_login

    record = {
        "tool": "MS-Emblator",
        "msemblator_commit": git_commit(root),
        "mode": args.mode,
        "input": args.input.name,
        "input_sha256": sha256(args.input),
        "parameter_file_sha256": sha256(root / "config" / "msemblator_parameter_file.yaml"),
        "sirius_login": "stored",
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    argv = ["msemblator", "--input", str(args.input), "--output", str(output), "--mode", str(args.mode)]
    if args.mode in (2, 3):
        argv += ["--sirius_user", PLACEHOLDER, "--sirius_pass", PLACEHOLDER]
    sys.argv = argv
    try:
        exit_code = cli.main() or 0
    except Exception as error:  # the record must say how the run ended, whatever ended it
        record["error"] = f"{type(error).__name__}: {error}"
        exit_code = 2
    finally:
        record["finished_at"] = datetime.now(timezone.utc).isoformat()
        if work.exists():
            shutil.move(str(work), str(output / "work"))
    record["exit_code"] = exit_code
    record["outputs"] = {path.name: sha256(path) for path in sorted(output.glob("*.csv"))}
    (output / "msemblator_run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
