"""Read-only, reproducible selections from campaign ledgers (C6)."""
from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from . import __version__

SMOKE = ("e1b72294985e17dcbd3c", "640de497e5f37da80f92")
SMALL = ("0742be0e247ed8f9243d", "99aaf1937a4742e83f18",
         "a70036a8148907c06520", "dbc4a221851c4acbe20c")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_digest(value):
    if isinstance(value, dict) and "campaigns" in value:
        # File mtimes are operational observations, not selection identity:
        # runner heartbeats can update them without changing a transition.
        def without_mtimes(item):
            if isinstance(item, dict):
                return {k: without_mtimes(v) for k, v in item.items()
                        if k not in {"ledger_mtime", "wal_mtime"}}
            return item
        value = dict(value, campaigns=[dict(c, ledger_state=without_mtimes(c.get("ledger_state")))
                                      for c in value["campaigns"]])
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    separators=(",", ":")).encode("utf-8")).hexdigest()


def git_commit():
    return subprocess.check_output(["git", "-C", str(Path(__file__).resolve().parents[2]),
                                    "rev-parse", "HEAD"], text=True).strip()


def refuse_analysis_output(path):
    # Windows paths must also be recognised by POSIX-hosted callers.
    resolved = str(Path(path).resolve()).replace("\\", "/").casefold()
    if "/msdial_public_reanalysis/analysis/" in resolved + "/":
        raise ValueError("refusing output under the analysis root")


def ledger_state(path):
    def stat(p):
        try:
            s = p.stat()
            return s.st_size, s.st_mtime_ns
        except FileNotFoundError:
            return 0, None
    size, mtime = stat(path)
    wal_size, wal_mtime = stat(Path(str(path) + "-wal"))
    return dict(ledger_bytes=size, ledger_mtime=mtime, wal_bytes=wal_size, wal_mtime=wal_mtime)


def read_done_units(path, *, snapshot=False):
    path = Path(path)
    before = ledger_state(path)
    if before["wal_bytes"] and not snapshot:
        raise ValueError("non-empty ledger WAL requires --snapshot")
    uri = path.resolve().as_uri() + "?mode=ro"
    def read(uri):
        connection = sqlite3.connect(uri, uri=True)
        connection.row_factory = sqlite3.Row
        try:
            connection.execute("BEGIN")
            transitions = connection.execute(
                "SELECT unit_key, to_state, seq, at FROM transition ORDER BY seq").fetchall()
            units = connection.execute("SELECT * FROM unit ORDER BY unit_key").fetchall()
            logical_state = dict(max_transition_seq=transitions[-1][2] if transitions else None,
                                 latest_transition_at=transitions[-1][3] if transitions else None)
            return dict((r[0], r[1]) for r in transitions), [dict(r) for r in units], logical_state
        finally:
            connection.close()
    read_mode = "ro-snapshot" if snapshot else "ro"
    try:
        latest, units, logical_state = read(uri)
    except sqlite3.OperationalError:
        if snapshot or ledger_state(path)["wal_bytes"]:
            raise
        latest, units, logical_state = read(uri + "&immutable=1")
        read_mode = "ro-immutable"
    after = ledger_state(path)
    if before != after:
        raise ValueError("ledger or WAL changed during selection; retry")
    selected = []
    known = {u["unit_key"] for u in units}
    if set(latest) - known:
        raise ValueError("transition references an absent unit")
    for unit in units:
        state = latest.get(unit["unit_key"])
        if state is None and unit.get("state") == "done":
            raise ValueError("done unit has no transition")
        if state is not None and unit.get("state") != state:
            raise ValueError("unit.state disagrees with latest transition")
        if state == "done":
            selected.append(unit)
    return selected, read_mode, dict(before, after=after, **logical_state) if snapshot else None


def describe_unit(campaign, unit):
    workspace = Path(unit["workspace"])
    output = workspace / "output"
    files = sorted(output.glob("AlignResult*.mdmsp"))
    if len(files) != 1:
        raise ValueError("expected exactly one alignment mdmsp per unit")
    msp = files[0]
    # Count headers only; annotation identities never enter this report.
    records = nonempty = 0
    with msp.open(encoding="utf-8-sig") as stream:
        for line in stream:
            key, sep, value = line.partition(":")
            if sep and key.strip().casefold() == "num peaks":
                records += 1
                nonempty += int(int(value.strip()) > 0)
    acquisitions = {}
    with (output / "analysis_files.csv").open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            acquisition = row["acquisition_type"]
            acquisitions[acquisition] = acquisitions.get(acquisition, 0) + 1
    result = {key: unit.get(key) for key in (
        "unit_key", "catalog_unit_id", "role", "parent_unit_key", "repository",
        "accession", "instrument", "ion_mode", "workspace", "handoff_path")}
    result.update(campaign_id=campaign, output_dir=str(output),
                  alignment_mdmsp=dict(file_name=msp.name, sha256=sha256_file(msp),
                                       bytes=msp.stat().st_size, records=records,
                                       records_nonempty=nonempty),
                  sample_count=len(list(output.glob("*.mdpeak"))),
                  acquisition_types=dict(sorted(acquisitions.items())),
                  dia=any(a.casefold() != "dda" for a in acquisitions),
                  ion_mode_both=str(unit.get("ion_mode", "")).casefold() == "both")
    return result


def build_devset(campaigns_root, campaigns, out, *, exclusions=None, snapshot=False,
                 name=None, force_new_date=False, command=None):
    if snapshot and not name:
        raise ValueError("--snapshot requires --name")
    name = name or "mbpost-devset"
    out = Path(out)
    refuse_analysis_output(out)
    if out.exists() and not force_new_date:
        raise FileExistsError("dev-set file already exists")
    exclusions = exclusions or {}
    campaign_rows, units, seen = [], [], set()
    for campaign in sorted(set(campaigns)):
        if Path(campaign).name != campaign or "/" in campaign or "\\" in campaign:
            raise ValueError("campaign must be a label")
        done, mode, state = read_done_units(Path(campaigns_root) / campaign / "ledger.sqlite",
                                             snapshot=snapshot)
        campaign_rows.append(dict(campaign_id=campaign,
                                  ledger_label=f"_campaigns/{campaign}/ledger.sqlite",
                                  read_mode=mode, done_units=len(done), ledger_state=state))
        for unit in done:
            key = unit["unit_key"]
            if key in seen:
                raise ValueError("unit occurs in multiple campaigns")
            seen.add(key)
            if key not in exclusions:
                units.append(describe_unit(campaign, unit))
    if set(exclusions) - seen:
        raise ValueError("exclusions do not match selected done units")
    units.sort(key=lambda u: (u["campaign_id"], u["unit_key"]))
    content = dict(campaigns=campaign_rows,
                   exclusions=[dict(unit_key=k, reason=v) for k, v in sorted(exclusions.items())],
                   units=units)
    if name == "mbpost-devset" and not snapshot:
        content["subsets"] = {"smoke": [k for k in SMOKE if k in seen and k not in exclusions],
                              "small": [k for k in SMALL if k in seen and k not in exclusions]}
    result = dict(schema="msdial-annotation-devset.v1", name=name,
                  created_at=datetime.now(timezone.utc).isoformat(),
                  generator=dict(package_version=__version__, git_commit=git_commit(), command=command),
                  **content,
                  totals=dict(units=len(units), samples=sum(u["sample_count"] for u in units),
                              representative_spectra_nonempty=sum(
                                  u["alignment_mdmsp"]["records_nonempty"] for u in units)),
                  content_sha256=canonical_digest(content))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w" if force_new_date else "x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    return result
