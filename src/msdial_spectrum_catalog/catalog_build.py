"""Sequential local builds, measurements and counts-only privacy verification."""
from __future__ import annotations

from contextlib import closing

import json
import platform
import re
import sqlite3
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from .devset import canonical_digest, git_commit, refuse_analysis_output, sha256_file
from .linkage import ingest_unit, read_unit_libraries
from .parsers import read_tsv
from .storage import connect
from .validate import validate_run


def table_counts(connection):
    tables = [row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name")]
    return {table: connection.execute('SELECT COUNT(*) FROM "' + table.replace('"', '""') + '"').fetchone()[0]
            for table in tables}


def privacy_scan(database):
    pattern = re.compile(r"(?:[A-Za-z]:[\\/]|\\\\[^\\/]+[\\/]|//[^/]+/)[^\t\r\n]*\.(msp|lbm2?)", re.I)
    exemptions = {("analysis_run", "output_directory"), ("artifact", "source_path"), ("sample", "raw_file_path")}
    hits = {}
    with closing(connect(database)) as connection:
        for table in table_counts(connection):
            quoted = '"' + table.replace('"', '""') + '"'
            for col in connection.execute(f"PRAGMA table_info({quoted})"):
                name, kind = col[1], col[2]
                if "TEXT" not in kind.upper() or (table, name) in exemptions:
                    continue
                column = '"' + name.replace('"', '""') + '"'
                count = sum(bool(pattern.search(str(row[0]))) for row in connection.execute(
                    f"SELECT {column} FROM {quoted} WHERE {column} IS NOT NULL"))
                if count:
                    hits[table + "." + name] = count
    return {"valid": not hits, "hits": sum(hits.values()), "columns": hits}


def _write_json(path, value):
    path = Path(path)
    refuse_analysis_output(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def ingest_devset(database, devset, *, granularity="representative", members=True, mztab="all"):
    database, devset = Path(database), Path(devset)
    refuse_analysis_output(database)
    selection = json.loads(devset.read_text(encoding="utf-8-sig"))
    if selection.get("schema") != "msdial-annotation-devset.v1":
        raise ValueError("unsupported development-set schema")
    content = {key: selection[key] for key in ("campaigns", "exclusions", "subsets", "units") if key in selection}
    if canonical_digest(content) != selection["content_sha256"]:
        raise ValueError("development-set content digest mismatch")
    # All preflight checks happen before mkdir/initialize/run.json creation.
    for unit in selection["units"]:
        try:
            libraries = read_unit_libraries(unit["workspace"])
        except (OSError, ValueError, KeyError, TypeError):
            if database.name.endswith(".local-vs21.sqlite"):
                # Naming is safe regardless of the unreadable unit's libraries.
                # Its actual failure belongs to the per-unit report below.
                continue
            raise ValueError("cannot verify private-library naming before build") from None
        if any(lib.get("private") is True or str(lib.get("distribution", "")).lower() == "private" for lib in libraries):
            if not database.name.endswith(".local-vs21.sqlite"):
                raise ValueError("private-library catalog requires .local-vs21.sqlite suffix")
    if database.exists() or (database.parent / "run.json").exists():
        raise FileExistsError("build output already exists")
    code_root = Path(__file__).resolve().parents[2]
    result = dict(devset=str(devset), content_sha256=selection["content_sha256"],
                  devset_sha256=sha256_file(devset), git_commit=git_commit(), pin_worktree=str(code_root),
                  python_version=platform.python_version(),
                  parameters=dict(granularity=granularity, members=members, mztab=mztab),
                  units=[], local_only=True,
                  contains="MS-DIAL annotation outputs including private-library hit identities")
    for unit in selection["units"]:
        entry = {"unit_key": unit["unit_key"]}
        try:
            report = ingest_unit(database, unit["workspace"], granularity=granularity,
                                 members=members, mztab=mztab, handoff=unit.get("handoff_path"),
                                 campaign_id=unit["campaign_id"])
            validation = validate_run(database, report.run_id)
            entry.update(report=asdict(report), validation=asdict(validation),
                         valid=report.valid and validation.valid)
        except Exception as error:
            # Exception text can include private row values or input paths. Report type only.
            entry.update(valid=False, error_type=type(error).__name__)
        result["units"].append(entry)
        result["db_bytes"] = database.stat().st_size if database.exists() else 0
        _write_json(database.parent / "run.json", result)
    result["failures"] = sum(not entry["valid"] for entry in result["units"])
    if database.exists():
        with closing(connect(database)) as connection:
            result["table_counts"] = table_counts(connection)
        result["database_sha256"] = sha256_file(database)
    _write_json(database.parent / "run.json", result)
    return result


def size_probe(unit_dirs, out):
    out = Path(out)
    refuse_analysis_output(out)
    if out.exists():
        raise FileExistsError("size report already exists")
    out.parent.mkdir(parents=True, exist_ok=True)
    result = dict(git_commit=git_commit(), pin_worktree=str(Path(__file__).resolve().parents[2]),
                  python_version=platform.python_version(), units=[])
    variants = (("full", True, "all"), ("representative", True, "all"),
                ("representative", False, "all"), ("representative", False, "none"))
    for unit_dir in unit_dirs:
        unit = dict(unit_key=Path(unit_dir).name, variants=[])
        for index, (granularity, members, mztab) in enumerate(variants):
            database = out.parent / f"probe-{Path(unit_dir).name}-{index}.local-vs21.sqlite"
            if database.exists():
                raise FileExistsError("probe scratch database already exists")
            started = time.monotonic()
            try:
                report = ingest_unit(database, unit_dir, granularity=granularity, members=members, mztab=mztab)
                with closing(connect(database)) as connection:
                    connection.execute("VACUUM")
                    counts = table_counts(connection)
                    try:
                        sizes = dict(connection.execute("SELECT name, SUM(pgsize) FROM dbstat GROUP BY name"))
                    except sqlite3.OperationalError:
                        sizes = None
                unit["variants"].append(dict(granularity=granularity, members=members, mztab=mztab,
                    valid=report.valid, error_count=len(report.errors), bytes=database.stat().st_size,
                    wall_seconds=round(time.monotonic() - started, 3), row_counts=counts, dbstat_bytes=sizes))
            finally:
                # Only this probe's exact scratch paths, never directories or input files.
                for suffix in ("", "-wal", "-shm", "-journal"):
                    Path(str(database) + suffix).unlink(missing_ok=True)
        result["units"].append(unit)
        _write_json(out, result)
    return result


def build_report(database, *, devset=None, out=None):
    """Counts-only handback and independent representative-sidecar census."""
    with closing(connect(database)) as connection:
        result = dict(db_bytes=Path(database).stat().st_size, database_sha256=sha256_file(database),
                      table_counts=table_counts(connection))
        result["metadata_matched"] = connection.execute("SELECT COUNT(*) FROM sample WHERE metadata_matched_by IS NOT NULL").fetchone()[0]
        result["samples_with_attributes"] = connection.execute("SELECT COUNT(DISTINCT sample_id) FROM sample_attribute").fetchone()[0]
        result["consensus_nonempty"] = connection.execute("SELECT COUNT(*) FROM spectrum WHERE spectrum_kind='alignment_consensus' AND peak_count>0").fetchone()[0]
        result["sample_spectra"] = connection.execute("SELECT COUNT(*) FROM spectrum WHERE sample_id IS NOT NULL").fetchone()[0]
        result["library_identity_unresolved_rows"] = connection.execute("SELECT COUNT(*) FROM run_annotator WHERE library_sha256 IS NULL").fetchone()[0]
        result["library_identity_unresolved_runs"] = connection.execute("SELECT COUNT(DISTINCT run_id) FROM run_annotator WHERE library_sha256 IS NULL").fetchone()[0]
        result["runs_with_annotators"] = connection.execute("SELECT COUNT(DISTINCT run_id) FROM run_annotator").fetchone()[0]
        result["console_tool_runs"] = connection.execute("SELECT COUNT(*) FROM annotation_tool_run WHERE tool_name='msdial.console'").fetchone()[0]
        result["console_runs_with_libraries"] = connection.execute("SELECT COUNT(DISTINCT t.run_id) FROM annotation_tool_run t JOIN annotation_tool_run_library l USING(tool_run_id) WHERE tool_name='msdial.console'").fetchone()[0]
        result["candidate_evidence_null"] = connection.execute("SELECT COUNT(*) FROM msdial_annotation_candidate WHERE subject_kind='alignment_feature' AND evidence_source IS NULL").fetchone()[0]
        runs = [row[0] for row in connection.execute("SELECT run_id FROM analysis_run ORDER BY run_id")]
        census = Counter({tuple(row[:3]): row[3] for row in connection.execute("""SELECT annotator_id,
            evidence_source, annotation_tag, COUNT(*) FROM msdial_annotation_candidate
            WHERE subject_kind='alignment_feature' AND is_representative=1
            GROUP BY annotator_id, evidence_source, annotation_tag""")})
        result["representative_candidates"] = [dict(annotator=key[0], evidence_source=key[1], tag=key[2], rows=value)
            for key, value in sorted(census.items(), key=lambda item: tuple(str(v) for v in item[0]))]
        result["annotator_tiers"] = [dict(row) for row in connection.execute(
            "SELECT tier, COUNT(*) AS rows FROM run_annotator GROUP BY tier ORDER BY tier")]
    result["validations"] = [dict(run_id=run, valid=validation.valid, counts=validation.counts,
                                  error_count=len(validation.errors), warning_count=len(validation.warnings))
                              for run in runs for validation in [validate_run(database, run)]]
    result["validation_failures"] = sum(not v["valid"] for v in result["validations"])
    result["privacy"] = privacy_scan(database)
    if devset:
        selection = json.loads(Path(devset).read_text(encoding="utf-8-sig"))
        source_census = Counter()
        for unit in selection["units"]:
            for path in sorted(Path(unit["output_dir"]).glob("AlignResult*.mdcandidate.tsv")):
                for _, row in read_tsv(path):
                    if row.get("is_representative", "").strip().lower() == "true":
                        evidence = row.get("evidence_source")
                        source_census[(row.get("annotator_id") or None,
                                       evidence if evidence and evidence.strip() else None,
                                       row.get("annotation_tag_vs1") or None)] += 1
        result["source_census_rows"] = sum(source_census.values())
        result["source_census_matches"] = source_census == census
        result["source_census_difference_rows"] = sum((source_census - census).values()) + sum((census - source_census).values())
    if out:
        _write_json(out, result)
    return result
