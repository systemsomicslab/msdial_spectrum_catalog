"""Allowlisted campaign metadata linkage; no manifest document is persisted."""
from __future__ import annotations

from contextlib import closing

import csv
import json
import ntpath
import re
from pathlib import Path

from .annotation import LibraryIdentity, ToolRunInput, record_tool_run
from .devset import refuse_analysis_output, sha256_file
from .ingest import ingest_run
from .parsers import read_analysis_csv
from .storage import connect, transaction

PROJECT_FIELDS = ("repository", "accession", "analysis_unit_id", "source_subrecord_id", "title",
                  "description", "public_url", "license", "separation", "acquisition_mode", "ion_mode")
OUTPUT_FIELDS = ("created_at", "analysis_type", "msdial_console_version", "msdial_interactive_version",
                 "annotation_pipeline_profile", "software_provenance_status")
CAMPAIGN_FIELDS = ("campaign_id", "catalog_unit_id", "unit_key", "parent_unit_key", "role", "state", "manifest_digest")
LIBRARY_FIELDS = ("name", "role", "sha256", "bytes", "library_kind", "distribution", "private", "version", "license")


def safe_value(value):
    if isinstance(value, str) and (re.search(r"^[A-Za-z]:[\\/]|^\\\\", value)
                                  or (re.search(r"\.msp|\.lbm", value, re.I) and re.search(r"[\\/]", value))):
        raise ValueError("path-like value in allowlisted metadata field")
    if isinstance(value, dict):
        return {safe_value(k): safe_value(v) for k, v in value.items()}
    if isinstance(value, list):
        return [safe_value(v) for v in value]
    return value


def _take(data, fields):
    return {field: safe_value(data[field]) for field in fields if field in data}


def _load(path):
    with Path(path).open(encoding="utf-8-sig") as stream:
        return json.load(stream)


def read_unit_metadata(unit_dir, *, handoff=None, campaign_id=None):
    unit_dir = Path(unit_dir)
    provenance = _load(unit_dir / "provenance/run-manifest.json")
    if provenance.get("schema") != "msdial-public-reanalysis-run.v1":
        raise ValueError("unsupported provenance manifest schema")
    project = _take(provenance["project"], PROJECT_FIELDS)
    campaign_raw = _load(unit_dir / "campaign-record.json")
    campaign = _take(campaign_raw, CAMPAIGN_FIELDS)
    if project["analysis_unit_id"] != campaign["unit_key"] or campaign["unit_key"] != unit_dir.name:
        raise ValueError("unit identifiers disagree")
    if campaign_id and campaign.get("campaign_id") != campaign_id:
        raise ValueError("campaign identifiers disagree")
    raw = _load(unit_dir / "output/run-manifest.json")
    output = _take(raw, OUTPUT_FIELDS)
    output["target_omics"] = _take(raw.get("target_omics", {}), ("repository_seed", "answered", "pinned", "applied"))
    output["console"] = _take(raw.get("console", {}), ("binary_sha256", "source_kind"))
    output["libraries"] = [_take(library, LIBRARY_FIELDS) for library in raw.get("libraries", [])]
    campaign["libraries_as_run"] = [_take(library, LIBRARY_FIELDS) for library in campaign_raw.get("libraries_as_run", [])]
    handoff_data = None
    if handoff:
        handoff_raw = _load(handoff)
        if handoff_raw.get("schema") != "msdial-repository-reanalysis-handoff.v1":
            raise ValueError("unsupported handoff schema")
        handoff_data = _take(handoff_raw, ("unit_attributes", "technical_settings"))
    return project, campaign, output, handoff_data


def read_annotator_settings(path):
    """The sole settings reader. Project named columns before materialising any row.

    In particular the path column is never accessed, retained, hashed or logged.
    """
    path = Path(path)
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, [])
        columns = {name: header.index(name) for name in ("annotator_id", "evidence_tier", "library_kind")
                   if name in header}
        if "annotator_id" not in columns:
            raise ValueError("settings lack annotator_id")
        result = []
        for cells in reader:
            row = {name: safe_value(cells[index] if index < len(cells) else "")
                   for name, index in columns.items()}
            evidence, kind = row.get("evidence_tier", ""), row.get("library_kind", "")
            tier = ("in_silico" if kind.lower() == "predicted" else "strict" if "high" in evidence.lower()
                    else "relaxed" if "low" in evidence.lower() else "unknown")
            row["tier"] = tier
            row["library_kind"] = kind or ("acquired" if tier in {"strict", "relaxed"} else "")
            result.append(row)
        return result


def match_annotator_library(annotator_id, libraries):
    """Resolve named annotators by role + stem; generic MSP names by role + kind.

    No guess is made when more than one manifest library matches.
    """
    generic_kind = {"msp_high_quality": "acquired", "msp_low_quality": "acquired",
                    "msp_in_silico": "predicted"}.get(annotator_id.casefold())
    if generic_kind:
        candidates = [lib for lib in libraries if lib.get("role", "").casefold() == "msp"
                      and lib.get("library_kind", "").casefold() == generic_kind]
        return candidates[0] if len(candidates) == 1 else None
    prefix, separator, name = annotator_id.partition(":")
    role = {"msp": "msp", "lbmdb": "lbm"}.get(prefix.casefold())
    if not separator or role is None:
        return None
    candidates = []
    for library in libraries:
        if str(library.get("role", "")).casefold() not in ({"lbm", "lbm2"} if role == "lbm" else {"msp"}):
            continue
        filename = library.get("name", "")
        stem = Path(filename).stem
        if name.casefold() in {filename.casefold(), stem.casefold(),
                              stem.casefold() + "_high", stem.casefold() + "_low",
                              stem.casefold() + "_high_quality", stem.casefold() + "_low_quality"}:
            candidates.append(library)
    return candidates[0] if len(candidates) == 1 else None


def _identity(row):
    return LibraryIdentity(row["name"], row["sha256"], role=row.get("role") or "",
                           library_kind=row.get("library_kind"), distribution=row.get("distribution"),
                           library_version=row.get("version"), license=row.get("license"))


def _attribute_value(value):
    safe_value(value)
    return value if isinstance(value, str) else json.dumps(value, ensure_ascii=False, sort_keys=True)


def ingest_unit(database, unit_dir, *, granularity="full", members=True, mztab="all",
                handoff=None, campaign_id=None, expected=None):
    refuse_analysis_output(database)
    unit_dir = Path(unit_dir)
    project, campaign, output, handoff_data = read_unit_metadata(unit_dir, handoff=handoff, campaign_id=campaign_id)
    if expected:
        for key, value in expected.items():
            actual = output.get(key) if key in OUTPUT_FIELDS else project.get(key)
            if value is not None and actual != value:
                raise ValueError("manual argument disagrees with unit manifest")
    # Preflight all identity and sample metadata before creating the database.
    manifests = {lib["sha256"]: lib for lib in output["libraries"]}
    as_run = campaign["libraries_as_run"]
    libraries = [dict(manifests.get(lib["sha256"], {}), **lib) for lib in as_run] if as_run else output["libraries"]
    identities = tuple(_identity(library) for library in libraries)
    settings_path = unit_dir / "output/msp_annotator_settings.tsv"
    annotators = read_annotator_settings(settings_path)
    sample_paths = sorted((unit_dir / "output").glob("*_sample_metadata_reviewed.tsv"))
    if len(sample_paths) > 1:
        raise ValueError("multiple reviewed sample metadata files")
    sample_rows = []
    if sample_paths:
        with sample_paths[0].open(encoding="utf-8-sig", newline="") as stream:
            sample_rows = list(csv.DictReader(stream, delimiter="\t"))
        # raw_file is used only for matching, never copied as manifest metadata.
        for row in sample_rows:
            for key, value in row.items():
                safe_value(key)
                if key != "raw_file":
                    safe_value(value)
    report = ingest_run(database, unit_dir / "output", project["repository"], project["accession"],
                        project["analysis_unit_id"], study_title=project.get("title"),
                        separation_type=project.get("separation"), ion_mode=project.get("ion_mode"),
                        acquisition_type=project.get("acquisition_mode"),
                        msdial_version=output.get("msdial_console_version"),
                        interactive_version=output.get("msdial_interactive_version"),
                        granularity=granularity, members=members, mztab=mztab)
    if "already full" in report.warnings:
        return report
    if as_run and set(manifests) != {lib["sha256"] for lib in as_run}:
        report.warnings.append("campaign and output library sha256 sets differ")
    if not settings_path.exists():
        report.warnings.append("annotator settings missing")
    if not sample_paths:
        report.warnings.append("reviewed sample metadata missing")
    analysis = read_analysis_csv(unit_dir / "output/analysis_files.csv")
    with transaction(database) as connection:
        run = connection.execute("SELECT * FROM analysis_run WHERE run_id=?", (report.run_id,)).fetchone()
        unit_id = run["analysis_unit_id"]
        connection.execute("UPDATE analysis_run SET campaign_id=?, manifest_digest=? WHERE run_id=?",
                           (campaign.get("campaign_id"), campaign.get("manifest_digest"), report.run_id))
        connection.execute("""UPDATE analysis_unit SET catalog_unit_id=?, parent_unit_key=?, unit_role=?,
            source_subrecord_id=?, target_omics_applied=? WHERE analysis_unit_id=?""",
            (campaign.get("catalog_unit_id"), campaign.get("parent_unit_key"), campaign.get("role"),
             project.get("source_subrecord_id"), output["target_omics"].get("applied"), unit_id))
        connection.execute("""UPDATE study SET description=?, public_url=?, license=? WHERE study_id=
            (SELECT study_id FROM analysis_unit WHERE analysis_unit_id=?)""",
            (project.get("description"), project.get("public_url"), project.get("license"), unit_id))
        sources = [("provenance_run_manifest", unit_dir / "provenance/run-manifest.json", "provenance/run-manifest.json"),
                   ("output_run_manifest", unit_dir / "output/run-manifest.json", "output/run-manifest.json"),
                   ("campaign_record", unit_dir / "campaign-record.json", "campaign-record.json")]
        if sample_paths:
            sources.append(("sample_metadata_reviewed", sample_paths[0], "output/" + sample_paths[0].name))
        if handoff:
            sources.append(("campaign_handoff", Path(handoff), "handoff.json"))
        for role, path, label in sources:
            connection.execute("INSERT OR REPLACE INTO run_metadata_source VALUES (?, ?, ?, ?, ?)",
                               (report.run_id, role, label, sha256_file(path), path.stat().st_size))
        groups = {"run_manifest.project": project, "run_manifest.target_omics": output["target_omics"]}
        if handoff_data:
            groups.update({"handoff." + key: value for key, value in handoff_data.items()})
        for source, values in groups.items():
            for key, value in values.items():
                connection.execute("INSERT OR REPLACE INTO analysis_unit_attribute VALUES (?, ?, ?, ?)",
                                   (unit_id, key, _attribute_value(value), source))
        matched = set()
        for sample in connection.execute("SELECT sample_id, sample_name FROM sample WHERE analysis_unit_id=?", (unit_id,)).fetchall():
            stem = sample["sample_name"]
            matches = [(i, row) for i, row in enumerate(sample_rows) if row.get("sample_id") == stem or row.get("raw_file") == stem]
            rule = "stem"
            if not matches:
                basename = ntpath.basename(analysis.get(stem, {}).get("file_path", ""))
                matches = [(i, row) for i, row in enumerate(sample_rows) if basename and row.get("raw_file") == basename]
                rule = "analysis_csv_basename"
            if len(matches) != 1:
                report.warnings.append("sample metadata unmatched or ambiguous")
                continue
            index, row = matches[0]
            matched.add(index)
            report.metadata_matched += 1
            connection.execute("""UPDATE sample SET repository_sample_id=?, class_id=?, source_name=?,
                raw_file_paired_by=?, metadata_matched_by=? WHERE sample_id=?""",
                (row.get("sample_id"), row.get("class_id"), row.get("source_name"),
                 row.get("raw_file_paired_by"), rule, sample["sample_id"]))
            for key, value in row.items():
                if " / " in key and value:
                    connection.execute("INSERT OR REPLACE INTO sample_attribute VALUES (?, ?, ?, 'sample_metadata_reviewed')",
                                       (sample["sample_id"], key, value))
        if len(matched) < len(sample_rows):
            report.warnings.append(f"{len(sample_rows) - len(matched)} metadata rows unmatched")
        for row in annotators:
            library = match_annotator_library(row["annotator_id"], output["libraries"])
            report.library_identity_unresolved += int(library is None)
            connection.execute("INSERT OR REPLACE INTO run_annotator VALUES (?, ?, ?, ?, ?, ?, ?)",
                               (report.run_id, row["annotator_id"], row.get("evidence_tier"), row.get("library_kind"),
                                row["tier"], library["name"] if library else None, library["sha256"] if library else None))
    record_tool_run(database, ToolRunInput("msdial.console", output.get("msdial_console_version"),
                    parameters={"analysis_type": output.get("analysis_type"),
                                "annotation_pipeline_profile": output.get("annotation_pipeline_profile"),
                                "target_omics_applied": output["target_omics"].get("applied"),
                                "msdial_interactive_version": output.get("msdial_interactive_version")},
                    provenance={"console_binary_sha256": output["console"].get("binary_sha256"),
                                "console_source_kind": output["console"].get("source_kind"),
                                "libraries_source": "libraries_as_run" if as_run else "output_manifest"},
                    input_fingerprint=run["run_fingerprint"], libraries=identities), run_id=report.run_id)
    report.db_bytes_after = Path(database).stat().st_size
    return report


def export_sample_attributes(database, out, *, unit=None):
    refuse_analysis_output(out)
    with closing(connect(database)) as connection:
        samples = connection.execute("""SELECT s.*, u.external_unit_id FROM sample s JOIN analysis_unit u
            USING(analysis_unit_id) WHERE (? IS NULL OR u.external_unit_id=? OR u.analysis_unit_id=?)
            ORDER BY u.external_unit_id, s.sample_name""", (unit, unit, unit)).fetchall()
        values = {row["sample_id"]: {} for row in samples}
        for row in connection.execute("SELECT * FROM sample_attribute ORDER BY source"):
            if row["sample_id"] in values:
                values[row["sample_id"]][row["attribute_key"]] = row["attribute_value"]
    keys = sorted({key for attrs in values.values() for key in attrs})
    with Path(out).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["unit", "sample", "class_id"] + keys)
        for sample in samples:
            writer.writerow([sample["external_unit_id"], sample["sample_name"], sample["class_id"]] +
                            [values[sample["sample_id"]].get(key, "") for key in keys])
    return len(samples)
