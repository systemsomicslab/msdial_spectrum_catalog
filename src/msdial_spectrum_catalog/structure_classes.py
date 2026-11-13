"""Load MS-FINDER's three-column direct-parent lookup without retaining its path."""
from __future__ import annotations

from contextlib import closing

import csv
import json
import re
from pathlib import Path

from .annotation import LibraryIdentity, ToolRunInput, record_tool_run
from .classes import LOOKUP_LEVEL
from .devset import refuse_analysis_output, sha256_file
from .identifiers import make_id
from .storage import initialize, transaction


def load_structure_classes(database, icd, *, update_run_json=False):
    """Load lookup rows atomically; updating adjacent run.json requires opt-in."""
    refuse_analysis_output(database)
    path = Path(icd)
    digest = sha256_file(path)
    library = LibraryIdentity(path.name, digest, role="class_lookup")
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        expected = {"InChIKey", "Direct parent name", "Direct parent ID"}
        if not expected.issubset(reader.fieldnames or []):
            raise ValueError("unsupported ICD columns")
        for row in reader:
            key, name, native = (row[field].strip() for field in ("InChIKey", "Direct parent name", "Direct parent ID"))
            if not key or not name or not native or native.lower() == "null":
                continue
            if not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", key):
                raise ValueError("invalid ICD structure key")
            if not re.fullmatch(r"CHEMONTID:\d{7}", native):
                raise ValueError("invalid ICD direct-parent ID")
            rows.append((key, name, native))
    initialize(database)
    added = terms_added = 0
    with transaction(database) as connection:
        tool_id = record_tool_run(database, ToolRunInput("msfinder.icd_lookup", libraries=(library,),
                                  input_fingerprint=digest, parameters={"layout": "direct-parent.v1"}),
                                  connection=connection)
        connection.execute("INSERT OR IGNORE INTO class_ontology VALUES ('classyfire.chemont', 'unknown', ?)",
                           ("ICD direct parents; ontology version and hierarchy depth unreported",))
        for key, name, native in rows:
            term_id = make_id("class-term", "classyfire.chemont", "unknown", native)
            terms_added += connection.execute("""INSERT OR IGNORE INTO class_term VALUES
                (?, 'classyfire.chemont', 'unknown', ?, ?, ?)""", (term_id, LOOKUP_LEVEL, name, native)).rowcount
            added += connection.execute("INSERT OR IGNORE INTO structure_class_assignment VALUES (?, 'inchikey', ?, ?, 'lookup')",
                                        (key, term_id, tool_id)).rowcount
        total = connection.execute("SELECT COUNT(*) FROM structure_class_assignment WHERE tool_run_id=?", (tool_id,)).fetchone()[0]
    result = dict(library_name=path.name, sha256=digest, rows_read=len(rows), assignments_added=added,
                  terms_added=terms_added, assignments=total, parent_edges_added=0)
    run_path = Path(database).parent / "run.json"
    if update_run_json and run_path.exists():
        run = json.loads(run_path.read_text(encoding="utf-8-sig"))
        run["structure_classes"] = result
        run["db_bytes"] = Path(database).stat().st_size
        run["database_sha256"] = sha256_file(database)
        from .catalog_build import table_counts
        from .storage import connect
        with closing(connect(database)) as connection:
            run["table_counts"] = table_counts(connection)
        run_path.write_text(json.dumps(run, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result
