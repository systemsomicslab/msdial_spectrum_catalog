"""Shared allowlisted projection of annotator settings for linkage and export."""
import csv
import re
from pathlib import Path


def read_annotator_settings(path):
    """Read TSV lines, retaining only the three named allowlisted columns.

    Other cells are discarded without interpretation, validation or hashing.
    Missing files yield no rows. Selected values are checked for paths; tier
    and a default acquired kind are derived from the selected evidence fields.
    """
    path = Path(path)
    if not path.exists():
        return []
    allowed = ("annotator_id", "evidence_tier", "library_kind")
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, [])
        indices = {key: header.index(key) for key in allowed if key in header}
        if "annotator_id" not in indices:
            raise ValueError("settings lack annotator_id")
        for cells in reader:
            if not cells:
                continue
            row = {key: cells[index] if index < len(cells) else ""
                   for key, index in indices.items()}
            if any(re.search(r"^[A-Za-z]:[\\/]|^\\\\", value)
                   or (re.search(r"\.msp|\.lbm", value, re.I) and re.search(r"[\\/]", value))
                   for value in row.values()):
                raise ValueError("path-like value in allowlisted settings field")
            evidence, kind = row.get("evidence_tier", ""), row.get("library_kind", "")
            tier = ("in_silico" if kind.lower() == "predicted" else "strict" if "high" in evidence.lower()
                    else "relaxed" if "low" in evidence.lower() else "unknown")
            row["tier"] = tier
            row["library_kind"] = kind or ("acquired" if tier in {"strict", "relaxed"} else "")
            rows.append(row)
    return rows
