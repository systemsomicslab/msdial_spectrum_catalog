"""Shared allowlisted projection of annotator settings for linkage and export."""
import csv
from pathlib import Path


def read_annotator_settings(path):
    """Read TSV lines, retaining only the three named allowlisted columns.

    Other cells are discarded without interpretation, validation or hashing.
    Missing files yield no rows. Values are preserved for each caller's policy.
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
        for cells in reader:
            if not cells:
                continue
            rows.append({key: cells[index] if index < len(cells) else ""
                         for key, index in indices.items()})
    return rows
