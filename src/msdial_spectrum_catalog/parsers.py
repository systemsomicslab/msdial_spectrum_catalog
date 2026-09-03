from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


COMMENT_TOKEN = re.compile(r"(?:^|\|)([A-Za-z0-9_]+)=([^|]*)")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_number(value: str | None, kind=float):
    if value is None or not value.strip() or value.strip().lower() == "null":
        return None
    try:
        return kind(value)
    except (TypeError, ValueError):
        return None


def fold_headers(row: dict[str, str]) -> dict[str, str]:
    """Index one table row by case-folded header."""
    return {str(key).strip().casefold(): (value or "").strip() for key, value in row.items()}


def folded_value(folded: dict[str, str], *headers: str) -> str | None:
    """Return the first header present in a folded row, or None when absent, empty or 'null'."""
    for header in headers:
        found = folded.get(header.strip().casefold())
        if found is None:
            continue
        if not found or found.lower() == "null":
            return None
        return found
    return None


def split_refs(value: str | None) -> list[str]:
    """Split an mzTab-M ID reference list, which MS-DIAL writes '|'- or ','-separated."""
    if not value or value.strip().lower() == "null":
        return []
    return [piece.strip() for piece in value.replace("|", ",").split(",") if piece.strip()]


def read_tsv(path: Path, header_prefix: str | None = None) -> Iterator[tuple[int, dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        header_line = 1
        if header_prefix:
            for header_line, line in enumerate(handle, start=1):
                if line.startswith(header_prefix):
                    header = line.rstrip("\r\n").split("\t")
                    break
            else:
                raise ValueError(f"Header {header_prefix!r} was not found in {path}")
            reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        else:
            reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=header_line + 1):
            yield row_number, {str(key): value or "" for key, value in row.items() if key is not None}


def read_analysis_csv(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        result = {}
        for row in reader:
            normalized = {str(key).strip().lower(): value or "" for key, value in row.items() if key is not None}
            name = normalized.get("file_name", "").strip()
            if name:
                result[name] = normalized
        return result


@dataclass(frozen=True)
class MspRecord:
    index: int
    fields: dict[str, str]
    comment_tokens: dict[str, str]
    peaks: list[tuple[float, float]]


def read_msp(path: Path) -> Iterator[MspRecord]:
    fields: dict[str, str] = {}
    peaks: list[tuple[float, float]] = []
    reading_peaks = False
    record_index = 0
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                if fields:
                    comment = fields.get("COMMENT", "")
                    yield MspRecord(record_index, fields, dict(COMMENT_TOKEN.findall(comment)), peaks)
                    record_index += 1
                    fields, peaks, reading_peaks = {}, [], False
                continue
            if reading_peaks:
                pieces = line.replace("\t", " ").split()
                if len(pieces) >= 2:
                    mz = parse_number(pieces[0])
                    intensity = parse_number(pieces[1])
                    if mz is not None and intensity is not None:
                        peaks.append((mz, intensity))
                continue
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            fields[key.strip().upper()] = value.strip()
            if key.strip().upper().replace(" ", "") == "NUMPEAKS":
                reading_peaks = True
        if fields:
            comment = fields.get("COMMENT", "")
            yield MspRecord(record_index, fields, dict(COMMENT_TOKEN.findall(comment)), peaks)


def read_mztab(path: Path) -> Iterator[tuple[int, str, dict[str, str]]]:
    headers: dict[str, list[str]] = {}
    header_types = {"SMH": "SML", "SFH": "SMF", "SEH": "SME"}
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for row_number, raw in enumerate(handle, start=1):
            columns = raw.rstrip("\r\n").split("\t")
            if not columns:
                continue
            row_type = columns[0]
            if row_type in header_types:
                headers[header_types[row_type]] = columns[1:]
                continue
            if row_type not in {"SML", "SMF", "SME"} or row_type not in headers:
                continue
            yield row_number, row_type, dict(zip(headers[row_type], columns[1:]))
