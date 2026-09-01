from __future__ import annotations

import hashlib
from urllib.parse import quote


def _part(value: object) -> str:
    return quote(str(value).strip(), safe="._-")


def make_id(kind: str, *parts: object) -> str:
    return "urn:msdial:" + _part(kind) + ":" + ":".join(_part(part) for part in parts)


def short_hash(value: str, length: int = 16) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:length]

