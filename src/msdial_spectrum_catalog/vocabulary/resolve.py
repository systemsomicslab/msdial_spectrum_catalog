from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from ..identifiers import make_id
from .model import (
    AXES,
    CLAIM_AXIS,
    EVIDENCE_AXIS,
    TERM_STATUSES,
    ClaimReading,
    Term,
    Vocabulary,
)


REGISTRY_DIR = Path(__file__).resolve().parent / "registry"
VOCABULARY_FAMILY = "shin-massbank-annotation"
DEFAULT_VOCABULARY = "smb-v2-consensus"

MIGRATIONS_FILE = "migrations"
USE_CASES_FILE = "use_cases"
_NON_VOCABULARY_FILES = frozenset({MIGRATIONS_FILE, USE_CASES_FILE})


def _read(name: str) -> dict:
    path = REGISTRY_DIR / f"{name}.json"
    if not path.is_file():
        raise ValueError(f"Registry file {path.name} is not available in {REGISTRY_DIR}")
    return json.loads(path.read_text(encoding="utf-8"))


@lru_cache(maxsize=None)
def available_versions() -> tuple[str, ...]:
    return tuple(
        sorted(path.stem for path in REGISTRY_DIR.glob("*.json") if path.stem not in _NON_VOCABULARY_FILES)
    )


def _term(payload: dict) -> Term:
    status = str(payload["status"])
    if status not in TERM_STATUSES:
        raise ValueError(f"Term {payload.get('token')!r} has unsupported status {status!r}")
    rank = payload.get("specificity_rank")
    return Term(
        token=str(payload["token"]),
        concept_id=str(payload["concept_id"]),
        label=str(payload["label"]),
        definition=str(payload["definition"]),
        definition_source=str(payload["definition_source"]),
        status=status,
        specificity_rank=None if rank is None else int(rank),
        criteria_required=bool(payload["criteria_required"]),
        value_fields=tuple(str(value) for value in payload.get("value_fields", ())),
        notes=str(payload.get("notes", "")),
    )


def _axis(payload: dict, axis: str, version: str) -> tuple[Term, ...]:
    terms = tuple(_term(entry) for entry in payload[axis])
    tokens = [term.token for term in terms]
    concepts = [term.concept_id for term in terms]
    if len(set(tokens)) != len(tokens):
        raise ValueError(f"Vocabulary {version} repeats a token on axis {axis}")
    if len(set(concepts)) != len(concepts):
        raise ValueError(f"Vocabulary {version} repeats a concept id on axis {axis}")
    return terms


@lru_cache(maxsize=None)
def load_vocabulary(version: str = DEFAULT_VOCABULARY) -> Vocabulary:
    if version not in available_versions():
        raise ValueError(
            f"Unknown vocabulary version {version!r}; available versions are {', '.join(available_versions())}"
        )
    data = _read(version)
    if data.get("version") != version:
        raise ValueError(f"Registry file {version}.json declares version {data.get('version')!r}")
    expected_id = make_id("vocabulary", data.get("vocabulary_family", VOCABULARY_FAMILY), version)
    if data.get("vocabulary_id") != expected_id:
        raise ValueError(
            f"Registry file {version}.json declares vocabulary_id {data.get('vocabulary_id')!r} "
            f"but the identifier policy produces {expected_id!r}"
        )
    axes_payload = data["axes"]
    if set(axes_payload) != set(AXES):
        raise ValueError(f"Vocabulary {version} must declare exactly the axes {', '.join(AXES)}")
    return Vocabulary(
        vocabulary_id=expected_id,
        version=version,
        status=str(data["status"]),
        provenance=dict(data["provenance"]),
        axes={axis: _axis(axes_payload, axis, version) for axis in AXES},
        combination_rules=tuple(data.get("combination_rules", ())),
        decision_rules=tuple(data.get("decision_rules", ())),
        open_issues=tuple(data.get("open_issues", ())),
    )


@lru_cache(maxsize=None)
def load_migrations() -> tuple[dict, ...]:
    return tuple(_read(MIGRATIONS_FILE)["migrations"])


@lru_cache(maxsize=None)
def load_use_cases() -> dict:
    return _read(USE_CASES_FILE)


@lru_cache(maxsize=None)
def find_migration(from_version: str, to_version: str) -> dict:
    for migration in load_migrations():
        if migration["from_version"] == from_version and migration["to_version"] == to_version:
            _check_migration(migration)
            return migration
    raise ValueError(
        f"No migration is registered from {from_version!r} to {to_version!r}; "
        "a migration must be added to registry/migrations.json before readings can be moved"
    )


def _check_migration(migration: dict) -> None:
    source = load_vocabulary(migration["from_version"])
    target = load_vocabulary(migration["to_version"])
    for axis, token_map in migration.get("map", {}).items():
        for from_token, to_token in token_map.items():
            from_term = source.term(axis, from_token)
            to_term = target.term(axis, to_token)
            if from_term is None or to_term is None:
                raise ValueError(
                    f"Migration {migration['migration_id']} maps {axis} {from_token!r} -> {to_token!r} "
                    "but one of those tokens is not registered"
                )
            if from_term.concept_id != to_term.concept_id:
                raise ValueError(
                    f"Migration {migration['migration_id']} maps {axis} {from_token!r} -> {to_token!r} "
                    f"across different concepts ({from_term.concept_id} vs {to_term.concept_id})"
                )
    for token in migration.get("dropped_terms", ()):
        for axis in AXES:
            term = source.term(axis, token)
            if term is not None and target.by_concept(axis, term.concept_id) is not None:
                raise ValueError(
                    f"Migration {migration['migration_id']} declares {token!r} dropped "
                    f"but {term.concept_id} still exists in {migration['to_version']}"
                )


def _map_concepts(
    source: Vocabulary, target: Vocabulary, axis: str, concept_ids: tuple[str, ...]
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    kept: list[str] = []
    lost: list[str] = []
    for concept_id in concept_ids:
        if target.by_concept(axis, concept_id) is not None:
            kept.append(concept_id)
            continue
        source_term = source.by_concept(axis, concept_id)
        lost.append(source_term.token if source_term is not None else concept_id)
    return tuple(kept), tuple(lost)


def migrate_reading(reading: ClaimReading, to_version: str) -> ClaimReading:
    if to_version == reading.vocab_version:
        return reading
    source = load_vocabulary(reading.vocab_version)
    target = load_vocabulary(to_version)
    migration = find_migration(reading.vocab_version, to_version)
    claims, lost_claims = _map_concepts(source, target, CLAIM_AXIS, reading.claim_concept_ids)
    evidence, lost_evidence = _map_concepts(source, target, EVIDENCE_AXIS, reading.evidence_concept_ids)
    lost = lost_claims + lost_evidence
    declared = set(migration.get("dropped_terms", ()))
    undeclared = [token for token in lost if token not in declared]
    if undeclared:
        raise ValueError(
            f"Migration {migration['migration_id']} would silently drop {', '.join(undeclared)}; "
            "list the term in dropped_terms before migrating readings that carry it"
        )
    return ClaimReading(
        level=reading.level,
        claim_concept_ids=claims,
        evidence_concept_ids=evidence,
        vocab_version=to_version,
        notation_verbatim=reading.notation_verbatim,
        unknown_tokens=tuple(reading.unknown_tokens) + lost,
        unresolved=reading.unresolved or bool(lost),
    )


CANONICAL_EVIDENCE_ORDER: tuple[str, ...] = tuple(
    term.concept_id for term in load_vocabulary(DEFAULT_VOCABULARY).axes[EVIDENCE_AXIS]
)
