from __future__ import annotations

from dataclasses import dataclass


LEVEL_AXIS = "annotation_level"
CLAIM_AXIS = "level3_claim"
EVIDENCE_AXIS = "evidence"
AXES = (LEVEL_AXIS, CLAIM_AXIS, EVIDENCE_AXIS)

TERM_STATUSES = frozenset({"accepted", "proposed", "superseded"})
PARSE_MODES = ("strict", "permissive", "quarantine")


class NotationError(ValueError):
    """Raised when a notation string is outside the proposal grammar or its tokens do not resolve."""


class AmbiguousClaimTokenError(NotationError):
    def __init__(self, token: str, candidates) -> None:
        self.token = token
        self.candidates = tuple(candidates)
        rendered = "; ".join(f"{version} -> {concept_id}" for version, concept_id in self.candidates)
        super().__init__(
            f"Claim token {token!r} is bound to different concepts by different vocabulary versions "
            f"({rendered}); an explicit vocab_version is required"
        )


@dataclass(frozen=True)
class Term:
    token: str
    concept_id: str
    label: str
    definition: str
    definition_source: str
    status: str
    specificity_rank: int | None
    criteria_required: bool
    value_fields: tuple[str, ...]
    notes: str


@dataclass(frozen=True)
class Vocabulary:
    vocabulary_id: str
    version: str
    status: str
    provenance: dict
    axes: dict[str, tuple[Term, ...]]
    combination_rules: tuple[dict, ...]
    decision_rules: tuple[dict, ...]
    open_issues: tuple[dict, ...]

    def term(self, axis: str, token: str) -> Term | None:
        for term in self.axes.get(axis, ()):
            if term.token == token:
                return term
        return None

    def by_concept(self, axis: str, concept_id: str) -> Term | None:
        for term in self.axes.get(axis, ()):
            if term.concept_id == concept_id:
                return term
        return None


@dataclass(frozen=True)
class ClaimReading:
    level: str
    claim_concept_ids: tuple[str, ...]
    evidence_concept_ids: tuple[str, ...]
    vocab_version: str
    notation_verbatim: str
    unknown_tokens: tuple[str, ...] = ()
    unresolved: bool = False
