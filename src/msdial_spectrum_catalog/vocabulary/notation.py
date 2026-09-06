from __future__ import annotations

import re

from .model import (
    CLAIM_AXIS,
    EVIDENCE_AXIS,
    LEVEL_AXIS,
    PARSE_MODES,
    AmbiguousClaimTokenError,
    ClaimReading,
    NotationError,
    Vocabulary,
)
from .resolve import (
    CANONICAL_EVIDENCE_ORDER,
    available_versions,
    load_vocabulary,
    migrate_reading,
)


CLAIM_LEVEL = "L3"

_NOTATION = re.compile(r"^(?P<level>L[1-5])(?:-(?P<claim>[A-Z]{2}(?:\+[A-Z]{2})*))?(?:\[(?P<evidence>[A-Z,]*)\])?$")
# These three patterns require surrounding whitespace so that a hypothetical two-letter token is never
# mistaken for an English word.
_SLASH = re.compile(r"/")
_OR = re.compile(r"\s+or\s+", re.IGNORECASE)
_PROSE = re.compile(r";|\s+if\s+|\s+otherwise\b|\s+unless\s+", re.IGNORECASE)

_ALLOWED_STATUSES = {"strict": frozenset({"accepted"}), "permissive": frozenset({"accepted", "proposed"})}


def _reject_out_of_grammar(notation: str) -> None:
    if _SLASH.search(notation):
        raise NotationError(
            f"{notation!r} uses slash alternation; alternative claims must be recorded as separate "
            "candidate readings, not inside one notation string"
        )
    if _OR.search(notation):
        raise NotationError(
            f"{notation!r} uses 'or' alternation; alternative claims must be recorded as separate "
            "candidate readings, not inside one notation string"
        )
    if _PROSE.search(notation):
        raise NotationError(
            f"{notation!r} carries a prose conditional; the condition belongs in a criteria set and a "
            "curation comment, and the resulting decision belongs in structured fields"
        )


def _split_tokens(raw: str | None, separator: str, notation: str, what: str) -> tuple[str, ...]:
    if raw is None or raw == "":
        return ()
    tokens = raw.split(separator)
    if any(token == "" for token in tokens):
        raise NotationError(f"{notation!r} contains an empty {what} token")
    if len(set(tokens)) != len(tokens):
        raise NotationError(f"{notation!r} repeats a {what} token")
    return tuple(tokens)


def _split(notation: str) -> tuple[str, tuple[str, ...], tuple[str, ...]]:
    if not isinstance(notation, str) or not notation.strip():
        raise NotationError("An empty notation cannot be parsed")
    text = notation.strip()
    _reject_out_of_grammar(text)
    match = _NOTATION.match(text)
    if match is None:
        raise NotationError(
            f"{notation!r} does not match the grammar "
            "level ['-' claim {'+' claim}] ['[' [evidence {',' evidence}] ']']"
        )
    level = match.group("level")
    claims = _split_tokens(match.group("claim"), "+", notation, "claim")
    evidence = _split_tokens(match.group("evidence"), ",", notation, "evidence")
    if claims and level != CLAIM_LEVEL:
        raise NotationError(f"{notation!r} attaches a claim tag to {level}; only {CLAIM_LEVEL} accepts a claim tag")
    return level, claims, evidence


def resolve_claim_token(token: str, vocab_version: str | None = None) -> tuple[str, str]:
    """Resolve a claim token to (vocab_version, concept_id), refusing to guess a re-bound token."""
    if vocab_version is not None:
        term = load_vocabulary(vocab_version).term(CLAIM_AXIS, token)
        if term is None:
            raise NotationError(f"Claim token {token!r} is not registered in {vocab_version}")
        return vocab_version, term.concept_id
    candidates = []
    for version in available_versions():
        term = load_vocabulary(version).term(CLAIM_AXIS, token)
        if term is not None:
            candidates.append((version, term.concept_id))
    if not candidates:
        raise NotationError(f"Claim token {token!r} is not registered in any vocabulary version")
    if len({concept_id for _, concept_id in candidates}) > 1:
        raise AmbiguousClaimTokenError(token, candidates)
    return candidates[0]


def _resolve(
    vocabulary: Vocabulary, axis: str, tokens: tuple[str, ...], mode: str, notation: str
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    allowed = _ALLOWED_STATUSES.get(mode)
    concepts: list[str] = []
    unknown: list[str] = []
    for token in tokens:
        term = vocabulary.term(axis, token)
        if term is None:
            if mode == "strict":
                raise NotationError(
                    f"{notation!r} uses token {token!r}, which is not registered on axis {axis} "
                    f"in {vocabulary.version}"
                )
            unknown.append(token)
            continue
        if allowed is not None and term.status not in allowed:
            if mode == "strict":
                raise NotationError(
                    f"{notation!r} uses token {token!r}, whose status in {vocabulary.version} is "
                    f"{term.status!r} and not accepted"
                )
            unknown.append(token)
            continue
        concepts.append(term.concept_id)
    return tuple(concepts), tuple(unknown)


def parse_notation(notation: str, vocab_version: str, *, mode: str = "strict") -> ClaimReading:
    if mode not in PARSE_MODES:
        raise ValueError(f"Unknown parse mode {mode!r}; expected one of {', '.join(PARSE_MODES)}")
    if not vocab_version:
        raise NotationError(
            "A vocabulary version is required: the claim token CP is bound to different concepts by "
            "different versions, so a notation can never be resolved without one"
        )
    vocabulary = load_vocabulary(vocab_version)
    level, claim_tokens, evidence_tokens = _split(notation)
    if vocabulary.term(LEVEL_AXIS, level) is None:
        raise NotationError(f"{notation!r} uses level {level!r}, which is not registered in {vocab_version}")
    claims, unknown_claims = _resolve(vocabulary, CLAIM_AXIS, claim_tokens, mode, notation)
    evidence, unknown_evidence = _resolve(vocabulary, EVIDENCE_AXIS, evidence_tokens, mode, notation)
    return ClaimReading(
        level=level,
        claim_concept_ids=claims,
        evidence_concept_ids=evidence,
        vocab_version=vocab_version,
        notation_verbatim=notation,
        unknown_tokens=unknown_claims + unknown_evidence,
        unresolved=mode == "quarantine",
    )


def parse_notation_any(
    notation: str, candidate_versions: tuple[str, ...] | None = None, *, mode: str = "strict"
) -> tuple[ClaimReading, ...]:
    _split(notation)
    readings = []
    for version in candidate_versions if candidate_versions is not None else available_versions():
        try:
            readings.append(parse_notation(notation, version, mode=mode))
        except NotationError:
            continue
    return tuple(readings)


def _claim_sort_key(vocabulary: Vocabulary, concept_id: str) -> tuple[int, str]:
    term = vocabulary.by_concept(CLAIM_AXIS, concept_id)
    if term is None:
        raise NotationError(f"Claim concept {concept_id!r} is not registered in {vocabulary.version}")
    return (len(vocabulary.axes[CLAIM_AXIS]) if term.specificity_rank is None else term.specificity_rank, term.token)


def _evidence_sort_key(vocabulary: Vocabulary, concept_id: str) -> tuple[int, str]:
    term = vocabulary.by_concept(EVIDENCE_AXIS, concept_id)
    if term is None:
        raise NotationError(f"Evidence concept {concept_id!r} is not registered in {vocabulary.version}")
    order = len(CANONICAL_EVIDENCE_ORDER)
    if concept_id in CANONICAL_EVIDENCE_ORDER:
        order = CANONICAL_EVIDENCE_ORDER.index(concept_id)
    return (order, term.token)


def emit_notation(reading: ClaimReading, vocab_version: str | None = None, *, delimiter: str = ",") -> str:
    if delimiter == "|":
        raise ValueError("'|' cannot be used as an evidence delimiter because mzTab reads it as a multi-value separator")
    target = vocab_version or reading.vocab_version
    if target != reading.vocab_version:
        reading = migrate_reading(reading, target)
    vocabulary = load_vocabulary(target)
    claims = sorted(reading.claim_concept_ids, key=lambda concept: _claim_sort_key(vocabulary, concept))
    evidence = sorted(reading.evidence_concept_ids, key=lambda concept: _evidence_sort_key(vocabulary, concept))
    # Unknown tokens are never emitted: they carry no concept, so emitting them would republish an
    # unresolved token as if this version had accepted it.
    claim_part = ""
    if claims:
        claim_part = "-" + "+".join(vocabulary.by_concept(CLAIM_AXIS, concept).token for concept in claims)
    evidence_part = delimiter.join(vocabulary.by_concept(EVIDENCE_AXIS, concept).token for concept in evidence)
    return f"{reading.level}{claim_part}[{evidence_part}]"


def validate_combination(reading: ClaimReading) -> tuple[str, ...]:
    try:
        vocabulary = load_vocabulary(reading.vocab_version)
    except (ValueError, OSError):
        return ()
    claims = frozenset(reading.claim_concept_ids)
    if len(claims) < 2:
        return ()
    warnings = []
    described = False
    for rule in vocabulary.combination_rules:
        if frozenset(rule.get("claim_concept_ids", ())) != claims:
            continue
        described = True
        status = rule.get("status")
        if status in {"discouraged", "avoid"}:
            tokens = "+".join(rule.get("claim_tokens", ()))
            warnings.append(f"Claim combination {tokens} is {status} in {vocabulary.version}: {rule.get('note', '')}")
    if not described:
        tokens = []
        for concept_id in claims:
            term = vocabulary.by_concept(CLAIM_AXIS, concept_id)
            tokens.append(term.token if term is not None else concept_id)
        warnings.append(
            f"Claim combination {'+'.join(sorted(tokens))} is not described by the combination rules of "
            f"{vocabulary.version}; record a curation comment before publishing it"
        )
    return tuple(warnings)
