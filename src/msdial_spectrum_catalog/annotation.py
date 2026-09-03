from __future__ import annotations

import json
from dataclasses import dataclass, field

from .identifiers import make_id, short_hash
from .storage import connect, initialize, transaction
from .vocabulary import (
    CLAIM_AXIS,
    DEFAULT_VOCABULARY,
    EVIDENCE_AXIS,
    ClaimReading,
    emit_notation,
    load_vocabulary,
    parse_notation,
    validate_combination,
)

LEVELS = ("L1", "L2", "L3", "L4", "L5")
CLAIM_LEVEL = "L3"
MODES = ("strict", "permissive", "quarantine")
STRUCTURE_CLAIM = "smb:claim/structure"
SPECTRAL_LIBRARY_EVIDENCE = "smb:evidence/spectral_library"
IN_SILICO_LIBRARY_KIND = "in_silico_predicted"


@dataclass(frozen=True)
class ToolRunInput:
    tool_name: str
    tool_version: str | None = None
    provenance: dict | None = None
    parameters: dict | None = None
    input_fingerprint: str | None = None
    status: str = "completed"
    started_at: str | None = None
    finished_at: str | None = None


@dataclass(frozen=True)
class CandidateInput:
    rank: int
    compound_name: str | None = None
    formula: str | None = None
    inchikey: str | None = None
    smiles: str | None = None
    cxsmiles: str | None = None
    external_db_ref: str | None = None
    ontology: str | None = None
    cxsmiles_validated: bool = False
    score: float | None = None
    score_type: str | None = None
    score_gap_to_next: float | None = None
    rank_is_positional: bool = False
    reference_spectrum_id: str | None = None
    tool_run_id: str | None = None
    exclusion_status: str | None = None
    exclusion_reason: str | None = None


@dataclass(frozen=True)
class EvidenceInput:
    tag: str
    subtype: str | None = None
    metric: str | None = None
    measured_value: float | None = None
    measured_unit: str | None = None
    comparison: str | None = None
    threshold_value: float | None = None
    passed: bool | None = None
    value: dict | None = None
    source_uri: str | None = None
    source_spectrum_id: str | None = None
    source_reference_spectrum_id: str | None = None
    source_tool_run_id: str | None = None
    criteria_rule_id: str | None = None
    out_of_distribution: bool = False
    ood_reason: str | None = None


@dataclass
class AnnotationReport:
    run_id: str = ""
    assertions: int = 0
    candidates: int = 0
    evidence: int = 0
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def _canonical(payload: object) -> str:
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _json_or_none(payload: object | None) -> str | None:
    return None if payload is None else _canonical(payload)


def _flag(value: bool | None) -> int | None:
    return None if value is None else int(bool(value))


def _row_to_dict(row) -> dict:
    return {key: row[key] for key in row.keys()}


def _skeleton(inchikey: str | None) -> str | None:
    return inchikey.split("-")[0] if inchikey else None


def _allowed_statuses(mode: str) -> frozenset[str] | None:
    if mode == "strict":
        return frozenset({"accepted"})
    if mode == "permissive":
        return frozenset({"accepted", "proposed"})
    return None


def _resolve(vocabulary, axis: str, token: str, mode: str):
    term = vocabulary.term(axis, token)
    if term is None:
        if mode == "strict":
            raise ValueError(
                f"Unknown {axis} token '{token}' in vocabulary '{vocabulary.version}'"
            )
        return None
    statuses = _allowed_statuses(mode)
    if statuses is not None and term.status not in statuses:
        raise ValueError(
            f"{axis} token '{token}' has status '{term.status}' which mode '{mode}' does not allow"
        )
    return term


def record_tool_run(database, tool_run: ToolRunInput, *, run_id: str | None = None) -> str:
    initialize(database)
    identity = _canonical(
        {
            "tool_name": tool_run.tool_name,
            "tool_version": tool_run.tool_version,
            "run_id": run_id,
            "parameters": tool_run.parameters,
            "provenance": tool_run.provenance,
            "input_fingerprint": tool_run.input_fingerprint,
        }
    )
    tool_run_id = make_id(
        "tool-run", tool_run.tool_name, tool_run.tool_version or "unversioned", short_hash(identity)
    )
    with transaction(database) as connection:
        connection.execute(
            """INSERT INTO annotation_tool_run(
                   tool_run_id, run_id, tool_name, tool_version, tool_provenance_json,
                   parameters_json, input_fingerprint, started_at, finished_at, status)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
               ON CONFLICT(tool_run_id) DO UPDATE SET
                   run_id = excluded.run_id,
                   tool_version = excluded.tool_version,
                   tool_provenance_json = excluded.tool_provenance_json,
                   parameters_json = excluded.parameters_json,
                   input_fingerprint = excluded.input_fingerprint,
                   started_at = excluded.started_at,
                   finished_at = excluded.finished_at,
                   status = excluded.status""",
            (
                tool_run_id,
                run_id,
                tool_run.tool_name,
                tool_run.tool_version,
                _json_or_none(tool_run.provenance),
                _json_or_none(tool_run.parameters),
                tool_run.input_fingerprint,
                tool_run.started_at,
                tool_run.finished_at,
                tool_run.status,
            ),
        )
    return tool_run_id


def record_criteria_set(
    database,
    criteria_set_id: str,
    name: str,
    *,
    version: str | None = None,
    description: str | None = None,
    vocab_version: str = DEFAULT_VOCABULARY,
    rules: list[dict] | None = None,
) -> str:
    initialize(database)
    vocabulary = load_vocabulary(vocab_version)
    rule_list = list(rules or [])
    resolved: list[tuple[str, str, str, dict]] = []
    seen: set[tuple[str, str]] = set()
    for rule in rule_list:
        tag = rule.get("tag")
        if not tag:
            raise ValueError("Every criteria rule needs a 'tag'")
        term = _resolve(vocabulary, EVIDENCE_AXIS, str(tag), "permissive")
        if term is None:
            raise ValueError(
                f"Unknown evidence token '{tag}' in vocabulary '{vocab_version}'; "
                "a criteria rule must resolve to a concept"
            )
        metric = str(rule.get("metric") or "")
        key = (term.concept_id, metric)
        if key in seen:
            raise ValueError(
                f"Duplicate criteria rule for evidence '{tag}' and metric '{metric}'"
            )
        seen.add(key)
        resolved.append((term.concept_id, str(tag), metric, rule))

    with transaction(database) as connection:
        connection.execute(
            """INSERT INTO criteria_set(criteria_set_id, name, version, description, rules_json, vocab_version)
               VALUES (?, ?, ?, ?, ?, ?)
               ON CONFLICT(criteria_set_id) DO UPDATE SET
                   name = excluded.name,
                   version = excluded.version,
                   description = excluded.description,
                   rules_json = excluded.rules_json,
                   vocab_version = excluded.vocab_version""",
            (criteria_set_id, name, version, description, _canonical(rule_list), vocab_version),
        )
        connection.execute("DELETE FROM criteria_rule WHERE criteria_set_id = ?", (criteria_set_id,))
        for concept_id, tag, metric, rule in resolved:
            connection.execute(
                """INSERT INTO criteria_rule(
                       criteria_rule_id, criteria_set_id, evidence_concept_id, evidence_token,
                       operational_criterion, metric, comparison, threshold_value, threshold_unit,
                       scope_json, example, notes)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    make_id("criteria-rule", criteria_set_id, concept_id, metric),
                    criteria_set_id,
                    concept_id,
                    tag,
                    rule.get("operational_criterion"),
                    metric,
                    rule.get("comparison"),
                    rule.get("threshold_value"),
                    rule.get("threshold_unit"),
                    _json_or_none(rule.get("scope")),
                    rule.get("example"),
                    rule.get("notes"),
                ),
            )
    return criteria_set_id


def criteria_rule_id_for(
    criteria_set_id: str, concept_id: str, metric: str | None = None
) -> str:
    return make_id("criteria-rule", criteria_set_id, concept_id, metric or "")


def record_assertion(
    database,
    *,
    spectrum_id: str,
    level: str,
    claim_tokens: tuple[str, ...] = (),
    evidence: list[EvidenceInput] | None = None,
    candidates: list[CandidateInput] | None = None,
    vocab_version: str = DEFAULT_VOCABULARY,
    criteria_set_id: str | None = None,
    tool_run_id: str | None = None,
    ambiguity_class_id: str | None = None,
    curation_comment: str | None = None,
    compound_name: str | None = None,
    formula: str | None = None,
    structure_id: str | None = None,
    alignment_feature_id: str | None = None,
    notation: str | None = None,
    mode: str = "strict",
) -> str:
    if mode not in MODES:
        raise ValueError(f"Unknown mode '{mode}'; expected one of {MODES}")
    if level not in LEVELS:
        raise ValueError(f"Unknown annotation level '{level}'; expected one of {LEVELS}")
    claim_list = tuple(claim_tokens)
    if claim_list and level != CLAIM_LEVEL:
        raise ValueError(
            f"Claim tokens {claim_list} are only legal at level {CLAIM_LEVEL}, not '{level}'"
        )
    if len(set(claim_list)) != len(claim_list):
        raise ValueError(f"Duplicate claim token in {claim_list}")

    evidence_list = list(evidence or [])
    candidate_list = list(candidates or [])
    ranks = [candidate.rank for candidate in candidate_list]
    if sorted(ranks) != list(range(1, len(ranks) + 1)):
        raise ValueError(
            f"Candidate ranks must be a contiguous 1..{len(ranks)} with no duplicates, got {sorted(ranks)}"
        )

    initialize(database)
    vocabulary = load_vocabulary(vocab_version)

    unknown_tokens: list[str] = []
    claim_terms = []
    for token in claim_list:
        term = _resolve(vocabulary, CLAIM_AXIS, token, mode)
        if term is None:
            unknown_tokens.append(token)
        else:
            claim_terms.append(term)
    # Most specific first: specificity_rank 1 is the strongest structural claim, and the primary
    # concept written onto the assertion must be that one, not whichever the caller listed first.
    claim_terms.sort(key=lambda term: (term.specificity_rank or 0, term.token))
    claim_concept_ids = tuple(term.concept_id for term in claim_terms)
    primary_claim = claim_concept_ids[0] if claim_concept_ids else None

    evidence_terms = []
    for item in evidence_list:
        term = _resolve(vocabulary, EVIDENCE_AXIS, item.tag, mode)
        if term is None:
            unknown_tokens.append(item.tag)
        evidence_terms.append(term)
    evidence_concept_ids: list[str] = []
    for term in evidence_terms:
        if term is not None and term.concept_id not in evidence_concept_ids:
            evidence_concept_ids.append(term.concept_id)

    reading_unresolved = bool(unknown_tokens) and mode == "quarantine"
    claim_unresolved = reading_unresolved or (level == CLAIM_LEVEL and not claim_list)

    reading = ClaimReading(
        level=level,
        claim_concept_ids=claim_concept_ids,
        evidence_concept_ids=tuple(evidence_concept_ids),
        vocab_version=vocab_version,
        notation_verbatim=notation or "",
        unknown_tokens=tuple(unknown_tokens),
        unresolved=reading_unresolved,
    )
    canonical_notation = emit_notation(reading, vocab_version)
    if notation is None:
        notation_verbatim = canonical_notation
    else:
        parsed = parse_notation(notation, vocab_version, mode=mode)
        if parsed.level != level:
            raise ValueError(
                f"Notation '{notation}' declares level '{parsed.level}' but level '{level}' was recorded"
            )
        if set(parsed.claim_concept_ids) != set(claim_concept_ids):
            raise ValueError(
                f"Notation '{notation}' declares claims {parsed.claim_concept_ids} "
                f"but {claim_concept_ids} were recorded"
            )
        if set(parsed.evidence_concept_ids) != set(evidence_concept_ids):
            raise ValueError(
                f"Notation '{notation}' declares evidence {parsed.evidence_concept_ids} "
                f"but {tuple(evidence_concept_ids)} were recorded"
            )
        notation_verbatim = notation

    claim_slug = "+".join(term.token for term in claim_terms) or "none"
    identity = _canonical(
        {
            "spectrum_id": spectrum_id,
            "level": level,
            "claims": list(claim_concept_ids) + sorted(unknown_tokens),
            "vocab_version": vocab_version,
            "subject": structure_id or compound_name or formula or "",
            "tool_run_id": tool_run_id or "",
        }
    )
    assertion_id = make_id("assertion", spectrum_id, level, claim_slug, short_hash(identity))
    subject_kind = "alignment_feature" if alignment_feature_id else "spectrum"
    annotation_claim = "+".join(claim_list) if claim_list else None

    with transaction(database) as connection:
        connection.execute(
            """INSERT INTO annotation_assertion(
                   assertion_id, spectrum_id, annotation_level, annotation_claim, compound_name,
                   formula, structure_id, criteria_set_id, curation_comment, subject_kind,
                   alignment_feature_id, vocab_version, claim_concept_id, notation_verbatim,
                   claim_unresolved, candidate_count, ambiguity_class_id, tool_run_id)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
               ON CONFLICT(assertion_id) DO UPDATE SET
                   annotation_level = excluded.annotation_level,
                   annotation_claim = excluded.annotation_claim,
                   compound_name = excluded.compound_name,
                   formula = excluded.formula,
                   structure_id = excluded.structure_id,
                   criteria_set_id = excluded.criteria_set_id,
                   curation_comment = excluded.curation_comment,
                   subject_kind = excluded.subject_kind,
                   alignment_feature_id = excluded.alignment_feature_id,
                   vocab_version = excluded.vocab_version,
                   claim_concept_id = excluded.claim_concept_id,
                   notation_verbatim = excluded.notation_verbatim,
                   claim_unresolved = excluded.claim_unresolved,
                   candidate_count = excluded.candidate_count,
                   ambiguity_class_id = excluded.ambiguity_class_id,
                   tool_run_id = excluded.tool_run_id""",
            (
                assertion_id,
                spectrum_id,
                level,
                annotation_claim,
                compound_name,
                formula,
                structure_id,
                criteria_set_id,
                curation_comment,
                subject_kind,
                alignment_feature_id,
                vocab_version,
                primary_claim,
                notation_verbatim,
                int(claim_unresolved),
                len(candidate_list),
                ambiguity_class_id,
                tool_run_id,
            ),
        )
        connection.execute(
            "DELETE FROM annotation_claim_component WHERE assertion_id = ?", (assertion_id,)
        )
        connection.execute("DELETE FROM annotation_evidence WHERE assertion_id = ?", (assertion_id,))
        connection.execute("DELETE FROM annotation_candidate WHERE assertion_id = ?", (assertion_id,))

        for ordinal, term in enumerate(claim_terms, start=1):
            connection.execute(
                """INSERT INTO annotation_claim_component(
                       claim_component_id, assertion_id, ordinal, claim_concept_id)
                   VALUES (?, ?, ?, ?)""",
                (
                    make_id("claim-component", assertion_id, ordinal, term.token),
                    assertion_id,
                    ordinal,
                    term.concept_id,
                ),
            )

        for ordinal, (item, term) in enumerate(zip(evidence_list, evidence_terms), start=1):
            connection.execute(
                """INSERT INTO annotation_evidence(
                       evidence_id, assertion_id, evidence_tag, evidence_value_json, source_uri,
                       evidence_concept_id, evidence_subtype, metric, measured_value, measured_unit,
                       comparison, threshold_value, passed, criteria_rule_id, source_spectrum_id,
                       source_reference_spectrum_id, source_tool_run_id, out_of_distribution, ood_reason)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    make_id("evidence", assertion_id, ordinal, item.tag),
                    assertion_id,
                    item.tag,
                    _json_or_none(item.value),
                    item.source_uri,
                    None if term is None else term.concept_id,
                    item.subtype,
                    item.metric,
                    item.measured_value,
                    item.measured_unit,
                    item.comparison,
                    item.threshold_value,
                    _flag(item.passed),
                    item.criteria_rule_id,
                    item.source_spectrum_id,
                    item.source_reference_spectrum_id,
                    item.source_tool_run_id,
                    int(bool(item.out_of_distribution)),
                    item.ood_reason,
                ),
            )

        for candidate in candidate_list:
            connection.execute(
                """INSERT INTO annotation_candidate(
                       candidate_id, assertion_id, rank, rank_is_positional, compound_name, formula,
                       inchikey, inchikey_skeleton, smiles, cxsmiles, cxsmiles_validated,
                       external_db_ref, ontology, score, score_type, score_gap_to_next, tool_run_id,
                       reference_spectrum_id, exclusion_status, exclusion_reason)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    make_id("candidate", assertion_id, candidate.rank),
                    assertion_id,
                    candidate.rank,
                    int(bool(candidate.rank_is_positional)),
                    candidate.compound_name,
                    candidate.formula,
                    candidate.inchikey,
                    _skeleton(candidate.inchikey),
                    candidate.smiles,
                    candidate.cxsmiles,
                    int(bool(candidate.cxsmiles_validated)),
                    candidate.external_db_ref,
                    candidate.ontology,
                    candidate.score,
                    candidate.score_type,
                    candidate.score_gap_to_next,
                    candidate.tool_run_id or tool_run_id,
                    candidate.reference_spectrum_id,
                    candidate.exclusion_status,
                    candidate.exclusion_reason,
                ),
            )
    return assertion_id


def _reading_from_row(connection, row) -> ClaimReading:
    claim_concept_ids = tuple(
        item[0]
        for item in connection.execute(
            "SELECT claim_concept_id FROM annotation_claim_component WHERE assertion_id = ? ORDER BY ordinal",
            (row["assertion_id"],),
        )
    )
    evidence_concept_ids: list[str] = []
    for item in connection.execute(
        "SELECT evidence_concept_id FROM annotation_evidence WHERE assertion_id = ? ORDER BY evidence_id",
        (row["assertion_id"],),
    ):
        if item[0] and item[0] not in evidence_concept_ids:
            evidence_concept_ids.append(item[0])
    return ClaimReading(
        level=row["annotation_level"],
        claim_concept_ids=claim_concept_ids,
        evidence_concept_ids=tuple(evidence_concept_ids),
        vocab_version=row["vocab_version"] or DEFAULT_VOCABULARY,
        notation_verbatim=row["notation_verbatim"] or "",
        unknown_tokens=(),
        unresolved=bool(row["claim_unresolved"]),
    )


def notation_for(database, assertion_id: str, *, vocab_version: str | None = None) -> str:
    connection = connect(database)
    try:
        row = connection.execute(
            "SELECT * FROM annotation_assertion WHERE assertion_id = ?", (assertion_id,)
        ).fetchone()
        if row is None:
            raise ValueError(f"Assertion '{assertion_id}' does not exist")
        if row["vocab_version"] is None and vocab_version is None:
            raise ValueError(f"Assertion '{assertion_id}' has no stored vocab_version")
        return emit_notation(_reading_from_row(connection, row), vocab_version)
    finally:
        connection.close()


def load_assertion(database, assertion_id: str) -> dict | None:
    connection = connect(database)
    try:
        row = connection.execute(
            "SELECT * FROM annotation_assertion WHERE assertion_id = ?", (assertion_id,)
        ).fetchone()
        if row is None:
            return None
        result = _row_to_dict(row)
        result["claim_components"] = [
            _row_to_dict(item)
            for item in connection.execute(
                "SELECT * FROM annotation_claim_component WHERE assertion_id = ? ORDER BY ordinal",
                (assertion_id,),
            )
        ]
        result["evidence"] = [
            _row_to_dict(item)
            for item in connection.execute(
                "SELECT * FROM annotation_evidence WHERE assertion_id = ? ORDER BY evidence_id",
                (assertion_id,),
            )
        ]
        result["candidates"] = [
            _row_to_dict(item)
            for item in connection.execute(
                "SELECT * FROM annotation_candidate WHERE assertion_id = ? ORDER BY rank",
                (assertion_id,),
            )
        ]
        result["notation"] = emit_notation(_reading_from_row(connection, row), row["vocab_version"])
        return result
    finally:
        connection.close()


def list_assertions(database, *, run_id=None, alignment_feature_id=None, spectrum_id=None) -> list[dict]:
    clauses: list[str] = []
    parameters: list[str] = []
    if run_id is not None:
        clauses.append(
            "(a.spectrum_id IN (SELECT spectrum_id FROM spectrum WHERE run_id = ?)"
            " OR a.alignment_feature_id IN (SELECT alignment_feature_id FROM alignment_feature WHERE run_id = ?))"
        )
        parameters.extend([run_id, run_id])
    if alignment_feature_id is not None:
        clauses.append("a.alignment_feature_id = ?")
        parameters.append(alignment_feature_id)
    if spectrum_id is not None:
        clauses.append("a.spectrum_id = ?")
        parameters.append(spectrum_id)
    where = f" WHERE {' AND '.join(clauses)}" if clauses else ""
    connection = connect(database)
    try:
        rows = connection.execute(
            f"SELECT a.* FROM annotation_assertion a{where} ORDER BY a.assertion_id", parameters
        ).fetchall()
        results = []
        for row in rows:
            record = _row_to_dict(row)
            record["evidence_count"] = connection.execute(
                "SELECT COUNT(*) FROM annotation_evidence WHERE assertion_id = ?", (row["assertion_id"],)
            ).fetchone()[0]
            record["claim_component_count"] = connection.execute(
                "SELECT COUNT(*) FROM annotation_claim_component WHERE assertion_id = ?",
                (row["assertion_id"],),
            ).fetchone()[0]
            record["notation"] = emit_notation(_reading_from_row(connection, row), row["vocab_version"])
            results.append(record)
        return results
    finally:
        connection.close()


def validate_annotations(database, run_id: str) -> AnnotationReport:
    report = AnnotationReport(run_id=run_id)
    connection = connect(database)
    try:
        if connection.execute("SELECT 1 FROM analysis_run WHERE run_id = ?", (run_id,)).fetchone() is None:
            report.errors.append("Run does not exist")
            return report
        rows = connection.execute(
            """SELECT a.* FROM annotation_assertion a
               WHERE a.spectrum_id IN (SELECT spectrum_id FROM spectrum WHERE run_id = ?)
                  OR a.alignment_feature_id IN (SELECT alignment_feature_id FROM alignment_feature WHERE run_id = ?)
               ORDER BY a.assertion_id""",
            (run_id, run_id),
        ).fetchall()
        report.assertions = len(rows)
        vocabularies: dict[str, object] = {}
        for row in rows:
            assertion_id = row["assertion_id"]
            if connection.execute(
                "SELECT 1 FROM spectrum WHERE spectrum_id = ?", (row["spectrum_id"],)
            ).fetchone() is None:
                report.errors.append(
                    f"{assertion_id}: spectrum '{row['spectrum_id']}' is missing from spectrum"
                )

            version = row["vocab_version"]
            vocabulary = None
            if not version:
                report.errors.append(f"{assertion_id}: no vocab_version is stored")
            else:
                if version not in vocabularies:
                    try:
                        vocabularies[version] = load_vocabulary(version)
                    except Exception as error:
                        vocabularies[version] = None
                        report.errors.append(f"{assertion_id}: vocabulary '{version}' is unusable ({error})")
                vocabulary = vocabularies[version]

            components = connection.execute(
                "SELECT ordinal, claim_concept_id FROM annotation_claim_component WHERE assertion_id = ? ORDER BY ordinal",
                (assertion_id,),
            ).fetchall()
            if row["annotation_level"] == CLAIM_LEVEL and not components:
                report.warnings.append(
                    f"{assertion_id}: level {CLAIM_LEVEL} carries no claim component, so the claim is incomplete"
                )
            claim_concept_ids = tuple(item["claim_concept_id"] for item in components)
            if vocabulary is not None:
                for concept_id in claim_concept_ids:
                    if vocabulary.by_concept(CLAIM_AXIS, concept_id) is None:
                        report.errors.append(
                            f"{assertion_id}: claim concept '{concept_id}' is absent from vocabulary '{version}'"
                        )

            candidate_ranks = [
                item[0]
                for item in connection.execute(
                    "SELECT rank FROM annotation_candidate WHERE assertion_id = ? ORDER BY rank",
                    (assertion_id,),
                )
            ]
            report.candidates += len(candidate_ranks)
            if candidate_ranks and candidate_ranks != list(range(1, len(candidate_ranks) + 1)):
                report.errors.append(
                    f"{assertion_id}: candidate ranks are not contiguous 1..N: {candidate_ranks}"
                )
            if len(candidate_ranks) != row["candidate_count"]:
                report.warnings.append(
                    f"{assertion_id}: candidate_count is {row['candidate_count']} "
                    f"but {len(candidate_ranks)} candidate rows are stored"
                )
            if (
                row["candidate_count"] > 1
                and row["claim_concept_id"] == STRUCTURE_CLAIM
                and not row["ambiguity_class_id"]
            ):
                report.warnings.append(
                    f"{assertion_id}: {row['candidate_count']} candidates carry a single-structure claim "
                    "without an ambiguity_class_id, so the 'A or B' has no stated ambiguity analysis"
                )

            evidence_rows = connection.execute(
                "SELECT * FROM annotation_evidence WHERE assertion_id = ? ORDER BY evidence_id",
                (assertion_id,),
            ).fetchall()
            report.evidence += len(evidence_rows)
            evidence_concept_ids: list[str] = []
            for item in evidence_rows:
                if item["evidence_concept_id"] is None:
                    report.warnings.append(
                        f"{assertion_id}: evidence '{item['evidence_tag']}' has no resolved concept"
                    )
                elif item["evidence_concept_id"] not in evidence_concept_ids:
                    evidence_concept_ids.append(item["evidence_concept_id"])
                if item["passed"] == 0:
                    report.warnings.append(
                        f"{assertion_id}: evidence '{item['evidence_tag']}' did not pass its criterion "
                        "but is asserted anyway"
                    )
                if (
                    item["evidence_concept_id"] == SPECTRAL_LIBRARY_EVIDENCE
                    and item["source_reference_spectrum_id"]
                ):
                    library_kind = connection.execute(
                        """SELECT rl.library_kind FROM reference_spectrum rs
                           JOIN reference_library rl USING(library_id)
                           WHERE rs.reference_spectrum_id = ?""",
                        (item["source_reference_spectrum_id"],),
                    ).fetchone()
                    if library_kind is not None and library_kind[0] == IN_SILICO_LIBRARY_KIND:
                        report.errors.append(
                            f"{assertion_id}: evidence '{item['evidence_tag']}' cites reference spectrum "
                            f"'{item['source_reference_spectrum_id']}' from an {IN_SILICO_LIBRARY_KIND} "
                            "library, which must never be reported as spectral-library evidence"
                        )

            if vocabulary is not None:
                reading = ClaimReading(
                    level=row["annotation_level"],
                    claim_concept_ids=claim_concept_ids,
                    evidence_concept_ids=tuple(evidence_concept_ids),
                    vocab_version=version,
                    notation_verbatim=row["notation_verbatim"] or "",
                    unknown_tokens=(),
                    unresolved=bool(row["claim_unresolved"]),
                )
                for warning in validate_combination(reading):
                    report.warnings.append(f"{assertion_id}: {warning}")
        return report
    finally:
        connection.close()
