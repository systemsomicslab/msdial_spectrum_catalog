"""Tool-linked class predictions, shared by all class-producing adapters."""
from __future__ import annotations

from contextlib import closing

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

from .identifiers import make_id
from .storage import connect, initialize, transaction

LEVELS = {"classyfire.chemont": ("kingdom", "superclass", "class", "subclass") + tuple(
    "level" + str(i) for i in range(5, 12)), "npclassifier": ("pathway", "superclass", "class")}

# ICD reports only the direct parent, with no depth. Do not invent a hierarchy.
LOOKUP_LEVEL = "direct_parent"


@dataclass(frozen=True)
class ClassTermInput:
    ontology_id: str
    ontology_version: str
    level: str
    name: str
    native_id: str | None = None


@dataclass(frozen=True)
class ClassPredictionInput:
    subject_kind: str
    subject_id: str
    lineage: Sequence[ClassTermInput]
    posteriors: Mapping[str, float]
    formula_context: str | None = None
    formula_rank: int | None = None
    rank: int = 1
    source_row_ref: str | None = None


def _validate_lineage(terms):
    previous = -1
    ontology = None
    for term in terms:
        levels = LEVELS.get(term.ontology_id, ())
        if term.level not in levels or not term.name or not term.ontology_version:
            raise ValueError("unknown class level, empty name or missing explicit ontology version")
        identity = (term.ontology_id, term.ontology_version)
        if ontology is not None and ontology != identity:
            raise ValueError("a lineage cannot mix ontologies or versions")
        index = levels.index(term.level)
        if index <= previous:
            raise ValueError("lineage is out of level order")
        previous, ontology = index, identity


def _upsert_lineage(connection, terms):
    _validate_lineage(terms)
    ids = []
    for term in terms:
        connection.execute("INSERT OR IGNORE INTO class_ontology(ontology_id, ontology_version) VALUES (?, ?)",
                           (term.ontology_id, term.ontology_version))
        term_id = make_id("class-term", term.ontology_id, term.ontology_version,
                          term.native_id or term.name.casefold())
        connection.execute("INSERT OR IGNORE INTO class_term VALUES (?, ?, ?, ?, ?, ?)",
                           (term_id, term.ontology_id, term.ontology_version, term.level, term.name, term.native_id))
        row = connection.execute("SELECT * FROM class_term WHERE term_id=?", (term_id,)).fetchone()
        if row is None or (row["level"], row["name"], row["native_id"]) != (term.level, term.name, term.native_id):
            raise ValueError("conflicting class term identity")
        if ids:
            connection.execute("INSERT OR IGNORE INTO class_term_parent VALUES (?, ?)", (term_id, ids[-1]))
        ids.append(term_id)
    return ids


def upsert_lineage(database, terms: Sequence[ClassTermInput]) -> list[str]:
    _validate_lineage(terms)
    initialize(database)
    with transaction(database) as connection:
        return _upsert_lineage(connection, terms)


def record_class_predictions(database, tool_run_id, predictions):
    initialize(database)
    count = 0
    with transaction(database) as connection:
        if not connection.execute("SELECT 1 FROM annotation_tool_run WHERE tool_run_id=?", (tool_run_id,)).fetchone():
            raise ValueError("class prediction tool run does not exist")
        for prediction in predictions:
            if prediction.subject_kind not in {"spectrum", "alignment_feature"}:
                raise ValueError("unsupported class prediction subject")
            table = prediction.subject_kind
            if not connection.execute(f"SELECT 1 FROM {table} WHERE {table}_id=?", (prediction.subject_id,)).fetchone():
                raise ValueError("class prediction subject does not exist")
            if prediction.rank < 1 or (prediction.formula_rank is not None and prediction.formula_rank < 1):
                raise ValueError("prediction ranks must be positive")
            for posterior in prediction.posteriors.values():
                if not math.isfinite(posterior) or not 0 <= posterior <= 1:
                    raise ValueError("posterior must be between zero and one")
            terms = list(prediction.lineage)
            if not terms or set(prediction.posteriors) - {t.level for t in terms}:
                raise ValueError("prediction posteriors must refer to its nonempty lineage")
            ids = _upsert_lineage(connection, terms)
            for index, (term, term_id) in enumerate(zip(terms, ids)):
                prediction_id = make_id("class-prediction", tool_run_id, prediction.subject_kind,
                                        prediction.subject_id, term_id, prediction.formula_context or "")
                cursor = connection.execute("""INSERT INTO class_prediction VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(class_prediction_id) DO UPDATE SET posterior=excluded.posterior,
                    rank=excluded.rank, is_most_specific=excluded.is_most_specific,
                    formula_rank=excluded.formula_rank, source_row_ref=excluded.source_row_ref""",
                    (prediction_id, tool_run_id, prediction.subject_kind, prediction.subject_id, term_id,
                     term.level, prediction.posteriors.get(term.level), prediction.rank,
                     int(index == len(terms) - 1), prediction.formula_context,
                     prediction.formula_rank, prediction.source_row_ref))
                count += cursor.rowcount
    return count


def list_class_predictions(database, subject_id):
    with closing(connect(database)) as connection:
        return [dict(row) for row in connection.execute("""SELECT p.*, t.ontology_id,
            t.ontology_version, t.name, t.native_id FROM class_prediction p JOIN class_term t USING(term_id)
            WHERE subject_id=? ORDER BY tool_run_id, rank, is_most_specific, level, term_id""", (subject_id,))]
