"""Why one alignment feature ended up with several candidates, and what that does to its claim.

MS-DIAL keeps up to three threshold-passing results per annotator and the catalog now stores all of
them. A set of several candidates is not one thing, though, and treating it as one would be the same
mistake as publishing only the winner. There are at least four reasons a set holds more than one
entry, and they have opposite consequences for what may be claimed:

  no product-ion evidence at all -- the multiplicity is a precursor-mass-window fact and says nothing
  about structure either way. On the reference run this is 1,071 of 1,919 annotated features.

  a rule-based lipid annotation -- the LBM library's registered spectra are barely used for the
  match; the name and its resolution come from diagnostic-ion logic in MS-DIAL's own source, so the
  answer is already at the resolution the fragments support (PC 16:0_18:1, never PC 16:0/18:1(9Z)).
  Spectral indistinguishability is not the right question to ask of it.

  a mass-window artifact -- the candidates differ in molecular formula and their reference masses are
  further apart than the search's own tolerance, so accurate mass separates them. Not chemistry.

  and only then the two that matter: the library genuinely cannot separate these references, which no
  re-acquisition fixes and which lowers the claim from a structure to a substructure or class; or the
  library can separate them and this measurement did not, which is a limitation of the run.

Anything that cannot be decided is `not_assessed` with the reason, never one of the decided states.
"never tested" and "tested and separable" look identical in SQL and mean opposite things.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Any

from .ambiguity import ClassDefinition
from .identifiers import make_id, short_hash
from .storage import connect, initialize, transaction

SUBJECT_KIND = "alignment_feature"

# What the assessment concluded. Exactly one is stored per set.
STATE_NO_MS2_EVIDENCE = "no_ms2_evidence"
STATE_SINGLE_CANDIDATE = "single_candidate"
STATE_RULE_BASED_LIBRARY = "rule_based_library"
STATE_MASS_SEPARABLE = "mass_separable"
STATE_LIBRARY_INDISTINGUISHABLE = "library_indistinguishable"
STATE_QUERY_NON_DISCRIMINATION = "query_non_discrimination"
STATE_NOT_ASSESSED = "not_assessed"

# What the state does to the honest claim ceiling. Never derived into a Level here: assigning a level
# is a curation decision, and this only reports what the evidence permits.
CEILING_NONE = "no_effect"
CEILING_LOWERS = "lowers_to_substructure_or_class"
CEILING_RUN_LIMITED = "run_limited"
CEILING_UNKNOWN = "unknown"

# MS-DIAL's VS1.0 annotation codes for a lipid resolved by diagnostic-ion rules rather than by
# spectral similarity: 400 acyl position, 410 chains, 420 class. Only CompareMS2LipidomicsScanProperties
# sets the flags behind them, so the tag is evidence from the code path rather than from a name.
RULE_BASED_ANNOTATION_TAGS = frozenset({"400", "410", "420"})

# MS-DIAL names its rule-based lipid database this in every export. A secondary signal, used because
# a lipid whose rules resolved only to "matched" carries tag 430 like any other MS/MS match.
DEFAULT_RULE_BASED_DATABASE_IDS = frozenset({"LbmDB"})


@dataclass(frozen=True)
class AssessmentDefinition:
    """Every rule this assessment applied, so a verdict can be recomputed or superseded later."""

    label: str = "candidate-ambiguity-v1"
    # Below this the accurate mass would not actually separate the candidates, so claiming it would
    # overstate what the reference masses support. The run's own MS1 tolerance is not stored in the
    # catalog, so this is a stated convention rather than a reading of the run.
    mass_separation_mda: float = 1.0
    # Within one set every contender was scored against the same query spectrum by the same annotator
    # on the same convention, so this comparison is legitimate even though cross-row ones are not.
    # Above this spread the product-ion evidence did choose, whatever the candidate list looks like.
    discrimination_margin: float = 0.01
    rule_based_database_ids: tuple[str, ...] = tuple(sorted(DEFAULT_RULE_BASED_DATABASE_IDS))
    allow_skeleton_inference: bool = False

    def as_rules(self) -> dict:
        return asdict(self)

    @property
    def rules_sha256(self) -> str:
        """A digest of the rules, with the label excluded, for the reason ClassDefinition does.

        The label travels beside the digest in every identity and every stored row, so folding it in
        would make the digest answer "which run was this" rather than "which rules were these", and
        the whole point is that a label does not change when a threshold does.
        """
        rules = {key: value for key, value in asdict(self).items() if key != "label"}
        return short_hash(
            json.dumps(rules, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
        )


@dataclass
class AssessmentReport:
    definition_label: str = ""
    subjects: int = 0
    states: dict[str, int] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def _number(value: object) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _contending_subset(candidates: list[dict]) -> tuple[list[dict], str]:
    """The candidates the product-ion evidence actually left in play, and the rule that chose them.

    Rank is a lexicographic order on (evidence tier, score), not a score order: a real MS/MS match
    outranks a higher-scoring precursor-only suggestion, which is the intended precedence. So a set
    holding one match beside two suggestions is not a three-way ambiguity, and reading it as one would
    overstate the ambiguity roughly fourfold on the reference run.
    """
    matched = [item for item in candidates if item.get("is_spectrum_match")]
    if matched:
        return matched, "spectrum_matched"
    compared = [item for item in candidates if item.get("is_spectrum_comparison_performed")]
    if compared:
        return compared, "compared_not_matched"
    return list(candidates), "all_kept_no_comparison"


def _is_rule_based(candidates: list[dict], definition: AssessmentDefinition) -> bool:
    rule_based_ids = {item.casefold() for item in definition.rule_based_database_ids}
    for item in candidates:
        if str(item.get("annotation_tag") or "").strip() in RULE_BASED_ANNOTATION_TAGS:
            return True
        if str(item.get("database_id") or "").strip().casefold() in rule_based_ids:
            return True
    return False


def _mass_spread_mda(contenders: list[dict]) -> float | None:
    masses = [_number(item.get("reference_mz")) for item in contenders]
    known = [value for value in masses if value is not None]
    if len(known) < 2:
        return None
    return (max(known) - min(known)) * 1000.0


def _discrimination_margin(contenders: list[dict]) -> float | None:
    scores = [_number(item.get("weighted_dot_product")) for item in contenders]
    known = [value for value in scores if value is not None]
    if len(known) < 2:
        return None
    return max(known) - min(known)


def assess_candidate_set(
    candidates: list[dict],
    *,
    definition: AssessmentDefinition | None = None,
    class_definition: ClassDefinition | None = None,
    covering_class: dict | None = None,
    library_available: bool = False,
) -> dict[str, Any]:
    """Decide why one alignment feature holds the candidates it holds.

    `covering_class` is the single ambiguity class containing every contender, when one was found.
    `library_available` says whether the library side could be consulted at all; without it the two
    library-dependent verdicts are unreachable and the set is `not_assessed` rather than separable.

    The cascade runs run-side first on purpose. A fact the run establishes on its own -- no product
    ions were compared, the contenders differ by more than the mass tolerance, only one contender
    survives -- is never displaced by the library being unavailable.
    """
    definition = definition or AssessmentDefinition()
    contenders, rule = _contending_subset(candidates)
    margin = _discrimination_margin(contenders)
    spread = _mass_spread_mda(contenders)
    verdict: dict[str, Any] = {
        "candidate_count": len(candidates),
        "contender_count": len(contenders),
        "contender_selection_rule": rule,
        "discrimination_margin": margin,
        "max_reference_mz_spread_mda": spread,
        "ambiguity_class_id": None,
        "linked_contender_count": 0,
        "unlinked_contender_count": len(contenders),
        "not_assessed_reason": None,
    }

    if not any(item.get("is_spectrum_comparison_performed") for item in candidates):
        return {
            **verdict,
            "state": STATE_NO_MS2_EVIDENCE,
            "ceiling_effect": CEILING_UNKNOWN,
            "state_reason": (
                "No product-ion spectrum was compared for any candidate, so the multiplicity is a "
                "precursor-mass-window fact. Neither library indistinguishability nor a limitation of "
                "this run is established, and a structure claim rests on no spectral evidence at all."
            ),
        }

    if _is_rule_based(candidates, definition):
        return {
            **verdict,
            "state": STATE_RULE_BASED_LIBRARY,
            "ceiling_effect": CEILING_NONE,
            "state_reason": (
                "The annotation came from MS-DIAL's rule-based lipid path, where the registered "
                "reference spectra serve as a threshold gate and the name and its resolution come "
                "from diagnostic-ion logic in the source. The result is already reported at the "
                "resolution the fragments support, so spectral indistinguishability between library "
                "entries is not the question that applies to it."
            ),
        }

    if len(contenders) == 1:
        return {
            **verdict,
            "state": STATE_SINGLE_CANDIDATE,
            "ceiling_effect": CEILING_NONE,
            "state_reason": (
                "One candidate survives the product-ion evidence. Note that MS-DIAL filters by "
                "threshold before capping the list, so this is not by itself evidence that the search "
                "distinguished it from everything else in the library."
            ),
        }

    if spread is not None and spread > definition.mass_separation_mda:
        formulas = {
            str(item.get("formula") or "").strip() for item in contenders if item.get("formula")
        }
        if len(formulas) > 1:
            return {
                **verdict,
                "state": STATE_MASS_SEPARABLE,
                "ceiling_effect": CEILING_NONE,
                "state_reason": (
                    f"The contenders differ in molecular formula and their reference masses span "
                    f"{spread:.2f} mDa, more than the {definition.mass_separation_mda:.2f} mDa this "
                    "assessment treats as separable. That is an artifact of the search's precursor "
                    "window, resolvable by accurate mass rather than by any spectrum."
                ),
            }

    if margin is not None and margin > definition.discrimination_margin:
        return {
            **verdict,
            "state": STATE_SINGLE_CANDIDATE,
            "ceiling_effect": CEILING_NONE,
            "state_reason": (
                f"The product-ion evidence separates the contenders: their weighted dot products span "
                f"{margin:.4f}, above the {definition.discrimination_margin:.4f} margin. Every "
                "contender was scored against the same query spectrum by the same annotator on the "
                "same convention, so this comparison is meaningful even though cross-feature ones are "
                "not. The list is ranked, not tied."
            ),
        }

    if not library_available:
        return {
            **verdict,
            "state": STATE_NOT_ASSESSED,
            "not_assessed_reason": "library_unavailable",
            "ceiling_effect": CEILING_UNKNOWN,
            "state_reason": (
                "The product-ion evidence did not separate the contenders, and whether the library "
                "can is untested: no reference library covering these records is ingested, so no "
                "ambiguity class exists to consult. This is not the same as the library being able "
                "to separate them, and must not be read as one."
            ),
        }

    if covering_class is None:
        return {
            **verdict,
            "state": STATE_QUERY_NON_DISCRIMINATION,
            "ceiling_effect": CEILING_RUN_LIMITED,
            "state_reason": (
                "The library distinguishes these reference records under the stated rule set and this "
                "measurement did not, so the limitation belongs to this run and its search settings "
                "rather than to chemistry. A better spectrum could resolve it."
            ),
        }

    members = list(covering_class.get("members") or [])
    discriminating = list(covering_class.get("discriminating_mz") or [])
    return {
        **verdict,
        "state": STATE_LIBRARY_INDISTINGUISHABLE,
        "ceiling_effect": CEILING_LOWERS,
        "ambiguity_class_id": covering_class.get("ambiguity_class_id"),
        "linked_contender_count": len(contenders),
        "unlinked_contender_count": 0,
        "state_reason": (
            "One ambiguity class contains every contender, so under "
            f"{(class_definition or ClassDefinition()).definition_id} these reference records were not "
            f"separable from one another. "
            + (
                f"The class carries {len(discriminating)} discriminating product ions, so this is a "
                "threshold effect that a higher-resolution or differently fragmented spectrum could "
                "still break."
                if discriminating
                else "The class carries no discriminating product ion, which is the strongest "
                "statement the stored library data supports: no product ion in it separates these "
                "records."
            )
            + f" The class holds {len(members)} members, and the library scope it was computed over "
            "is a systematically incomplete lower bound."
        ),
    }


_CANDIDATE_COLUMNS = """candidate_rank, is_spectrum_match, is_spectrum_comparison_performed,
    annotation_tag, database_id, formula, inchikey, reference_mz, weighted_dot_product, name"""


def _covering_class(connection, contenders: list[dict], definition_id: str) -> dict | None:
    """The one ambiguity class of this definition that contains every contender, if there is one.

    Every contender must resolve to exactly one reference record, and all of them must be members of
    the same class. A partial overlap is not a covering class: the contender that did not resolve
    could be the one that changes the answer, in either direction.
    """
    keys = [str(item.get("inchikey") or "").strip() for item in contenders]
    if not all(keys):
        return None
    resolved: list[str] = []
    for key in keys:
        rows = connection.execute(
            "SELECT reference_spectrum_id FROM reference_spectrum WHERE inchikey = ?", (key,)
        ).fetchall()
        if len(rows) != 1:
            return None
        resolved.append(rows[0]["reference_spectrum_id"])

    placeholders = ", ".join("?" for _ in resolved)
    rows = connection.execute(
        f"""SELECT c.ambiguity_class_id, c.member_count, c.discriminating_mz_json
            FROM ambiguity_class c
            JOIN ambiguity_class_member m USING(ambiguity_class_id)
            WHERE c.definition_id = ? AND m.reference_spectrum_id IN ({placeholders})
            GROUP BY c.ambiguity_class_id
            HAVING COUNT(DISTINCT m.reference_spectrum_id) = ?
            ORDER BY c.member_count, c.ambiguity_class_id""",
        (definition_id, *resolved, len(set(resolved))),
    ).fetchall()
    if not rows:
        return None
    # The tightest true statement: the smallest class that still contains all of them.
    row = rows[0]
    members = connection.execute(
        "SELECT reference_spectrum_id, record_name, inchikey FROM ambiguity_class_member "
        "WHERE ambiguity_class_id = ?",
        (row["ambiguity_class_id"],),
    ).fetchall()
    try:
        discriminating = json.loads(row["discriminating_mz_json"] or "[]")
    except ValueError:
        discriminating = []
    return {
        "ambiguity_class_id": row["ambiguity_class_id"],
        "members": [dict(item) for item in members],
        "discriminating_mz": discriminating,
    }


def assess_run_candidates(
    database,
    run_id: str,
    *,
    definition: AssessmentDefinition | None = None,
    class_definition: ClassDefinition | None = None,
) -> AssessmentReport:
    """Assess every candidate set of one run and store the verdicts."""
    definition = definition or AssessmentDefinition()
    class_definition = class_definition or ClassDefinition()
    report = AssessmentReport(definition_label=definition.label)
    initialize(database)
    with transaction(database) as connection:
        exists = connection.execute(
            "SELECT 1 FROM analysis_run WHERE run_id = ?", (run_id,)
        ).fetchone()
        if exists is None:
            report.errors.append(f"Run {run_id} does not exist in this catalog.")
            return report
        library_available = bool(
            connection.execute("SELECT 1 FROM reference_spectrum LIMIT 1").fetchone()
        )
        if not library_available:
            report.warnings.append(
                "No reference library is ingested, so every set that reaches the library stage is "
                "recorded as not_assessed rather than as separable."
            )
        subjects = [
            row["subject_id"]
            for row in connection.execute(
                "SELECT DISTINCT subject_id FROM msdial_annotation_candidate "
                "WHERE run_id = ? AND subject_kind = ? ORDER BY subject_id",
                (run_id, SUBJECT_KIND),
            )
        ]
        for subject_id in subjects:
            candidates = [
                dict(row)
                for row in connection.execute(
                    f"SELECT {_CANDIDATE_COLUMNS} FROM msdial_annotation_candidate "
                    "WHERE run_id = ? AND subject_kind = ? AND subject_id = ? "
                    "ORDER BY candidate_rank",
                    (run_id, SUBJECT_KIND, subject_id),
                )
            ]
            if not candidates:
                continue
            covering = None
            if library_available:
                contenders, _ = _contending_subset(candidates)
                covering = _covering_class(connection, contenders, class_definition.definition_id)
            verdict = assess_candidate_set(
                candidates,
                definition=definition,
                class_definition=class_definition,
                covering_class=covering,
                library_available=library_available,
            )
            connection.execute(
                """INSERT INTO candidate_set_assessment(
                       candidate_set_assessment_id, run_id, subject_kind, subject_id,
                       definition_label, rules_sha256, state, not_assessed_reason, ceiling_effect,
                       candidate_count, contender_count, contender_selection_rule,
                       discrimination_margin, max_reference_mz_spread_mda, ambiguity_class_id,
                       linked_contender_count, unlinked_contender_count, state_reason)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                   ON CONFLICT(run_id, subject_kind, subject_id, definition_label, rules_sha256)
                   DO UPDATE SET
                       state = excluded.state,
                       not_assessed_reason = excluded.not_assessed_reason,
                       ceiling_effect = excluded.ceiling_effect,
                       candidate_count = excluded.candidate_count,
                       contender_count = excluded.contender_count,
                       contender_selection_rule = excluded.contender_selection_rule,
                       discrimination_margin = excluded.discrimination_margin,
                       max_reference_mz_spread_mda = excluded.max_reference_mz_spread_mda,
                       ambiguity_class_id = excluded.ambiguity_class_id,
                       linked_contender_count = excluded.linked_contender_count,
                       unlinked_contender_count = excluded.unlinked_contender_count,
                       state_reason = excluded.state_reason""",
                (
                    make_id(
                        "candidate-assessment",
                        run_id,
                        subject_id,
                        definition.label,
                        definition.rules_sha256,
                    ),
                    run_id, SUBJECT_KIND, subject_id, definition.label, definition.rules_sha256,
                    verdict["state"], verdict["not_assessed_reason"], verdict["ceiling_effect"],
                    verdict["candidate_count"], verdict["contender_count"],
                    verdict["contender_selection_rule"], verdict["discrimination_margin"],
                    verdict["max_reference_mz_spread_mda"], verdict["ambiguity_class_id"],
                    verdict["linked_contender_count"], verdict["unlinked_contender_count"],
                    verdict["state_reason"],
                ),
            )
            report.subjects += 1
            report.states[verdict["state"]] = report.states.get(verdict["state"], 0) + 1
    return report
