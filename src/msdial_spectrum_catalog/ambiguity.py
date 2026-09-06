"""Library-side ambiguity classes: which other reference entries a match cannot be told apart from.

Given that a query matched reference entry A, the ambiguity class of A is the set of entries that would
be indistinguishable from A under the study's measurement conditions. It is a property of the library and
is computable without seeing any query, which matters because MS-DIAL applies its score threshold before
capping a candidate list at three hits per annotator: an entry that is near-identical to A but narrowly
fails the threshold never appears in the per-query candidate list at all.

Four properties of this computation are deliberate and must not be "simplified" away.

ANCHORED, NOT PARTITIONED. N(A) = {A} union {X : sim(A, X) >= threshold}. The neighbourhood is defined
relative to A, so there is no seed and the result does not depend on iteration order. Greedy clustering
would be seed-dependent. Similarity is not transitive, so no transitive closure is taken; whether N(A)
happens to be a clique is recorded instead of being forced.

THREE OUTCOMES, NOT TWO. A pair can be indistinguishable, distinguishable, or not assessable. A pair that
fails the admissibility gate is `insufficient_evidence` and a pair measured under different conditions is
`condition_mismatch`; neither asserts indistinguishability and neither may be merged into a class.

THE FORMULA GUARD SPLITS, IT DOES NOT REPORT. Two entries inside the precursor window whose formulas
differ are separable by formula and mass evidence alone, so they are an artifact of the window width and
are tagged `isobaric_not_isomeric` rather than reported as ambiguity.

THE THRESHOLD IS A CONVENTION. There is no ground-truth set of isomer pairs known to be separable by
retention time, mobility or an authentic standard, so no error rate can be attached to any threshold. The
full rule set therefore travels with every class it produced, and so does the library scope: a class
computed on public spectra alone is systematically incomplete and says so.
"""

from __future__ import annotations

import json
import math
import zlib
from dataclasses import asdict, dataclass, field
from typing import Iterator

from .identifiers import make_id, short_hash
from .similarity import Peaks, compare, normalize_to_base_peak
from .storage import connect, initialize, transaction

SUBJECT_KIND = "reference_spectrum"
METHOD = "ambiguity_weighted_cosine"
SECONDARY_METHOD = "entropy_similarity"
SCORE_CONVENTION = "cosine"

# Collision energy is binned rather than compared exactly: parseable values in the lab's public release
# mix '20', '20.0' and '45HCD', and a 1 eV difference is not a different measurement condition.
CE_BIN_WIDTH = 10.0

# A block larger than this is reported. The blocking step is what keeps the computation sparse, so a
# pathological block -- a precursor window holding thousands of records -- must be visible rather than
# silently quadratic.
LARGE_BLOCK = 500

# Below these spreads the evidence would not actually break the tie, so claiming it would overstate what
# the rows support.
RT_SEPARATION_MIN = 0.05
CCS_SEPARATION = 1.0

_UNKNOWN = "unknown"


@dataclass(frozen=True)
class ClassDefinition:
    """Every parameter of one ambiguity computation, versioned so a class can be recomputed later."""

    definition_id: str = "ambiguity-v1"
    weighted_cosine_threshold: float = 0.90
    entropy_similarity_threshold: float = 0.85
    mz_tolerance_da: float = 0.025
    minimum_informative_peaks: int = 6
    minimum_matched_peaks: int = 4
    relative_floor: float = 0.01
    precursor_tolerance_da: float = 0.01
    precursor_tolerance_ppm: float = 10.0
    require_formula_agreement: bool = True
    require_condition_match: bool = True
    linkage: str = "anchored_neighbourhood"

    def as_rules(self) -> dict:
        """Return every parameter as a plain dict, so `ClassDefinition(**rules)` reconstructs it."""
        return asdict(self)

    @property
    def rules_sha256(self) -> str:
        """A digest of every threshold and rule, not of the free-text label they ran under.

        A class asserts that its members could not be told apart *under a stated rule set*. The
        identifier used to encode only definition_id, which defaults to "ambiguity-v1" and is a label
        the caller may leave untouched while changing --weighted-cosine, --entropy or --tolerance on
        the command line. Two runs at different thresholds therefore produced the same identifiers and
        the upsert replaced the earlier rows in place: same id, different meaning, and nothing stored
        that said so. Deriving the identifier from the rules makes a recomputation at a different
        threshold a different class rather than a silent overwrite of the old one.

        definition_id is excluded from the digest because it travels beside it in every identity and
        every stored row. Folding it in would make the digest answer "which run was this" rather than
        "which rules were these", and the point is that the label does not change when a threshold
        does.
        """
        rules = {key: value for key, value in asdict(self).items() if key != "definition_id"}
        return short_hash(_canonical(rules))


@dataclass
class AmbiguityReport:
    definition_id: str = ""
    blocks: int = 0
    pairs_compared: int = 0
    pairs_insufficient_evidence: int = 0
    pairs_condition_mismatch: int = 0
    pairs_isobaric_not_isomeric: int = 0
    edges: int = 0
    classes: int = 0
    singletons: int = 0
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def _canonical(payload: object) -> str:
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _row_to_dict(row) -> dict:
    return {key: row[key] for key in row.keys()}


def _text(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _mode_label(value: object | None) -> str:
    return (_text(value) or _UNKNOWN).lower()


def _skeleton(row: dict) -> str | None:
    """Prefer the stored first block, falling back to the InChIKey's own first block."""
    stored = _text(row.get("inchikey_skeleton"))
    if stored:
        return stored
    inchikey = _text(row.get("inchikey"))
    return inchikey.split("-")[0] if inchikey else None


def instrument_class_label(value: object | None) -> str:
    """Normalize a stored instrument class for condition comparison."""
    return (_text(value) or _UNKNOWN).lower()


def collision_energy_bin(value: float | None) -> str:
    """Bin a numeric collision energy, or report that the condition is unestablished."""
    if value is None:
        return _UNKNOWN
    index = math.floor(float(value) / CE_BIN_WIDTH)
    return f"{index * CE_BIN_WIDTH:g}-{(index + 1) * CE_BIN_WIDTH:g}"


def condition_key(row: dict) -> tuple[str, str]:
    """The measurement condition of one reference record: instrument class and collision-energy bin."""
    return (
        instrument_class_label(row.get("instrument_class")),
        collision_energy_bin(row.get("collision_energy_value")),
    )


def blocking_key(row: dict, *, definition: ClassDefinition | None = None) -> str:
    """Label the comparison block a record belongs to: ion mode, precursor type and m/z window."""
    definition = definition or ClassDefinition()
    mode = _mode_label(row.get("ion_mode"))
    precursor_type = _text(row.get("precursor_type")) or _UNKNOWN
    mz = row.get("precursor_mz")
    window = _UNKNOWN if mz is None else str(math.floor(float(mz) / definition.precursor_tolerance_da))
    return f"{mode}|{precursor_type}|{window}"


def _precursor_tolerance(mz: float, definition: ClassDefinition) -> float:
    return max(definition.precursor_tolerance_da, abs(mz) * definition.precursor_tolerance_ppm / 1e6)


def _within_precursor_window(row_a: dict, row_b: dict, definition: ClassDefinition) -> bool:
    mz_a, mz_b = row_a.get("precursor_mz"), row_b.get("precursor_mz")
    if mz_a is None or mz_b is None:
        return False
    tolerance = max(_precursor_tolerance(mz_a, definition), _precursor_tolerance(mz_b, definition))
    return abs(float(mz_a) - float(mz_b)) <= tolerance


def _peaks_from_payload(payload: object) -> Peaks:
    raw = payload.get("peaks", []) if isinstance(payload, dict) else payload
    peaks: Peaks = []
    for entry in raw or []:
        try:
            mz, intensity = float(entry[0]), float(entry[1])
        except (TypeError, ValueError, IndexError):
            continue
        peaks.append((mz, intensity))
    peaks.sort()
    return peaks


def _load_peaks(connection, payload_sha256: str) -> Peaks:
    row = connection.execute(
        "SELECT compression, payload FROM spectrum_blob WHERE payload_sha256 = ?", (payload_sha256,)
    ).fetchone()
    if row is None:
        return []
    payload = bytes(row["payload"])
    if row["compression"] == "zlib-json":
        payload = zlib.decompress(payload)
    return _peaks_from_payload(json.loads(payload.decode("utf-8")))


_ROW_COLUMNS = """reference_spectrum_id, library_id, library_record_index, record_name, inchikey,
    inchikey_skeleton, smiles, formula, ontology, precursor_mz, precursor_type, ion_mode,
    instrument_type, instrument_class, collision_energy_raw, collision_energy_value, rt_min, ccs,
    peak_count, payload_sha256"""


def iter_blocks(
    database, definition: ClassDefinition, *, library_ids: list[str] | None = None
) -> Iterator[list[dict]]:
    """Yield the comparison blocks of a reference library, one anchored precursor window at a time.

    A block is the anchor plus every record sharing its ion mode and precursor type whose precursor m/z
    lies within max(precursor_tolerance_da, precursor_tolerance_ppm) of it. Anything outside that window
    is separable by MS1 mass alone, so it is never compared. The stream is walked with a sliding window
    and blocks with identical membership are yielded once, so nothing resembling a dense N x N matrix is
    ever built: at the scale of the lab's full positive library a dense double matrix would be about
    46 TB, while the admitted edges are sparse enough to hold as a list.
    """
    if library_ids is not None and not library_ids:
        return
    clause = ""
    parameters: list[str] = []
    if library_ids is not None:
        clause = " AND library_id IN (" + ", ".join("?" for _ in library_ids) + ")"
        parameters = list(library_ids)
    connection = connect(database)
    try:
        cursor = connection.cursor()
        cursor.execute(
            f"""SELECT {_ROW_COLUMNS} FROM reference_spectrum
                WHERE precursor_mz IS NOT NULL{clause}
                ORDER BY lower(COALESCE(ion_mode, '')), COALESCE(precursor_type, ''),
                         precursor_mz, reference_spectrum_id""",
            parameters,
        )
        seen: set[frozenset[str]] = set()
        group: tuple[str, str] | None = None
        window: list[dict] = []
        anchor = 0
        for raw in cursor:
            row = _row_to_dict(raw)
            key = (_mode_label(row.get("ion_mode")), _text(row.get("precursor_type")) or _UNKNOWN)
            if key != group:
                yield from _drain(connection, window, anchor, definition, seen)
                # Two blocks can only share a membership when their windows overlap, which cannot happen
                # across ion modes or precursor types, so the dedup set need not outlive the group.
                group, window, anchor, seen = key, [], 0, set()
            window.append(row)
            # Once the newest row lies beyond the anchor's window, every right-hand neighbour of that
            # anchor has been read, so its block is complete and the window can move on.
            while anchor < len(window):
                centre = float(window[anchor]["precursor_mz"])
                if float(row["precursor_mz"]) - centre <= _precursor_tolerance(centre, definition):
                    break
                block = _block_at(connection, window, anchor, definition, seen)
                if block is not None:
                    yield block
                anchor = _trim(window, anchor + 1, definition)
        yield from _drain(connection, window, anchor, definition, seen)
    finally:
        connection.close()


def _drain(
    connection,
    window: list[dict],
    anchor: int,
    definition: ClassDefinition,
    seen: set[frozenset[str]],
) -> Iterator[list[dict]]:
    """Emit the blocks of the anchors still pending when a group or the stream ends."""
    while anchor < len(window):
        block = _block_at(connection, window, anchor, definition, seen)
        if block is not None:
            yield block
        anchor += 1


def _trim(window: list[dict], anchor: int, definition: ClassDefinition) -> int:
    """Drop rows that have fallen behind the current anchor's window, keeping the window bounded."""
    if anchor >= len(window):
        return anchor
    centre = float(window[anchor]["precursor_mz"])
    lower = centre - _precursor_tolerance(centre, definition)
    while window and float(window[0]["precursor_mz"]) < lower:
        window.pop(0)
        anchor -= 1
    return anchor


def _block_at(
    connection,
    window: list[dict],
    anchor: int,
    definition: ClassDefinition,
    seen: set[frozenset[str]],
) -> list[dict] | None:
    """Collect the anchor's precursor window, loading peaks only for blocks that can produce a pair."""
    centre = float(window[anchor]["precursor_mz"])
    span = _precursor_tolerance(centre, definition)
    block = [row for row in window if abs(float(row["precursor_mz"]) - centre) <= span]
    if len(block) < 2:
        return None
    members = frozenset(row["reference_spectrum_id"] for row in block)
    if members in seen:
        return None
    seen.add(members)
    for row in block:
        if "peaks" not in row:
            row["peaks"] = _load_peaks(connection, row["payload_sha256"])
    return block


def is_clique(members, edges) -> bool:
    """Report whether every pair of members is mutually above threshold.

    Similarity is not transitive, so an anchored neighbourhood need not be a clique: A and C can both be
    indistinguishable from B while being distinguishable from each other. A non-clique neighbourhood is
    still an honest statement about the anchor, only a weaker one -- "A or X or Y, though X and Y are
    mutually distinguishable" -- so which it is gets recorded rather than resolved by forcing a partition.
    `edges` is any iterable of id pairs, which includes a dict keyed by pairs.
    """
    present = {tuple(sorted((str(a), str(b)))) for a, b in edges}
    ids = sorted({str(member) for member in members})
    for index, first in enumerate(ids):
        for second in ids[index + 1:]:
            if (first, second) not in present:
                return False
    return True


def discriminating_ions(members: list[dict], *, tolerance: float, relative_floor: float = 0.01) -> list[dict]:
    """Report the product ions that would separate the members of a class.

    An ion carried by some members and absent from others is what turns "A or B" into "A or B, separable
    by m/z X", which is actionable; an ion shared by all members separates nothing. Ordered by how evenly
    the ion splits the class, because a 50/50 split is the most informative single measurement.
    """
    member_ids = [str(member["reference_spectrum_id"]) for member in members]
    if len(member_ids) < 2:
        return []
    entries: list[tuple[float, float, str]] = []
    for member in members:
        identifier = str(member["reference_spectrum_id"])
        for mz, intensity in normalize_to_base_peak(member.get("peaks") or []):
            if intensity >= relative_floor:
                entries.append((mz, intensity, identifier))
    entries.sort()
    groups: list[list[tuple[float, float, str]]] = []
    for entry in entries:
        if groups and entry[0] - groups[-1][0][0] <= tolerance:
            groups[-1].append(entry)
        else:
            groups.append([entry])
    total = len(member_ids)
    result: list[dict] = []
    for group in groups:
        carriers: dict[str, float] = {}
        for mz, intensity, identifier in group:
            carriers[identifier] = max(carriers.get(identifier, 0.0), intensity)
        if len(carriers) == total:
            continue
        present = sorted(carriers)
        result.append(
            {
                "mz": round(sum(mz for mz, _, _ in group) / len(group), 5),
                "present_in": present,
                "absent_from": sorted(identifier for identifier in member_ids if identifier not in carriers),
                "mean_relative_intensity": round(sum(carriers.values()) / len(carriers), 5),
                "split_fraction": round(len(present) / total, 5),
            }
        )
    result.sort(
        key=lambda ion: (
            abs(ion["split_fraction"] - 0.5),
            -ion["mean_relative_intensity"],
            ion["mz"],
        )
    )
    return result


def discriminating_evidence_needed(members: list[dict]) -> str:
    """Name the evidence that would break the tie, derived from what the member rows actually carry."""
    needed: list[str] = []
    if _separates(members, "rt_min", RT_SEPARATION_MIN):
        needed.append("RT")
    if _separates(members, "ccs", CCS_SEPARATION):
        needed.append("IM")
    # An authentic standard resolves any class, so RS is always available even when the rows carry
    # neither a retention time nor a collision cross-section.
    needed.append("RS")
    return ",".join(needed)


def _separates(members: list[dict], column: str, minimum_spread: float) -> bool:
    values = [float(member[column]) for member in members if member.get(column) is not None]
    return len(values) >= 2 and max(values) - min(values) >= minimum_spread


def _agreement(values: list[str | None], same: str, different: str, unknown: str) -> str:
    if any(value is None for value in values):
        # Agreement cannot be established from a missing field, and reporting "different" would be a
        # claim the rows do not support.
        return unknown
    distinct = {value for value in values}
    if len(distinct) == 1:
        return same
    if len(distinct) == len(values):
        return different
    return "mixed"


def _library_scope(connection, library_ids: list[str] | None) -> dict:
    clause = ""
    parameters: list[str] = []
    if library_ids is not None:
        clause = " WHERE library_id IN (" + ", ".join("?" for _ in library_ids) + ")" if library_ids else " WHERE 0"
        parameters = list(library_ids)
    libraries = [
        {
            "library_id": row["library_id"],
            "library_name": row["library_name"],
            "library_version": row["library_version"],
            "library_kind": row["library_kind"],
            "record_count": row["record_count"],
        }
        for row in connection.execute(
            "SELECT library_id, library_name, library_version, library_kind, record_count "
            f"FROM reference_library{clause} ORDER BY library_id",
            parameters,
        )
    ]
    return {
        "compared_libraries": libraries,
        "completeness": "systematically_incomplete",
        "note": (
            "Only the listed libraries were compared. An entry absent from them cannot be excluded, "
            "so the class is a lower bound on the true ambiguity."
        ),
    }


def compute_ambiguity_classes(
    database,
    *,
    definition: ClassDefinition | None = None,
    library_ids: list[str] | None = None,
    tool_run_id: str | None = None,
    progress=None,
) -> AmbiguityReport:
    """Compute and persist ambiguity classes over the reference library (blocking through evidence)."""
    definition = definition or ClassDefinition()
    initialize(database)
    report = AmbiguityReport(definition_id=definition.definition_id)

    edges: dict[tuple[str, str], dict] = {}
    adjacency: dict[str, set[str]] = {}
    rows: dict[str, dict] = {}
    compared: set[tuple[str, str]] = set()
    compared_ids: set[str] = set()
    missing_formula = 0
    largest_block = 0

    for block in iter_blocks(database, definition, library_ids=library_ids):
        report.blocks += 1
        largest_block = max(largest_block, len(block))
        if progress is not None:
            progress(
                {
                    "blocks": report.blocks,
                    "block_size": len(block),
                    "blocking_key": blocking_key(block[0], definition=definition),
                    "pairs_compared": report.pairs_compared,
                    "edges": len(edges),
                }
            )
        for index, row_a in enumerate(block):
            for row_b in block[index + 1:]:
                id_a, id_b = row_a["reference_spectrum_id"], row_b["reference_spectrum_id"]
                key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
                if key in compared:
                    continue
                if not _within_precursor_window(row_a, row_b, definition):
                    continue
                compared.add(key)
                compared_ids.update(key)
                report.pairs_compared += 1

                if definition.require_condition_match and condition_key(row_a) != condition_key(row_b):
                    # A different instrument class or collision energy is a third outcome: the pair was
                    # not shown to be distinguishable, and indistinguishability cannot be asserted either.
                    report.pairs_condition_mismatch += 1
                    continue

                formula_a, formula_b = _text(row_a.get("formula")), _text(row_b.get("formula"))
                if formula_a is None or formula_b is None:
                    missing_formula += 1
                    if definition.require_formula_agreement:
                        report.pairs_isobaric_not_isomeric += 1
                        continue
                elif formula_a != formula_b:
                    report.pairs_isobaric_not_isomeric += 1
                    if definition.require_formula_agreement:
                        continue

                comparison = compare(
                    row_a["peaks"],
                    row_b["peaks"],
                    tolerance=definition.mz_tolerance_da,
                    minimum_informative_peaks=definition.minimum_informative_peaks,
                    minimum_matched_peaks=definition.minimum_matched_peaks,
                    relative_floor=definition.relative_floor,
                )
                if not comparison.comparable or comparison.entropy_similarity is None:
                    report.pairs_insufficient_evidence += 1
                    continue
                if (
                    comparison.weighted_cosine < definition.weighted_cosine_threshold
                    or comparison.entropy_similarity < definition.entropy_similarity_threshold
                ):
                    continue

                edges[key] = {
                    "score": comparison.weighted_cosine,
                    "secondary_score": comparison.entropy_similarity,
                    "matched_peak_count": comparison.matched_peak_count,
                }
                adjacency.setdefault(id_a, set()).add(id_b)
                adjacency.setdefault(id_b, set()).add(id_a)
                rows[id_a], rows[id_b] = row_a, row_b

    report.edges = len(edges)
    if largest_block >= LARGE_BLOCK:
        report.warnings.append(
            f"The largest precursor block held {largest_block} records; comparison inside a block is "
            "quadratic, so a block this size is worth inspecting"
        )
    if missing_formula:
        report.warnings.append(
            f"{missing_formula} pairs were refused because a FORMULA is missing, which is a different "
            "situation from a formula that disagrees"
        )

    with transaction(database) as connection:
        scoped = ""
        if library_ids is not None:
            scoped = (
                " AND library_id IN (" + ", ".join("?" for _ in library_ids) + ")"
                if library_ids
                else " AND 0"
            )
        total = connection.execute(
            f"SELECT COUNT(*) FROM reference_spectrum WHERE precursor_mz IS NOT NULL{scoped}",
            list(library_ids or []),
        ).fetchone()[0]
        report.singletons = total - len(adjacency)
        scope = _library_scope(connection, library_ids)
        # The label, plus a digest of the rules it stood for. The label alone was the stored version,
        # and it defaults to "ambiguity-v1" whatever thresholds the caller passed on the command line.
        method_version = f"{definition.definition_id}+{definition.rules_sha256}"

        for (id_a, id_b), edge in sorted(edges.items()):
            connection.execute(
                """INSERT INTO spectrum_similarity(
                       similarity_id, subject_kind_a, subject_id_a, subject_kind_b, subject_id_b,
                       method, method_version, score, score_convention, secondary_method,
                       secondary_score, matched_peak_count, mz_tolerance_da, tool_run_id)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                   ON CONFLICT(subject_kind_a, subject_id_a, subject_kind_b, subject_id_b, method)
                   DO UPDATE SET
                       method_version = excluded.method_version,
                       score = excluded.score,
                       score_convention = excluded.score_convention,
                       secondary_method = excluded.secondary_method,
                       secondary_score = excluded.secondary_score,
                       matched_peak_count = excluded.matched_peak_count,
                       mz_tolerance_da = excluded.mz_tolerance_da,
                       tool_run_id = excluded.tool_run_id""",
                (
                    make_id(
                        "similarity", METHOD, definition.rules_sha256, short_hash(f"{id_a}|{id_b}")
                    ),
                    SUBJECT_KIND, id_a, SUBJECT_KIND, id_b,
                    METHOD, method_version, edge["score"], SCORE_CONVENTION,
                    SECONDARY_METHOD, edge["secondary_score"], edge["matched_peak_count"],
                    definition.mz_tolerance_da, tool_run_id,
                ),
            )
        # An edge this run refused must not survive as if this run had admitted it, so stale edges
        # between records it actually compared are dropped -- but only within its own rule set. Rows
        # produced under a different threshold are a different answer to a different question, and now
        # that the identity encodes the rules they coexist rather than replacing one another.
        for stale_a, stale_b in connection.execute(
            "SELECT subject_id_a, subject_id_b FROM spectrum_similarity "
            "WHERE method = ? AND method_version = ?",
            (METHOD, method_version),
        ).fetchall():
            if (stale_a, stale_b) not in edges and stale_a in compared_ids and stale_b in compared_ids:
                connection.execute(
                    "DELETE FROM spectrum_similarity WHERE subject_kind_a = ? AND subject_id_a = ? "
                    "AND subject_kind_b = ? AND subject_id_b = ? AND method = ? AND method_version = ?",
                    (SUBJECT_KIND, stale_a, SUBJECT_KIND, stale_b, METHOD, method_version),
                )

        unestablished = 0
        for anchor_id in sorted(compared_ids):
            class_id = _class_id(definition, anchor_id)
            connection.execute(
                "DELETE FROM ambiguity_class_member WHERE ambiguity_class_id = ?", (class_id,)
            )
            neighbours = adjacency.get(anchor_id)
            if not neighbours:
                connection.execute(
                    "DELETE FROM ambiguity_class WHERE ambiguity_class_id = ?", (class_id,)
                )
                continue
            members = [rows[anchor_id]] + [rows[member] for member in sorted(neighbours)]
            record = _class_record(definition, anchor_id, members, edges, scope)
            if not record["condition_established"]:
                unestablished += 1
            connection.execute(
                """INSERT INTO ambiguity_class(
                       ambiguity_class_id, definition_id, library_scope_json, blocking_key,
                       condition_scope_json, member_count, formula_agreement, skeleton_agreement,
                       linkage_rule, min_pairwise_score, score_convention, discriminating_mz_json,
                       discriminating_evidence_needed, tool_run_id)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                   ON CONFLICT(ambiguity_class_id) DO UPDATE SET
                       definition_id = excluded.definition_id,
                       library_scope_json = excluded.library_scope_json,
                       blocking_key = excluded.blocking_key,
                       condition_scope_json = excluded.condition_scope_json,
                       member_count = excluded.member_count,
                       formula_agreement = excluded.formula_agreement,
                       skeleton_agreement = excluded.skeleton_agreement,
                       linkage_rule = excluded.linkage_rule,
                       min_pairwise_score = excluded.min_pairwise_score,
                       score_convention = excluded.score_convention,
                       discriminating_mz_json = excluded.discriminating_mz_json,
                       discriminating_evidence_needed = excluded.discriminating_evidence_needed,
                       tool_run_id = excluded.tool_run_id""",
                (
                    class_id, definition.definition_id, _canonical(record["library_scope"]),
                    record["blocking_key"], _canonical(record["condition_scope"]), len(members),
                    record["formula_agreement"], record["skeleton_agreement"], definition.linkage,
                    record["min_pairwise_score"], SCORE_CONVENTION,
                    _canonical(record["discriminating_mz"]), record["discriminating_evidence_needed"],
                    tool_run_id,
                ),
            )
            for ordinal, member in enumerate(members, start=1):
                connection.execute(
                    """INSERT INTO ambiguity_class_member(
                           ambiguity_class_member_id, ambiguity_class_id, reference_spectrum_id,
                           inchikey, inchikey_skeleton, smiles, formula, record_name)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        make_id("ambiguity-member", class_id, ordinal),
                        class_id,
                        member["reference_spectrum_id"],
                        _text(member.get("inchikey")),
                        _skeleton(member),
                        _text(member.get("smiles")),
                        _text(member.get("formula")),
                        _text(member.get("record_name")),
                    ),
                )
            report.classes += 1
        if unestablished:
            report.warnings.append(
                f"{unestablished} classes hold under an unestablished collision-energy bin, so the "
                "condition they are asserted under is weaker than the rest"
            )
    return report


def _class_id(definition: ClassDefinition, anchor_id: str) -> str:
    # The rules digest is part of the identity, not decoration. See ClassDefinition.rules_sha256.
    return make_id(
        "ambiguity-class",
        definition.definition_id,
        definition.rules_sha256,
        short_hash(anchor_id),
    )


def _class_record(
    definition: ClassDefinition,
    anchor_id: str,
    members: list[dict],
    edges: dict[tuple[str, str], dict],
    scope: dict,
) -> dict:
    member_ids = [member["reference_spectrum_id"] for member in members]
    inside = set(member_ids)
    internal = {key: edge for key, edge in edges.items() if key[0] in inside and key[1] in inside}
    instrument_classes = sorted({instrument_class_label(member.get("instrument_class")) for member in members})
    ce_bins = sorted({collision_energy_bin(member.get("collision_energy_value")) for member in members})
    condition_established = _UNKNOWN not in ce_bins and _UNKNOWN not in instrument_classes
    return {
        "blocking_key": blocking_key(members[0], definition=definition),
        "formula_agreement": _agreement(
            [_text(member.get("formula")) for member in members],
            "same_formula",
            "isobaric_not_isomeric",
            "unknown_formula",
        ),
        "skeleton_agreement": _agreement(
            [_skeleton(member) for member in members],
            "same_skeleton",
            "different_skeleton",
            "unknown_skeleton",
        ),
        # The honest minimum over what was actually established: for a clique that is every member pair,
        # and for a non-clique it is the anchor's edges plus whichever member pairs also passed.
        "min_pairwise_score": min(edge["score"] for edge in internal.values()) if internal else None,
        "discriminating_mz": discriminating_ions(
            members, tolerance=definition.mz_tolerance_da, relative_floor=definition.relative_floor
        ),
        "discriminating_evidence_needed": discriminating_evidence_needed(members),
        "condition_established": condition_established,
        "library_scope": {
            **scope,
            "member_libraries": sorted({member["library_id"] for member in members}),
        },
        "condition_scope": {
            "anchor_reference_spectrum_id": anchor_id,
            "clique": is_clique(member_ids, internal),
            "instrument_classes": instrument_classes,
            "collision_energy_bins": ce_bins,
            "collision_energy_bin_width": CE_BIN_WIDTH,
            "condition_established": condition_established,
            "require_condition_match": definition.require_condition_match,
            # The DDL has no column for the rule set, and a threshold that cannot be recovered cannot be
            # recalibrated, so the whole definition travels with the class it produced.
            "definition_rules": definition.as_rules(),
        },
    }


def ambiguity_class_for(database, reference_spectrum_id: str) -> dict | None:
    """Return the ambiguity class anchored on one reference entry, ready to become a candidate list.

    This is the report-time lookup: given the entry a query matched, it answers which other entries the
    match did not distinguish it from, under which condition that holds, which product ions would
    separate them, and what evidence would break the tie. Members come back anchor first so the caller
    can map them straight onto ranked candidates.
    """
    connection = connect(database)
    try:
        rows = connection.execute(
            """SELECT c.* FROM ambiguity_class c
               JOIN ambiguity_class_member m USING(ambiguity_class_id)
               WHERE m.reference_spectrum_id = ?
               ORDER BY c.created_at DESC, c.ambiguity_class_id""",
            (reference_spectrum_id,),
        ).fetchall()
        for row in rows:
            condition_scope = json.loads(row["condition_scope_json"] or "{}")
            if condition_scope.get("anchor_reference_spectrum_id") != reference_spectrum_id:
                continue
            members = [
                _row_to_dict(member)
                for member in connection.execute(
                    "SELECT * FROM ambiguity_class_member WHERE ambiguity_class_id = ? "
                    "ORDER BY ambiguity_class_member_id",
                    (row["ambiguity_class_id"],),
                )
            ]
            ordered = [member for member in members if member["reference_spectrum_id"] == reference_spectrum_id]
            ordered += sorted(
                (member for member in members if member["reference_spectrum_id"] != reference_spectrum_id),
                key=lambda member: member["reference_spectrum_id"],
            )
            result = _row_to_dict(row)
            result["condition_scope"] = condition_scope
            result["library_scope"] = json.loads(row["library_scope_json"] or "{}")
            result["discriminating_mz"] = json.loads(row["discriminating_mz_json"] or "[]")
            result["anchor_reference_spectrum_id"] = reference_spectrum_id
            result["clique"] = bool(condition_scope.get("clique"))
            result["members"] = [
                {
                    "rank": ordinal,
                    "reference_spectrum_id": member["reference_spectrum_id"],
                    "compound_name": member["record_name"],
                    "formula": member["formula"],
                    "inchikey": member["inchikey"],
                    "inchikey_skeleton": member["inchikey_skeleton"],
                    "smiles": member["smiles"],
                    "is_anchor": member["reference_spectrum_id"] == reference_spectrum_id,
                }
                for ordinal, member in enumerate(ordered, start=1)
            ]
            return result
        return None
    finally:
        connection.close()
