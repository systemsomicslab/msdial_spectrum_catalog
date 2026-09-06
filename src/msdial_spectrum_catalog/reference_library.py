"""Reference-library ingest, condition normalization and skeleton consensus spectra.

An ambiguity class asks which OTHER library entries would be indistinguishable from a matched entry
under the study's measurement conditions. Two things have to happen before that question can be posed,
and both are prerequisites rather than refinements.

CONDITION NORMALIZATION. A quarter of the records in the lab's public positive release carry an unusable
COLLISIONENERGY, and the parseable ones mix '20', '20.0' and '45HCD'; instrument strings mix vendor names
with analyzer abbreviations. Comparing spectra across unnormalized conditions compares fragmentation
regimes rather than compounds, so INSTRUMENTTYPE is mapped to an analyzer class and COLLISIONENERGY to a
value, a unit when the string states one, and a bin. A missing collision energy becomes CE_UNKNOWN and
never 0, because a missing energy is not an energy of zero.

SKELETON CONSENSUS. The same release holds roughly fifteen records per unique InChIKey first block.
Without merging replicates, a similarity survey is dominated by repeat measurements of one compound,
which is not the question being asked. Records are therefore grouped by (InChIKey first block, ion mode,
precursor type, instrument class, collision-energy bin) and merged into one consensus spectrum per group.

Nothing here decides whether two entries are indistinguishable; that is the similarity stage's work.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from .identifiers import make_id
from .parsers import MspRecord, parse_number, read_msp, sha256_file
from .similarity import normalize_to_base_peak
from .storage import connect, initialize, transaction

# reference_spectrum has no CHECK constraint, and the distinction is load-bearing: an in-silico predicted
# spectrum must never be reported as spectral-library evidence, so an unrecognised kind is refused.
EXPERIMENTAL_LIBRARY_KIND = "experimental_reference"
IN_SILICO_LIBRARY_KIND = "in_silico_predicted"
LIBRARY_KINDS = (EXPERIMENTAL_LIBRARY_KIND, IN_SILICO_LIBRARY_KIND)

UNKNOWN_INSTRUMENT_CLASS = "UNKNOWN"
INSTRUMENT_CLASS_NAMES = ("FT", "TOF", "IT", "QQQ", UNKNOWN_INSTRUMENT_CLASS)

UNKNOWN_COLLISION_ENERGY_BIN = "CE_UNKNOWN"
COLLISION_ENERGY_BIN_WIDTH = 10.0

CONSENSUS_MZ_BIN_WIDTH = 0.01
CONSENSUS_MINIMUM_MEMBER_FRACTION = 0.5

# Field values that state the absence of a value. 'NA' is deliberately absent: it is also a formula.
_MISSING_TOKENS = frozenset({"", "NAN", "NULL", "N/A", "NONE"})

# The InChIKey first block is a 14-character hash of the connectivity layer.
_SKELETON_LENGTH = 14

# INSTRUMENTTYPE is free text, so the table is scanned in order and the first informative substring wins.
# Order carries the meaning:
#   - 'ITFT' precedes 'IT' because an ion-trap-Orbitrap hybrid detects ions in the FT cell, so its spectra
#     behave like FT data; matching 'IT' first would misclassify every LC-ESI-ITFT record.
#   - every TOF rule precedes 'IT' because separators are squeezed out before matching, and 'ESI-TOF'
#     becomes 'ESITOF', which contains 'IT' across the join.
#   - bare 'FT' is last so a more specific analyzer name is preferred wherever one is stated.
INSTRUMENT_CLASSES: tuple[tuple[str, str], ...] = (
    ("ITFT", "FT"),
    ("QFT", "FT"),
    ("FTICR", "FT"),
    ("FTMS", "FT"),
    ("ORBITRAP", "FT"),
    ("EXACTIVE", "FT"),
    ("EXPLORIS", "FT"),
    ("QTOF", "TOF"),
    ("TOF", "TOF"),
    ("QQQ", "QQQ"),
    ("TRIPLEQUAD", "QQQ"),
    ("QQ", "QQQ"),
    ("IONTRAP", "IT"),
    ("LTQ", "IT"),
    ("IT", "IT"),
    ("FT", "FT"),
)

# COMMENT tokens consulted when COLLISIONENERGY itself is unusable, in order of preference.
COLLISION_ENERGY_COMMENT_KEYS = (
    "COLLISIONENERGY",
    "COLLISION_ENERGY",
    "CE",
    "NCE",
    "HCD",
    "COLLISION",
)

# IONMODE spellings. Grouping on the raw string would split 'Positive' from 'POSITIVE'.
ION_MODES: dict[str, str] = {
    "P": "Positive",
    "POS": "Positive",
    "POSITIVE": "Positive",
    "+": "Positive",
    "N": "Negative",
    "NEG": "Negative",
    "NEGATIVE": "Negative",
    "-": "Negative",
}

_NUMBER = re.compile(r"\d+(?:\.\d+)?")
_RANGE = re.compile(r"\d+(?:\.\d+)?\s*(?:-|--|–|TO)\s*\d+(?:\.\d+)?")
# Unit names are matched as whole tokens. 'EV' as a bare substring also occurs inside 'LEVEL' and
# 'SEVEN', which would report an eV energy for a string that never stated one.
_UNIT_PATTERNS: tuple[tuple[str, str], ...] = (
    (r"(?<![A-Z])NCE(?![A-Z])", "NCE"),
    (r"%|(?<![A-Z])PERCENT(?![A-Z])", "percent"),
    (r"(?<![A-Z])EV(?![A-Z])", "eV"),
)

# A knife-edge fraction such as 1/3 must not fail the presence test on binary rounding alone.
_FRACTION_EPSILON = 1e-9


@dataclass(frozen=True)
class ConsensusKey:
    """The condition group one consensus spectrum is built over."""

    inchikey_skeleton: str | None
    ion_mode: str | None
    precursor_type: str | None
    instrument_class: str
    collision_energy_bin: str


@dataclass(frozen=True)
class ReferenceRecord:
    """One reference-library record with its acquisition conditions already normalized."""

    record_index: int
    peaks: list[tuple[float, float]]
    record_name: str | None = None
    inchikey: str | None = None
    inchikey_skeleton: str | None = None
    smiles: str | None = None
    formula: str | None = None
    ontology: str | None = None
    precursor_mz: float | None = None
    precursor_type: str | None = None
    ion_mode: str | None = None
    instrument_type: str | None = None
    instrument_class: str = UNKNOWN_INSTRUMENT_CLASS
    collision_energy_raw: str | None = None
    collision_energy_value: float | None = None
    collision_energy_unit: str | None = None
    rt_min: float | None = None
    ccs: float | None = None


@dataclass
class ConsensusSpectrum:
    """One merged spectrum plus the identity of the records it was merged from."""

    key: ConsensusKey
    peaks: list[tuple[float, float]]
    member_count: int
    record_names: list[str] = field(default_factory=list)
    inchikeys: list[str] = field(default_factory=list)
    smiles: str | None = None
    formula: str | None = None
    precursor_mz: float | None = None
    ontology: str | None = None


@dataclass
class ReferenceIngestReport:
    """Outcome of one reference-library ingest."""

    library_id: str
    records_read: int = 0
    records_skipped: int = 0
    consensus_spectra: int = 0
    blobs_written: int = 0
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def _is_missing(raw: str | None) -> bool:
    """Report whether a free-text field states the absence of a value."""
    return raw is None or raw.strip().upper() in _MISSING_TOKENS


def _text(raw: str | None) -> str | None:
    """Strip a free-text field, or return None when it carries no value."""
    return None if _is_missing(raw) else raw.strip()


def _stripped(raw: str | None) -> str | None:
    """Strip a free-text field, keeping a stated missing-value token such as 'nan'."""
    if raw is None or not raw.strip():
        return None
    return raw.strip()


def _squeeze(raw: str | None) -> str:
    """Uppercase a free-text field and drop every separator, so 'Q-Exactive' matches 'QEXACTIVE'."""
    if raw is None:
        return ""
    return "".join(character for character in raw.upper() if character.isascii() and character.isalnum())


def normalize_instrument_type(raw: str | None) -> str:
    """Map a free-text INSTRUMENTTYPE onto an analyzer class, or UNKNOWN when nothing matches."""
    squeezed = _squeeze(raw)
    if not squeezed:
        return UNKNOWN_INSTRUMENT_CLASS
    for token, instrument_class in INSTRUMENT_CLASSES:
        if token in squeezed:
            return instrument_class
    return UNKNOWN_INSTRUMENT_CLASS


def _collision_energy_unit(text: str) -> str | None:
    for pattern, unit in _UNIT_PATTERNS:
        if re.search(pattern, text):
            return unit
    return None


def normalize_collision_energy(raw: str | None) -> tuple[float | None, str | None]:
    """Split a free-text COLLISIONENERGY into a numeric value and, when stated, its unit.

    A missing or unparseable energy gives (None, None) and never 0.0. 'HCD' names a dissociation method
    rather than an energy unit, so '45HCD' yields a value with no unit instead of an assumed NCE percent.
    """
    if _is_missing(raw):
        return None, None
    text = raw.strip().upper()
    numbers = _NUMBER.findall(text)
    if not numbers:
        return None, None
    ramp = _RANGE.search(text)
    if ramp is not None:
        low, high = (float(number) for number in _NUMBER.findall(ramp.group(0))[:2])
        value = (low + high) / 2.0
    else:
        value = float(numbers[0])
    return value, _collision_energy_unit(text)


def collision_energy_bin(value: float | None, width: float = COLLISION_ENERGY_BIN_WIDTH) -> str:
    """Label the collision-energy bin a value falls in, or CE_UNKNOWN when there is no value."""
    if value is None:
        return UNKNOWN_COLLISION_ENERGY_BIN
    if width <= 0:
        raise ValueError(f"collision-energy bin width must be positive, got {width!r}")
    lower = math.floor(value / width) * width
    return f"CE{lower:g}-{lower + width:g}"


def inchikey_skeleton(inchikey: str | None) -> str | None:
    """Return the uppercased first block of an InChIKey, or None when there is no usable one."""
    if _is_missing(inchikey):
        return None
    block = inchikey.strip().upper().split("-")[0]
    # A first block that is not 14 ASCII letters is not a connectivity hash. Grouping on it would merge
    # records that share nothing but a malformed field.
    if len(block) != _SKELETON_LENGTH or not block.isascii() or not block.isalpha():
        return None
    return block


def normalize_ion_mode(raw: str | None) -> str | None:
    """Map a free-text IONMODE onto 'Positive' or 'Negative', or None when nothing matches."""
    if _is_missing(raw):
        return None
    return ION_MODES.get(raw.strip().upper())


def _collision_energy_from_comment(comment_tokens: dict[str, str]) -> str | None:
    """Recover a collision energy from the COMMENT tokens when the field itself is unusable."""
    folded = {key.upper(): value for key, value in comment_tokens.items()}
    for key in COLLISION_ENERGY_COMMENT_KEYS:
        candidate = folded.get(key)
        if not _is_missing(candidate) and _NUMBER.search(candidate.upper()):
            return candidate.strip()
    return None


def reference_record_from_msp(record: MspRecord) -> ReferenceRecord:
    """Map one MSP record onto a ReferenceRecord, normalizing its acquisition conditions."""
    fields = record.fields
    # collision_energy_raw keeps the string the value was read from, and the stated token itself when
    # nothing could be read from it, so 'the field said nan' stays distinguishable from 'no field'.
    raw_energy = _stripped(fields.get("COLLISIONENERGY"))
    if _is_missing(raw_energy):
        recovered = _collision_energy_from_comment(record.comment_tokens)
        if recovered is not None:
            raw_energy = recovered
    energy_value, energy_unit = normalize_collision_energy(raw_energy)
    key = _text(fields.get("INCHIKEY"))
    return ReferenceRecord(
        record_index=record.index,
        peaks=list(record.peaks),
        record_name=_text(fields.get("NAME")),
        inchikey=key,
        inchikey_skeleton=inchikey_skeleton(key),
        smiles=_text(fields.get("SMILES")),
        formula=_text(fields.get("FORMULA")),
        ontology=_text(fields.get("ONTOLOGY")),
        precursor_mz=parse_number(fields.get("PRECURSORMZ")),
        precursor_type=_text(fields.get("PRECURSORTYPE")),
        ion_mode=normalize_ion_mode(fields.get("IONMODE")),
        instrument_type=_text(fields.get("INSTRUMENTTYPE")),
        instrument_class=normalize_instrument_type(fields.get("INSTRUMENTTYPE")),
        collision_energy_raw=raw_energy,
        collision_energy_value=energy_value,
        collision_energy_unit=energy_unit,
        rt_min=parse_number(fields.get("RETENTIONTIME")),
        ccs=parse_number(fields.get("CCS")),
    )


def consensus_key_for(record: ReferenceRecord) -> ConsensusKey:
    """Derive the condition group of one record."""
    return ConsensusKey(
        inchikey_skeleton=record.inchikey_skeleton,
        ion_mode=record.ion_mode,
        precursor_type=record.precursor_type,
        instrument_class=record.instrument_class,
        collision_energy_bin=collision_energy_bin(record.collision_energy_value),
    )


def _frames(
    points: list[tuple[float, float, int]], bin_width: float
) -> Iterator[list[tuple[float, float, int]]]:
    """Group m/z-sorted (m/z, intensity, member) points into frames no wider than bin_width.

    Each frame is bounded by its own lowest m/z rather than by the previous point, so a dense ladder of
    fragments cannot chain into one arbitrarily wide frame.
    """
    frame: list[tuple[float, float, int]] = []
    start = 0.0
    for point in points:
        if frame and point[0] - start > bin_width:
            yield frame
            frame = []
        if not frame:
            start = point[0]
        frame.append(point)
    if frame:
        yield frame


def _unique(values: list[str | None]) -> list[str]:
    """Keep the distinct non-empty values in first-seen order."""
    seen: dict[str, None] = {}
    for value in values:
        if value:
            seen.setdefault(value, None)
    return list(seen)


def _dominant(values: list[str | None]) -> str | None:
    """Return the most frequent non-empty value, breaking ties by first appearance.

    The mode rather than the first value, so one mis-annotated replicate cannot rename a whole group.
    """
    counts: dict[str, int] = {}
    first_seen: dict[str, int] = {}
    for index, value in enumerate(values):
        if not value:
            continue
        counts[value] = counts.get(value, 0) + 1
        first_seen.setdefault(value, index)
    if not counts:
        return None
    return max(counts, key=lambda value: (counts[value], -first_seen[value]))


def _consensus_of(
    key: ConsensusKey,
    members: list[ReferenceRecord],
    bin_width: float,
    minimum_member_fraction: float,
) -> ConsensusSpectrum:
    """Merge the member spectra of one condition group onto a shared m/z grid."""
    points: list[tuple[float, float, int]] = []
    for index, member in enumerate(members):
        # Replicate absolute intensities are not comparable, so every member is put on the same scale
        # before anything is averaged.
        for mz, intensity in normalize_to_base_peak(member.peaks):
            points.append((mz, intensity, index))
    points.sort()
    required = minimum_member_fraction * len(members) - _FRACTION_EPSILON
    merged: list[tuple[float, float]] = []
    for frame in _frames(points, bin_width):
        contributors: dict[int, float] = {}
        weight = 0.0
        weighted_mz = 0.0
        for mz, intensity, index in frame:
            contributors[index] = contributors.get(index, 0.0) + intensity
            weighted_mz += mz * intensity
            weight += intensity
        if len(contributors) < required:
            continue
        # Mean over the members that HAVE the peak, not over every member: dividing by the member count
        # would push a peak seen in a minority of replicates toward the noise floor, and would make the
        # consensus intensity depend on how many replicates happened to be deposited.
        merged.append(
            (
                weighted_mz / weight if weight > 0 else frame[0][0],
                sum(contributors.values()) / len(contributors),
            )
        )
    merged.sort()
    precursors = [member.precursor_mz for member in members if member.precursor_mz is not None]
    return ConsensusSpectrum(
        key=key,
        peaks=normalize_to_base_peak(merged),
        member_count=len(members),
        record_names=_unique([member.record_name for member in members]),
        inchikeys=_unique([member.inchikey for member in members]),
        smiles=_dominant([member.smiles for member in members]),
        formula=_dominant([member.formula for member in members]),
        precursor_mz=sum(precursors) / len(precursors) if precursors else None,
        ontology=_dominant([member.ontology for member in members]),
    )


def _group_identity(key: ConsensusKey, record: ReferenceRecord) -> tuple[ConsensusKey, int | None]:
    """Identify the group a record belongs to, keeping skeleton-less records apart.

    Two records without an InChIKey are not known to be the same compound, so merging them would
    fabricate a consensus spectrum across unrelated compounds.
    """
    return (key, record.record_index if key.inchikey_skeleton is None else None)


def build_consensus(
    records: Iterator[ReferenceRecord],
    *,
    bin_width: float = CONSENSUS_MZ_BIN_WIDTH,
    minimum_member_fraction: float = CONSENSUS_MINIMUM_MEMBER_FRACTION,
) -> Iterator[ConsensusSpectrum]:
    """Merge reference records into one consensus spectrum per condition group.

    A consensus peak is kept only when it appears in at least minimum_member_fraction of the members, so
    a single noisy replicate cannot invent a peak; a one-member group satisfies that test trivially and
    passes through unchanged apart from base-peak normalization.

    Members are buffered by group, which suits a bounded caller-supplied record set. A whole production
    library is grouped by iter_library_consensus instead, which streams the groups out of the database.
    """
    groups: dict[tuple[ConsensusKey, int | None], tuple[ConsensusKey, list[ReferenceRecord]]] = {}
    for record in records:
        key = consensus_key_for(record)
        groups.setdefault(_group_identity(key, record), (key, []))[1].append(record)
    for key, members in groups.values():
        yield _consensus_of(key, members, bin_width, minimum_member_fraction)


def _put_blob(connection, payload: dict) -> tuple[str, bool]:
    """Store one content-addressed zlib-JSON spectrum payload, reporting whether it was new.

    The serialization stays byte-identical to ingest._put_blob: both writers share spectrum_blob, and a
    differing separator or key order would store the same spectrum twice under two digests.
    """
    raw = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode("utf-8")
    payload_sha = hashlib.sha256(raw).hexdigest()
    cursor = connection.execute(
        "INSERT OR IGNORE INTO spectrum_blob(payload_sha256, compression, uncompressed_bytes, payload) "
        "VALUES (?, 'zlib-json', ?, ?)",
        (payload_sha, len(raw), zlib.compress(raw, level=9)),
    )
    return payload_sha, cursor.rowcount == 1


def _peaks_from_blob(connection, payload_sha256: str) -> list[tuple[float, float]]:
    row = connection.execute(
        "SELECT compression, payload FROM spectrum_blob WHERE payload_sha256 = ?", (payload_sha256,)
    ).fetchone()
    if row is None or row["compression"] != "zlib-json":
        return []
    payload = json.loads(zlib.decompress(row["payload"]).decode("utf-8"))
    return [(float(mz), float(intensity)) for mz, intensity in payload.get("peaks", [])]


# How much of a file to look at before deciding it is not MSP text. Large enough to cross the header
# comments some vendors prepend, small enough that a 700 MB library costs nothing to reject.
_FORMAT_SNIFF_BYTES = 64 * 1024

# Field names that identify a record as MSP. NAME and Num Peaks are the two an MSP record cannot omit
# and still be one; finding any of these settles the question.
_MSP_MARKER_FIELDS = ("name:", "num peaks:", "numpeaks:", "precursormz:")


def msp_format_problem(path: Path) -> str | None:
    """Report why a file is not an MSP text library, or None when it looks like one.

    This ingest is the only door onto the library half of the catalog, and it used to open for
    anything. Handed a binary MS-DIAL .lbm2 it hashed the file, streamed it through the MSP reader,
    obtained pseudo-records from the bytes, skipped every one of them, wrote a library row with no
    records and returned valid with no errors and exit status 0. Verified with 102,400 bytes of random
    binary: records_read 400, records_skipped 400, errors []. That silence is why the library half of a
    real catalog can stay empty without anyone noticing.
    """
    try:
        head = path.read_bytes()[:_FORMAT_SNIFF_BYTES]
    except OSError as error:
        return f"could not be read: {error}"
    if not head:
        return "is empty"
    if b"\x00" in head:
        return (
            "contains NUL bytes in its first 64 KB, so it is a binary file rather than MSP text. "
            "An MS-DIAL .lbm2 or .msp2 library is binary and cannot be ingested here; export it to "
            "text MSP first"
        )
    try:
        text = head.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = head.decode("latin-1", errors="replace")
    lowered = text.casefold()
    if not any(marker in lowered for marker in _MSP_MARKER_FIELDS):
        return (
            "carries none of the MSP record fields (NAME, PRECURSORMZ, Num Peaks) in its first 64 KB, "
            "so it does not look like an MSP library"
        )
    return None


def _skip_reason(record: MspRecord) -> str | None:
    """Report why a record cannot be ingested, or None when it can."""
    if parse_number(record.fields.get("NUM PEAKS"), int) is None:
        return "unparseable Num Peaks"
    if not record.peaks:
        return "no peaks"
    if parse_number(record.fields.get("PRECURSORMZ")) is None:
        return "no precursor m/z"
    return None


def _record_from_row(connection, row) -> ReferenceRecord:
    return ReferenceRecord(
        record_index=row["library_record_index"],
        peaks=_peaks_from_blob(connection, row["payload_sha256"]),
        record_name=row["record_name"],
        inchikey=row["inchikey"],
        inchikey_skeleton=row["inchikey_skeleton"],
        smiles=row["smiles"],
        formula=row["formula"],
        ontology=row["ontology"],
        precursor_mz=row["precursor_mz"],
        precursor_type=row["precursor_type"],
        ion_mode=row["ion_mode"],
        instrument_type=row["instrument_type"],
        instrument_class=row["instrument_class"] or UNKNOWN_INSTRUMENT_CLASS,
        collision_energy_raw=row["collision_energy_raw"],
        collision_energy_value=row["collision_energy_value"],
        collision_energy_unit=row["collision_energy_unit"],
        rt_min=row["rt_min"],
        ccs=row["ccs"],
    )


def _iter_consensus(
    connection,
    library_id: str,
    *,
    bin_width: float,
    minimum_member_fraction: float,
) -> Iterator[ConsensusSpectrum]:
    """Stream one library's consensus spectra, holding only the members of the current group.

    The rows are ordered by the group key so groups arrive contiguously. Ordering by
    collision_energy_value rather than by its bin label is enough: the bin is monotone in the value, so
    every row of one bin is still contiguous, and a NULL value sorts as one run of CE_UNKNOWN.
    """
    rows = connection.execute(
        """SELECT library_record_index, record_name, inchikey, inchikey_skeleton, smiles, formula,
                  ontology, precursor_mz, precursor_type, ion_mode, instrument_type, instrument_class,
                  collision_energy_raw, collision_energy_value, collision_energy_unit, rt_min, ccs,
                  payload_sha256
           FROM reference_spectrum
           WHERE library_id = ?
           ORDER BY inchikey_skeleton, ion_mode, precursor_type, instrument_class,
                    collision_energy_value, library_record_index""",
        (library_id,),
    )
    current: ConsensusKey | None = None
    members: list[ReferenceRecord] = []
    for row in rows:
        record = _record_from_row(connection, row)
        key = consensus_key_for(record)
        if members and (key != current or key.inchikey_skeleton is None):
            yield _consensus_of(current, members, bin_width, minimum_member_fraction)
            members = []
        current = key
        members.append(record)
    if members:
        yield _consensus_of(current, members, bin_width, minimum_member_fraction)


def iter_library_consensus(
    database: str | Path,
    library_id: str,
    *,
    bin_width: float = CONSENSUS_MZ_BIN_WIDTH,
    minimum_member_fraction: float = CONSENSUS_MINIMUM_MEMBER_FRACTION,
) -> Iterator[ConsensusSpectrum]:
    """Stream the skeleton consensus spectra of one ingested library."""
    connection = connect(database)
    try:
        yield from _iter_consensus(
            connection,
            library_id,
            bin_width=bin_width,
            minimum_member_fraction=minimum_member_fraction,
        )
    finally:
        connection.close()


def _consensus_payload(
    spectrum: ConsensusSpectrum, library_id: str, library_kind: str, bin_width: float, fraction: float
) -> dict:
    """Shape one consensus spectrum for the content-addressed blob store.

    Every parameter travels with the spectrum. The threshold and bin widths are conventions rather than
    measurements, so a later recalibration has to be able to tell which convention produced this payload.
    """
    return {
        "kind": "skeleton_consensus",
        "library_id": library_id,
        "library_kind": library_kind,
        "consensus_key": {
            "inchikey_skeleton": spectrum.key.inchikey_skeleton,
            "ion_mode": spectrum.key.ion_mode,
            "precursor_type": spectrum.key.precursor_type,
            "instrument_class": spectrum.key.instrument_class,
            "collision_energy_bin": spectrum.key.collision_energy_bin,
        },
        "member_count": spectrum.member_count,
        "record_names": spectrum.record_names,
        "inchikeys": spectrum.inchikeys,
        "smiles": spectrum.smiles,
        "formula": spectrum.formula,
        "precursor_mz": spectrum.precursor_mz,
        "ontology": spectrum.ontology,
        "peaks": [[mz, intensity] for mz, intensity in spectrum.peaks],
        "parameters": {
            "mz_bin_width": bin_width,
            "minimum_member_fraction": fraction,
            "collision_energy_bin_width": COLLISION_ENERGY_BIN_WIDTH,
        },
    }


def ingest_reference_library(
    database: str | Path,
    msp_path: str | Path,
    *,
    library_name: str,
    library_version: str | None = None,
    library_kind: str = EXPERIMENTAL_LIBRARY_KIND,
    source_uri: str | None = None,
    license: str | None = None,
    tool_run_id: str | None = None,
    consensus: bool = True,
    limit: int | None = None,
    precursor_mz_range: tuple[float, float] | None = None,
) -> ReferenceIngestReport:
    """Stream one MSP reference library into the catalog and build its skeleton consensus spectra.

    The file is streamed record by record and the consensus pass reads the ingested rows back in group
    order, so neither pass holds the library in memory.

    'limit' caps the number of in-scope records examined and 'precursor_mz_range' (inclusive) selects a
    mass window, both applied before any record is compressed or written, so a bounded pre-test on a
    large library stays bounded. A record outside the requested window is out of scope rather than
    rejected, so it is counted neither as read nor as skipped; a record with no peaks, no precursor m/z
    or an unparseable Num Peaks is counted as skipped. record_count on the library row is therefore what
    this call ingested, which is the whole file only when neither bound was given.
    """
    if library_kind not in LIBRARY_KINDS:
        raise ValueError(f"library_kind must be one of {LIBRARY_KINDS}, got {library_kind!r}")
    path = Path(msp_path)
    initialize(database)
    library_id = make_id("reference-library", library_name, library_version or "unversioned")
    report = ReferenceIngestReport(library_id=library_id)
    # Refused before anything is hashed, written or registered. A library row for a file that could
    # never yield a record is worse than no row: it makes an empty library look like an ingested one.
    problem = msp_format_problem(path)
    if problem is not None:
        report.errors.append(f"{path.name} {problem}.")
        return report
    file_sha256 = sha256_file(path)
    byte_size = path.stat().st_size
    written = 0
    with transaction(database) as connection:
        existing = connection.execute(
            "SELECT sha256 FROM reference_library WHERE library_id = ?", (library_id,)
        ).fetchone()
        if existing is not None and existing["sha256"] not in (None, file_sha256):
            report.warnings.append(
                f"library {library_id} was ingested from a file with digest {existing['sha256']}; "
                "records of the previous file are updated in place but not removed"
            )
        connection.execute(
            """INSERT INTO reference_library(
                   library_id, library_name, library_version, library_kind, source_uri, license,
                   sha256, byte_size, record_count, tool_run_id)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, 0, ?)
               ON CONFLICT(library_id) DO UPDATE SET
                   library_name = excluded.library_name,
                   library_version = excluded.library_version,
                   library_kind = excluded.library_kind,
                   source_uri = excluded.source_uri,
                   license = excluded.license,
                   sha256 = excluded.sha256,
                   byte_size = excluded.byte_size,
                   tool_run_id = excluded.tool_run_id""",
            (
                library_id,
                library_name,
                library_version,
                library_kind,
                source_uri,
                license,
                file_sha256,
                byte_size,
                tool_run_id,
            ),
        )
        for record in read_msp(path):
            if limit is not None and report.records_read >= limit:
                break
            precursor_mz = parse_number(record.fields.get("PRECURSORMZ"))
            if (
                precursor_mz_range is not None
                and precursor_mz is not None
                and not precursor_mz_range[0] <= precursor_mz <= precursor_mz_range[1]
            ):
                continue
            report.records_read += 1
            reason = _skip_reason(record)
            if reason is not None:
                report.records_skipped += 1
                continue
            reference = reference_record_from_msp(record)
            payload_sha256, is_new_blob = _put_blob(
                connection, {"fields": record.fields, "peaks": record.peaks}
            )
            if is_new_blob:
                report.blobs_written += 1
            reference_spectrum_id = make_id("reference-spectrum", library_id, record.index)
            connection.execute(
                """INSERT INTO reference_spectrum(
                       reference_spectrum_id, library_id, library_record_index, record_name, inchikey,
                       inchikey_skeleton, smiles, formula, ontology, precursor_mz, precursor_type,
                       ion_mode, instrument_type, instrument_class, collision_energy_raw,
                       collision_energy_value, collision_energy_unit, rt_min, ccs, peak_count,
                       payload_sha256)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                   ON CONFLICT(reference_spectrum_id) DO UPDATE SET
                       record_name = excluded.record_name,
                       inchikey = excluded.inchikey,
                       inchikey_skeleton = excluded.inchikey_skeleton,
                       smiles = excluded.smiles,
                       formula = excluded.formula,
                       ontology = excluded.ontology,
                       precursor_mz = excluded.precursor_mz,
                       precursor_type = excluded.precursor_type,
                       ion_mode = excluded.ion_mode,
                       instrument_type = excluded.instrument_type,
                       instrument_class = excluded.instrument_class,
                       collision_energy_raw = excluded.collision_energy_raw,
                       collision_energy_value = excluded.collision_energy_value,
                       collision_energy_unit = excluded.collision_energy_unit,
                       rt_min = excluded.rt_min,
                       ccs = excluded.ccs,
                       peak_count = excluded.peak_count,
                       payload_sha256 = excluded.payload_sha256""",
                (
                    reference_spectrum_id,
                    library_id,
                    record.index,
                    reference.record_name,
                    reference.inchikey,
                    reference.inchikey_skeleton,
                    reference.smiles,
                    reference.formula,
                    reference.ontology,
                    reference.precursor_mz,
                    reference.precursor_type,
                    reference.ion_mode,
                    reference.instrument_type,
                    reference.instrument_class,
                    reference.collision_energy_raw,
                    reference.collision_energy_value,
                    reference.collision_energy_unit,
                    reference.rt_min,
                    reference.ccs,
                    len(record.peaks),
                    payload_sha256,
                ),
            )
            written += 1
        connection.execute(
            "UPDATE reference_library SET record_count = ? WHERE library_id = ?", (written, library_id)
        )
        if consensus:
            for spectrum in _iter_consensus(
                connection,
                library_id,
                bin_width=CONSENSUS_MZ_BIN_WIDTH,
                minimum_member_fraction=CONSENSUS_MINIMUM_MEMBER_FRACTION,
            ):
                payload_sha256, is_new_blob = _put_blob(
                    connection,
                    _consensus_payload(
                        spectrum,
                        library_id,
                        library_kind,
                        CONSENSUS_MZ_BIN_WIDTH,
                        CONSENSUS_MINIMUM_MEMBER_FRACTION,
                    ),
                )
                if is_new_blob:
                    report.blobs_written += 1
                # The payload used to go into the content-addressed store with nothing pointing at it.
                # It was computed, compressed, counted in the report, and then unreachable: no table
                # named it and no query could find it. A row makes it retrievable and makes the
                # reported count mean something.
                key = spectrum.key
                connection.execute(
                    """INSERT INTO reference_consensus_spectrum(
                           reference_consensus_spectrum_id, library_id, inchikey_skeleton, ion_mode,
                           precursor_type, instrument_class, collision_energy_bin, member_count,
                           precursor_mz, formula, ontology, smiles, peak_count, mz_bin_width,
                           minimum_member_fraction, payload_sha256)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                       ON CONFLICT(library_id, inchikey_skeleton, ion_mode, precursor_type,
                                   instrument_class, collision_energy_bin) DO UPDATE SET
                           member_count = excluded.member_count,
                           precursor_mz = excluded.precursor_mz,
                           formula = excluded.formula,
                           ontology = excluded.ontology,
                           smiles = excluded.smiles,
                           peak_count = excluded.peak_count,
                           mz_bin_width = excluded.mz_bin_width,
                           minimum_member_fraction = excluded.minimum_member_fraction,
                           payload_sha256 = excluded.payload_sha256""",
                    (
                        make_id(
                            "reference-consensus",
                            library_id,
                            key.inchikey_skeleton or "",
                            key.ion_mode or "",
                            key.precursor_type or "",
                            key.instrument_class,
                            key.collision_energy_bin,
                        ),
                        library_id,
                        key.inchikey_skeleton,
                        key.ion_mode,
                        key.precursor_type,
                        key.instrument_class,
                        key.collision_energy_bin,
                        spectrum.member_count,
                        spectrum.precursor_mz,
                        spectrum.formula,
                        spectrum.ontology,
                        spectrum.smiles,
                        len(spectrum.peaks),
                        CONSENSUS_MZ_BIN_WIDTH,
                        CONSENSUS_MINIMUM_MEMBER_FRACTION,
                        payload_sha256,
                    ),
                )
                report.consensus_spectra += 1
    # A file that parsed but yielded nothing usable is a failed ingest, not a successful empty one.
    # The format sniff catches a wholly wrong file; this catches one that is MSP-shaped but whose
    # records all lack peaks or a precursor m/z, which leaves the same empty library behind.
    if report.records_read and report.records_read == report.records_skipped:
        report.errors.append(
            f"All {report.records_read} records of {path.name} were skipped, so the library is empty. "
            "Every record lacked peaks, a precursor m/z, or a parseable Num Peaks."
        )
    elif report.records_read == 0:
        report.errors.append(
            f"No records were read from {path.name}"
            + (" within the requested precursor m/z window." if precursor_mz_range else ".")
        )
    return report
