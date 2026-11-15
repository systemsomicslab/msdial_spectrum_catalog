"""Select representative MS/MS from one reanalysis unit and write them as MS-Emblator input.

MS-Emblator (systemsomicslab/msemblator) runs outside MS-DIAL, and its results rejoin the unit's
alignment by Peak ID. Three properties of its input decide whether that join is sound:

- MS-Emblator keys each spectrum by its NAME and a repeated NAME silently replaces the spectrum before
  it. MS-DIAL's NAME is the annotation ("Unknown", the same lipid twice), so every exported NAME is
  rebuilt from the unit and the alignment ID and is unique by construction.
- The representative record carries the identity of the library match that won it (FORMULA, INCHIKEY,
  SMILES, ONTOLOGY). An in-silico tool that is handed the answer is not evaluated, and a mode-3 run would
  take the library's formula as given, so none of those fields is written; the manifest keeps them, for
  the benchmark only.
- The source files are recorded by SHA-256, so a later ingest can refuse results that no longer match the
  outputs they were exported from.

An optional peak cleanup (max_entropy, in `spectrum_cleanup`) drops each spectrum's weakest ions until its
spectral entropy is at most a ceiling. The manifest keeps the raw and the written peak count and entropy, and
the selection record states the rule, so a run says which peaks the tools were given.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from .parsers import read_msp, read_tsv, sha256_file
from .annotator_settings import read_annotator_settings
from .spectrum_cleanup import clean_by_entropy, spectral_entropy

TIERS = ("LBM2", "VS21 strict", "VS21 relaxed", "GLACIER in-silico", "unassigned")
SELECTIONS = {
    # A spectrum a rule-based lipid annotation or a strict reference match already explains is not sent.
    "unresolved": frozenset({"VS21 relaxed", "GLACIER in-silico", "unassigned"}),
    # Strict reference matches, whose library identity is the truth a benchmark compares against.
    "benchmark": frozenset({"VS21 strict"}),
}
MINIMUM_PEAKS = 3
_UNSAFE_NAME = re.compile(r"[^A-Za-z0-9_.-]")


@dataclass(frozen=True)
class ExportedSpectrum:
    name: str
    unit_id: str
    alignment_id: int
    ion_mode: str
    precursor_mz: float
    precursor_type: str
    retention_time_min: float | None
    peak_count: int
    peaks_sha256: str
    tier: str
    raw_peak_count: int
    spectral_entropy_raw: float
    spectral_entropy: float
    reference_name: str
    reference_formula: str
    reference_inchikey: str
    reference_ontology: str


def _sampling_rank(unit_id: str, alignment_id: int) -> int:
    return int.from_bytes(hashlib.sha256(f"{unit_id}:{alignment_id}".encode()).digest()[:8], "big")


def annotator_tiers(output: Path) -> dict[str, str]:
    """Map each MSP annotator of the run to its evidence tier, from msp_annotator_settings.tsv.

    Only annotator_id, evidence_tier and library_kind are read: the file also holds each library's
    absolute path, which must not travel with the export.
    """
    mapping: dict[str, str] = {}
    path = output / "msp_annotator_settings.tsv"
    if not path.exists():
        return mapping
    for row in read_annotator_settings(path):
        evidence = row.get("evidence_tier", "").lower()
        if row.get("library_kind", "").lower() == "predicted":
            mapping[row["annotator_id"]] = "GLACIER in-silico"
        elif "high" in evidence:
            mapping[row["annotator_id"]] = "VS21 strict"
        elif "low" in evidence:
            mapping[row["annotator_id"]] = "VS21 relaxed"
    return mapping


def candidate_tier(row: dict[str, str], mapping: dict[str, str]) -> str:
    """The tier one representative candidate earns: its annotator alone is not enough."""
    tag = row.get("annotation_tag_vs1", "")
    evidence = row.get("evidence_source", "")
    if row.get("annotator_id", "").lower().startswith("lbmdb:") and tag in {"410", "420", "430"} \
            and evidence == "RuleBased":
        return "LBM2"
    tier = mapping.get(row.get("annotator_id", ""), "unassigned")
    if tier == "GLACIER in-silico" and tag == "450" and evidence == "InSilico":
        return tier
    if tier in {"VS21 strict", "VS21 relaxed"} and tag == "430" and evidence == "ReferenceSpectrum":
        return tier
    return "unassigned"


def latest_alignment_stem(output: Path) -> str:
    matches = sorted(output.glob("AlignResult-*.mdmsp"), key=lambda path: path.stat().st_mtime)
    if not matches:
        raise FileNotFoundError(f"No AlignResult-*.mdmsp in {output}")
    return matches[-1].name.removesuffix(".mdmsp")


def representative_tiers(output: Path, stem: str) -> dict[int, tuple[str, dict[str, str]]]:
    path = output / f"{stem}.mdcandidate.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}; the run must export annotation candidates")
    mapping = annotator_tiers(output)
    tiers: dict[int, tuple[str, dict[str, str]]] = {}
    for row_number, row in read_tsv(path):
        if row.get("is_representative", "").lower() != "true":
            continue
        spot = int(row["alignment_master_id"])
        if spot in tiers:
            raise ValueError(f"{path}:{row_number}: alignment {spot} has two representative candidates")
        tiers[spot] = (candidate_tier(row, mapping), row)
    return tiers


def cleanup_rule(max_entropy: float | None) -> dict:
    """The peak cleanup as the selection record and the plan state it."""
    if max_entropy is None:
        return {"rule": None}
    return {"rule": "spectral_entropy", "max_entropy": max_entropy, "logarithm": "natural",
            "drops": "the weakest ions, equal intensities together, until the entropy is at most max_entropy",
            "reference": "Li et al., Nat. Methods 18, 1524-1531 (2021)"}


def select_spectra(output: Path, unit_id: str, selection: str = "unresolved",
                   max_spectra: int | None = None,
                   max_entropy: float | None = None) -> tuple[list[tuple[ExportedSpectrum, list]], dict]:
    """Return the selected spectra with their peaks, and the counts that explain the selection.

    With max_entropy, each selected spectrum's weakest ions are dropped until its spectral entropy is at
    most max_entropy; the minimum-peak rule is then checked again."""
    if selection not in SELECTIONS:
        raise ValueError(f"selection is one of {sorted(SELECTIONS)}")
    if max_entropy is not None and max_entropy <= 0:
        raise ValueError("max_entropy must be positive")
    stem = latest_alignment_stem(output)
    tiers = representative_tiers(output, stem)
    wanted = SELECTIONS[selection]
    counts = {"records": 0, "no_peak_id": 0, "isotope": 0, "too_few_peaks": 0, "not_selected_tier": 0,
              "selected_before_cap": 0}
    cleanup = {**cleanup_rule(max_entropy), "spectra_cleaned": 0, "peaks_removed": 0}
    by_tier = {tier: 0 for tier in TIERS}
    chosen = []
    for record in read_msp(output / f"{stem}.mdmsp"):
        counts["records"] += 1
        peak_id = record.comment_tokens.get("PEAKID")
        if peak_id is None:
            counts["no_peak_id"] += 1
            continue
        if record.comment_tokens.get("ISOTOPE", "M+0") != "M+0":
            counts["isotope"] += 1
            continue
        if len(record.peaks) < MINIMUM_PEAKS:
            counts["too_few_peaks"] += 1
            continue
        spot = int(peak_id)
        tier, candidate = tiers.get(spot, ("unassigned", {}))
        by_tier[tier] += 1
        if tier not in wanted:
            counts["not_selected_tier"] += 1
            continue
        peaks = record.peaks
        entropy_raw = entropy = spectral_entropy(peaks)
        if max_entropy is not None:
            cleaned = clean_by_entropy(peaks, max_entropy)
            if cleaned.removed:
                cleanup["spectra_cleaned"] += 1
                cleanup["peaks_removed"] += cleaned.removed
            peaks, entropy = cleaned.peaks, cleaned.entropy_after
            if len(peaks) < MINIMUM_PEAKS:  # unreachable for a ceiling above ln 3; kept as a guard
                counts["too_few_peaks"] += 1
                continue
        peaks_text = "\n".join(f"{mz:.5f}\t{intensity:.0f}" for mz, intensity in peaks)
        fields = record.fields
        reference = candidate if selection == "benchmark" else {}
        chosen.append((ExportedSpectrum(
            name=_UNSAFE_NAME.sub("_", f"{unit_id}_{spot}"),
            unit_id=unit_id,
            alignment_id=spot,
            ion_mode=fields.get("IONMODE", ""),
            precursor_mz=float(fields["PRECURSORMZ"]),
            precursor_type=fields.get("PRECURSORTYPE", ""),
            retention_time_min=float(fields["RETENTIONTIME"]) if fields.get("RETENTIONTIME") else None,
            peak_count=len(peaks),
            peaks_sha256=hashlib.sha256(peaks_text.encode()).hexdigest(),
            tier=tier,
            raw_peak_count=len(record.peaks),
            spectral_entropy_raw=round(entropy_raw, 4),
            spectral_entropy=round(entropy, 4),
            reference_name=reference.get("name", ""),
            reference_formula=reference.get("formula", ""),
            reference_inchikey=reference.get("inchikey", ""),
            reference_ontology=reference.get("ontology", ""),
        ), peaks))
    counts["selected_before_cap"] = len(chosen)
    if max_spectra is not None and len(chosen) > max_spectra:
        chosen = sorted(chosen, key=lambda item: _sampling_rank(unit_id, item[0].alignment_id))[:max_spectra]
        chosen.sort(key=lambda item: item[0].alignment_id)
    counts["selected"] = len(chosen)
    summary = {
        "unit_id": unit_id,
        "selection": selection,
        "selected_tiers": sorted(wanted),
        "minimum_peaks": MINIMUM_PEAKS,
        "max_spectra": max_spectra,
        "peak_cleanup": cleanup,
        "alignment_stem": stem,
        "counts": counts,
        "representative_spectra_by_tier": by_tier,
        "sources": {
            name: sha256_file(output / name)
            for name in (f"{stem}.mdmsp", f"{stem}.mdcandidate.tsv")
            if (output / name).exists()
        },
    }
    return chosen, summary


def write_msp(path: Path, spectra: list[tuple[ExportedSpectrum, list]]) -> None:
    """Write MS-Emblator input: identity-free records named by unit and alignment ID."""
    names = [spectrum.name for spectrum, _ in spectra]
    if len(set(names)) != len(names):
        raise ValueError("Exported NAMEs must be unique; MS-Emblator would silently drop the repeats")
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for spectrum, peaks in spectra:
            handle.write(f"NAME: {spectrum.name}\n")
            handle.write(f"PRECURSORMZ: {spectrum.precursor_mz}\n")
            handle.write(f"PRECURSORTYPE: {spectrum.precursor_type}\n")
            handle.write(f"IONMODE: {spectrum.ion_mode}\n")
            if spectrum.retention_time_min is not None:
                handle.write(f"RETENTIONTIME: {spectrum.retention_time_min}\n")
            handle.write(f"COMMENT: |UNIT={spectrum.unit_id}|ALIGNMENTID={spectrum.alignment_id}|\n")
            handle.write(f"Num Peaks: {len(peaks)}\n")
            for mz, intensity in peaks:
                handle.write(f"{mz:.5f}\t{intensity:.0f}\n")
            handle.write("\n")


def manifest_name(unit_id: str, selection: str) -> str:
    """The manifest's file name. The benchmark manifest holds the private library's identities (the truth
    the benchmark compares against), so its name says it is local-only."""
    suffix = "_manifest.local-vs21.tsv" if selection == "benchmark" else "_manifest.tsv"
    return f"{unit_id}_{selection}{suffix}"


def export_unit(output: Path, unit_id: str, destination: Path, selection: str = "unresolved",
                max_spectra: int | None = None, max_entropy: float | None = None) -> dict:
    """Write <unit>_<selection>.msp, its manifest TSV and a selection JSON into destination."""
    spectra, summary = select_spectra(output, unit_id, selection, max_spectra, max_entropy)
    destination.mkdir(parents=True, exist_ok=True)
    base = destination / f"{unit_id}_{selection}"
    write_msp(base.with_suffix(".msp"), spectra)
    fields = list(ExportedSpectrum.__dataclass_fields__)
    summary["manifest"] = manifest_name(unit_id, selection)
    with (destination / summary["manifest"]).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for spectrum, _ in spectra:
            writer.writerow(asdict(spectrum))
    summary["msp"] = base.with_suffix(".msp").name
    summary["msp_sha256"] = sha256_file(base.with_suffix(".msp"))
    base.with_name(base.name + "_selection.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary
