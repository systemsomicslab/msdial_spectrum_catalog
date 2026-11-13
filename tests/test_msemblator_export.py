"""Tests for the MS-Emblator export: unique names, no leaked identity, and the selection rules."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.msemblator_export import export_unit
from msdial_spectrum_catalog.parsers import read_msp

STEM = "AlignResult-2026"
SETTINGS = (
    "annotator_id\tmsp_file_path\tevidence_tier\tlibrary_kind\n"
    "msp_high_quality\tX:\\private\\lib.msp\tMSP high quality\tacquired\n"
    "msp_low_quality\tX:\\private\\lib.msp\tMSP low quality candidate\tacquired\n"
    "msp_in_silico\tX:\\private\\glacier.msp\tIn-silico MSP (GLACIER)\tpredicted\n"
)
CANDIDATE_HEADER = ["alignment_master_id", "is_representative", "annotator_id", "annotation_tag_vs1",
                    "evidence_source", "name", "formula", "inchikey", "ontology"]
# Spot 1: strict reference match; 2: LBM2; 3: relaxed; 4: precursor-only (unassigned); 5: too few peaks;
# 6: an M+1 isotope; 7: no candidate at all.
CANDIDATES = [
    ["1", "true", "msp_high_quality", "430", "ReferenceSpectrum", "Oxalic acid", "C2H2O4",
     "MUBZPKHOEPUJKR-UHFFFAOYSA-N", "Dicarboxylic acids"],
    ["1", "false", "msp_low_quality", "430", "ReferenceSpectrum", "Other", "", "", ""],
    ["2", "true", "LbmDB: test", "410", "RuleBased", "PC 34:1", "", "", "PC"],
    ["3", "true", "msp_low_quality", "430", "ReferenceSpectrum", "Unknown", "", "", ""],
    ["4", "true", "msp_high_quality", "530", "PrecursorOnly", "Unknown", "", "", ""],
    ["5", "true", "msp_low_quality", "430", "ReferenceSpectrum", "Unknown", "", "", ""],
    ["6", "true", "msp_low_quality", "430", "ReferenceSpectrum", "Unknown", "", "", ""],
]


def _record(spot: int, peaks: int, isotope: str = "M+0") -> str:
    lines = [
        "NAME: Unknown",
        "PRECURSORMZ: 300.1",
        "PRECURSORTYPE: [M+H]+",
        "IONMODE: Positive",
        "RETENTIONTIME: 5.5",
        "FORMULA: C2H2O4",
        "INCHIKEY: MUBZPKHOEPUJKR-UHFFFAOYSA-N",
        "SMILES: OC(=O)C(O)=O",
        "ONTOLOGY: Dicarboxylic acids",
        f"COMMENT: |PEAKID={spot}|ISOTOPE={isotope}|EVIDENCE=Test",
        f"Num Peaks: {peaks}",
    ]
    lines += [f"{100 + i}.0\t{10 * (i + 1)}" for i in range(peaks)]
    return "\n".join(lines) + "\n\n"


class MsemblatorExportTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.output = self.root / "output"
        self.output.mkdir()
        (self.output / "msp_annotator_settings.tsv").write_text(SETTINGS, encoding="utf-8")
        with (self.output / f"{STEM}.mdcandidate.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(CANDIDATE_HEADER)
            writer.writerows(CANDIDATES)
        records = [_record(1, 5), _record(2, 5), _record(3, 5), _record(4, 5), _record(5, 2),
                   _record(6, 5, "M+1"), _record(7, 5)]
        (self.output / f"{STEM}.mdmsp").write_text("".join(records), encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def _export(self, selection, max_spectra=None):
        destination = self.root / "out"
        summary = export_unit(self.output, "unit1", destination, selection, max_spectra)
        records = list(read_msp(destination / f"unit1_{selection}.msp"))
        with (destination / f"unit1_{selection}_manifest.tsv").open(encoding="utf-8") as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        return summary, records, manifest

    def test_unresolved_takes_relaxed_unassigned_and_uncandidated_spectra(self):
        summary, records, manifest = self._export("unresolved")
        self.assertEqual([r.fields["NAME"] for r in records], ["unit1_3", "unit1_4", "unit1_7"])
        self.assertEqual([m["tier"] for m in manifest], ["VS21 relaxed", "unassigned", "unassigned"])
        self.assertEqual(summary["counts"]["too_few_peaks"], 1)
        self.assertEqual(summary["counts"]["isotope"], 1)
        self.assertEqual(summary["representative_spectra_by_tier"]["LBM2"], 1)

    def test_no_identity_field_reaches_the_tool(self):
        _, records, manifest = self._export("unresolved")
        for record in records:
            self.assertFalse({"FORMULA", "INCHIKEY", "SMILES", "ONTOLOGY"} & set(record.fields))
            self.assertEqual(record.comment_tokens["UNIT"], "unit1")
        self.assertTrue(all(m["reference_inchikey"] == "" for m in manifest))

    def test_benchmark_keeps_the_reference_only_in_the_manifest(self):
        _, records, manifest = self._export("benchmark")
        self.assertEqual([r.fields["NAME"] for r in records], ["unit1_1"])
        self.assertNotIn("INCHIKEY", records[0].fields)
        self.assertEqual(manifest[0]["reference_inchikey"], "MUBZPKHOEPUJKR-UHFFFAOYSA-N")

    def test_cap_is_deterministic_and_recorded(self):
        first, records_a, _ = self._export("unresolved", max_spectra=2)
        second, records_b, _ = self._export("unresolved", max_spectra=2)
        self.assertEqual([r.fields["NAME"] for r in records_a], [r.fields["NAME"] for r in records_b])
        self.assertEqual(first["counts"]["selected_before_cap"], 3)
        self.assertEqual(first["counts"]["selected"], 2)
        self.assertEqual(first["msp_sha256"], second["msp_sha256"])

    def test_no_library_path_in_the_selection_record(self):
        summary, _, _ = self._export("unresolved")
        self.assertNotIn("private", json.dumps(summary))


if __name__ == "__main__":
    unittest.main()
