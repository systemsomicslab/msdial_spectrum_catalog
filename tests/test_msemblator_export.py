"""Tests for the MS-Emblator export: unique names, no leaked identity, and the selection rules."""

from __future__ import annotations

import csv
import hashlib
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
        suffix = "_manifest.local-vs21.tsv" if selection == "benchmark" else "_manifest.tsv"
        self.assertEqual(summary["manifest"], f"unit1_{selection}{suffix}")
        with (destination / f"unit1_{selection}{suffix}").open(encoding="utf-8") as handle:
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

    def _noisy_spot_3(self, spot=3):
        """Give a spot (3 by default) three strong ions and 300 weak ones (spectral entropy above 4.5)."""
        text = (self.output / f"{STEM}.mdmsp").read_text(encoding="utf-8")
        noisy = _record(spot, 0).replace("Num Peaks: 0\n", "Num Peaks: 303\n").rstrip("\n") + "\n"
        noisy += "".join(f"{100.5 + k:.1f}\t{1000 + 50 * (k % 41)}\n" for k in range(300))
        noisy += "500.0\t100000\n600.0\t60000\n700.0\t30000\n\n"
        (self.output / f"{STEM}.mdmsp").write_text(text.replace(_record(spot, 5), noisy), encoding="utf-8")

    def test_entropy_cleanup_drops_weak_ions_and_records_it(self):
        self._noisy_spot_3()
        destination = self.root / "out"
        summary = export_unit(self.output, "unit1", destination, "unresolved", max_entropy=4.5)
        records = {r.fields["NAME"]: r for r in read_msp(destination / "unit1_unresolved.msp")}
        with (destination / "unit1_unresolved_manifest.tsv").open(encoding="utf-8") as handle:
            manifest = {m["name"]: m for m in csv.DictReader(handle, delimiter="\t")}
        cleaned, untouched = manifest["unit1_3"], manifest["unit1_4"]
        self.assertEqual(int(cleaned["raw_peak_count"]), 303)
        self.assertEqual(int(cleaned["peak_count"]), len(records["unit1_3"].peaks))
        self.assertLess(int(cleaned["peak_count"]), 303)
        self.assertIn((500.0, 100000.0), [(mz, float(i)) for mz, i in records["unit1_3"].peaks])
        self.assertGreater(float(cleaned["spectral_entropy_raw"]), 4.5)
        self.assertLessEqual(float(cleaned["spectral_entropy"]), 4.5)
        self.assertEqual(untouched["raw_peak_count"], untouched["peak_count"])
        self.assertEqual(untouched["spectral_entropy_raw"], untouched["spectral_entropy"])
        cleanup = summary["peak_cleanup"]
        self.assertEqual((cleanup["rule"], cleanup["max_entropy"]), ("spectral_entropy", 4.5))
        self.assertEqual(cleanup["spectra_cleaned"], 1)
        self.assertEqual(cleanup["peaks_removed"], 303 - int(cleaned["peak_count"]))
        self.assertEqual((cleaned["cleanup_status"], untouched["cleanup_status"]), ("cleaned", "unchanged"))
        kept = sorted(float(i) for _, i in records["unit1_3"].peaks)
        noise = [1000 + 50 * (k % 41) for k in range(300)]
        self.assertEqual(float(cleaned["kept_threshold_pct_of_base"]), round(100 * kept[0] / 100000, 2))
        self.assertEqual(float(cleaned["largest_removed_pct_of_base"]),
                         round(100 * max(i for i in noise if i < kept[0]) / 100000, 2))
        self.assertEqual(untouched["kept_threshold_pct_of_base"], "")
        # the hash describes the written peaks, not the raw ones
        written = "\n".join(f"{mz:.5f}\t{float(i):.0f}" for mz, i in records["unit1_3"].peaks)
        self.assertEqual(cleaned["peaks_sha256"], hashlib.sha256(written.encode()).hexdigest())

    def test_no_cleanup_by_default(self):
        self._noisy_spot_3()
        summary, records, manifest = self._export("unresolved")
        self.assertEqual(summary["peak_cleanup"], {"rule": None, "spectra_cleaned": 0, "peaks_removed": 0,
                                                   "spectra_above_ceiling": 0, "spectra_kept_raw_too_few_peaks": 0,
                                                   "non_positive_entries_written": 0})
        self.assertEqual(len(records[0].peaks), 303)
        self.assertEqual(manifest[0]["raw_peak_count"], manifest[0]["peak_count"])
        self.assertEqual(manifest[0]["cleanup_status"], "not_requested")

    def test_the_cleanup_never_removes_a_spectrum(self):
        # two strong ions over 200 equal weak ones: the weak level goes as a whole, leaving two ions
        text = (self.output / f"{STEM}.mdmsp").read_text(encoding="utf-8")
        record = _record(3, 0).replace("Num Peaks: 0\n", "Num Peaks: 202\n").rstrip("\n") + "\n"
        record += "".join(f"{100.5 + k:.1f}\t10\n" for k in range(200)) + "500.0\t100\n600.0\t50\n\n"
        (self.output / f"{STEM}.mdmsp").write_text(text.replace(_record(3, 5), record), encoding="utf-8")
        summary = export_unit(self.output, "unit1", self.root / "out", "unresolved", max_entropy=4.5)
        with (self.root / "out" / "unit1_unresolved_manifest.tsv").open(encoding="utf-8") as handle:
            manifest = {m["name"]: m for m in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(sorted(manifest), ["unit1_3", "unit1_4", "unit1_7"])  # the same spectra as without it
        self.assertEqual(manifest["unit1_3"]["cleanup_status"], "kept_raw_too_few_peaks")
        self.assertEqual(manifest["unit1_3"]["peak_count"], "202")
        cleanup = summary["peak_cleanup"]
        self.assertEqual((cleanup["spectra_kept_raw_too_few_peaks"], cleanup["spectra_cleaned"],
                          cleanup["peaks_removed"]), (1, 0, 0))
        self.assertEqual(sum(summary["representative_spectra_by_tier"].values()) + summary["counts"]["too_few_peaks"]
                         + summary["counts"]["isotope"] + summary["counts"]["no_peak_id"], summary["counts"]["records"])

    def test_zero_entries_are_written_in_both_arms(self):
        text = (self.output / f"{STEM}.mdmsp").read_text(encoding="utf-8")
        self._noisy_spot_3()
        noisy = (self.output / f"{STEM}.mdmsp").read_text(encoding="utf-8")
        with_zero = noisy.replace("Num Peaks: 303\n", "Num Peaks: 304\n").replace("500.0\t100000\n",
                                                                                  "450.0\t0\n500.0\t100000\n")
        with_zero = with_zero.replace(_record(4, 5), _record(4, 5).replace("Num Peaks: 5\n", "Num Peaks: 6\n")
                                      .replace("100.0\t10\n", "99.0\t0\n100.0\t10\n"))
        self.assertNotEqual(text, with_zero)
        (self.output / f"{STEM}.mdmsp").write_text(with_zero, encoding="utf-8")
        plain = export_unit(self.output, "unit1", self.root / "a", "unresolved")
        cleaned = export_unit(self.output, "unit1", self.root / "b", "unresolved", max_entropy=4.5)
        self.assertEqual(plain["peak_cleanup"]["non_positive_entries_written"], 2)
        self.assertEqual(cleaned["peak_cleanup"]["non_positive_entries_written"], 2)
        records = {r.fields["NAME"]: r for r in read_msp(self.root / "b" / "unit1_unresolved.msp")}
        self.assertIn(450.0, [mz for mz, _ in records["unit1_3"].peaks])

    def test_counts_describe_the_written_spectra_after_the_cap(self):
        self._noisy_spot_3(spot=4)  # a cap of 1 keeps spot 3 (lowest sampling rank), so the cleaned spot 4 goes
        for cap, expected in ((1, (0, 0)), (3, None)):
            destination = self.root / f"cap{cap}"
            summary = export_unit(self.output, "unit1", destination, "unresolved", max_spectra=cap, max_entropy=4.5)
            with (destination / "unit1_unresolved_manifest.tsv").open(encoding="utf-8") as handle:
                manifest = list(csv.DictReader(handle, delimiter="\t"))
            removed = sum(int(m["raw_peak_count"]) - int(m["peak_count"]) for m in manifest)
            cleanup = summary["peak_cleanup"]
            self.assertEqual(cleanup["peaks_removed"], removed)
            self.assertEqual(cleanup["spectra_cleaned"], sum(m["cleanup_status"] == "cleaned" for m in manifest))
            if expected is not None:
                self.assertEqual([m["name"] for m in manifest], ["unit1_3"])
                self.assertEqual((cleanup["spectra_cleaned"], cleanup["peaks_removed"]), expected)
            else:
                self.assertEqual(cleanup["spectra_cleaned"], 1)
                self.assertGreater(cleanup["peaks_removed"], 0)

    def test_the_too_few_peaks_fallback_counts_positive_ions(self):
        # two strong ions, 200 equal weak ones and 5 zero entries: 7 entries would remain, but 2 positive ions
        text = (self.output / f"{STEM}.mdmsp").read_text(encoding="utf-8")
        record = _record(3, 0).replace("Num Peaks: 0\n", "Num Peaks: 207\n").rstrip("\n") + "\n"
        record += "".join(f"{100.5 + k:.1f}\t10\n" for k in range(200))
        record += "".join(f"{400.5 + k:.1f}\t0\n" for k in range(5)) + "500.0\t100\n600.0\t50\n\n"
        (self.output / f"{STEM}.mdmsp").write_text(text.replace(_record(3, 5), record), encoding="utf-8")
        summary = export_unit(self.output, "unit1", self.root / "out", "unresolved", max_entropy=4.5)
        with (self.root / "out" / "unit1_unresolved_manifest.tsv").open(encoding="utf-8") as handle:
            manifest = {m["name"]: m for m in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(manifest["unit1_3"]["cleanup_status"], "kept_raw_too_few_peaks")
        self.assertEqual(manifest["unit1_3"]["peak_count"], "207")
        self.assertEqual(summary["peak_cleanup"]["non_positive_entries_written"], 5)

    def test_a_spectrum_above_the_ceiling_is_counted(self):
        text = (self.output / f"{STEM}.mdmsp").read_text(encoding="utf-8")
        record = _record(3, 0).replace("Num Peaks: 0\n", "Num Peaks: 201\n").rstrip("\n") + "\n"
        record += "".join(f"{100.5 + k:.1f}\t10\n" for k in range(200)) + "500.0\t1\n\n"
        (self.output / f"{STEM}.mdmsp").write_text(text.replace(_record(3, 5), record), encoding="utf-8")
        summary = export_unit(self.output, "unit1", self.root / "out", "unresolved", max_entropy=4.5)
        with (self.root / "out" / "unit1_unresolved_manifest.tsv").open(encoding="utf-8") as handle:
            manifest = {m["name"]: m for m in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(manifest["unit1_3"]["cleanup_status"], "above_ceiling")
        self.assertGreater(float(manifest["unit1_3"]["spectral_entropy"]), 4.5)
        cleanup = summary["peak_cleanup"]
        self.assertEqual((cleanup["spectra_above_ceiling"], cleanup["spectra_cleaned"], cleanup["peaks_removed"]),
                         (1, 1, 1))

    def test_a_ceiling_that_is_not_finite_and_positive_is_refused(self):
        for ceiling in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError, msg=repr(ceiling)):
                export_unit(self.output, "unit1", self.root / "out", "unresolved", max_entropy=ceiling)


if __name__ == "__main__":
    unittest.main()
