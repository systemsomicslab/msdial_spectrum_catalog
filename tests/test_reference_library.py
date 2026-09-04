from __future__ import annotations

import json
import tempfile
import unittest
import zlib
from pathlib import Path

from msdial_spectrum_catalog.reference_library import (
    ReferenceRecord,
    build_consensus,
    collision_energy_bin,
    inchikey_skeleton,
    ingest_reference_library,
    iter_library_consensus,
    normalize_collision_energy,
    normalize_instrument_type,
    normalize_ion_mode,
    reference_record_from_msp,
)
from msdial_spectrum_catalog.parsers import read_msp
from msdial_spectrum_catalog.storage import connect

SKELETON = "ABCDEFGHIJKLMN"
INCHIKEY = f"{SKELETON}-OPQRSTUVWX-Y"


def _msp(*records: dict) -> str:
    """Render MSP text; '_peaks' carries the peak list and '_num_peaks' overrides the declared count."""
    blocks = []
    for record in records:
        peaks = record.get("_peaks", [])
        lines = [f"{key}: {value}" for key, value in record.items() if not key.startswith("_")]
        lines.append(f"Num Peaks: {record.get('_num_peaks', len(peaks))}")
        lines.extend(f"{mz} {intensity}" for mz, intensity in peaks)
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks) + "\n"


def _record(name: str, **overrides) -> dict:
    record = {
        "NAME": name,
        "PRECURSORMZ": "300.1",
        "PRECURSORTYPE": "[M+H]+",
        "FORMULA": "C15H10O7",
        "INCHIKEY": INCHIKEY,
        "SMILES": "OC1=CC=CC=C1",
        "ONTOLOGY": "Flavonoids",
        "RETENTIONTIME": "2.5",
        "IONMODE": "Positive",
        "COLLISIONENERGY": "20",
        "INSTRUMENTTYPE": "LC-ESI-QTOF",
        "_peaks": [(100.0, 100), (200.0, 50)],
    }
    record.update(overrides)
    return record


def _reference(index: int, peaks: list[tuple[float, float]], **overrides) -> ReferenceRecord:
    values = {
        "record_name": f"compound {index}",
        "inchikey": INCHIKEY,
        "inchikey_skeleton": SKELETON,
        "formula": "C15H10O7",
        "precursor_mz": 300.1,
        "precursor_type": "[M+H]+",
        "ion_mode": "Positive",
        "instrument_class": "TOF",
        "collision_energy_value": 20.0,
    }
    values.update(overrides)
    return ReferenceRecord(record_index=index, peaks=peaks, **values)


def _mz_values(spectrum) -> list[float]:
    return [round(mz, 6) for mz, _ in spectrum.peaks]


class ConditionNormalizationTest(unittest.TestCase):
    def test_instrument_strings_map_to_analyzer_classes(self):
        cases = {
            "Orbitrap": "FT",
            "ESI-QFT": "FT",
            "LC-ESI-QFT": "FT",
            "Q-Exactive": "FT",
            "QTOF": "TOF",
            "LC-ESI-QTOF": "TOF",
            "TOF": "TOF",
            "IT": "IT",
            "LC-ESI-IT": "IT",
            "QQQ": "QQQ",
            "LC-ESI-QQ": "QQQ",
        }
        for raw, expected in cases.items():
            with self.subTest(raw=raw):
                self.assertEqual(normalize_instrument_type(raw), expected)

    def test_ion_trap_orbitrap_hybrid_resolves_to_ft(self):
        # The trap in substring matching: 'ITFT' contains 'IT', but detection happens in the FT cell.
        for raw in ("ITFT", "LC-ESI-ITFT", "lc-esi-itft"):
            with self.subTest(raw=raw):
                self.assertEqual(normalize_instrument_type(raw), "FT")

    def test_instrument_matching_is_case_insensitive(self):
        self.assertEqual(normalize_instrument_type("orbitrap"), "FT")
        self.assertEqual(normalize_instrument_type("lc-esi-qtof"), "TOF")

    def test_unstated_instrument_is_unknown_rather_than_guessed(self):
        for raw in (None, "", "   ", "nan", "NaN", "CI-B-EI", "ESI"):
            with self.subTest(raw=raw):
                self.assertEqual(normalize_instrument_type(raw), "UNKNOWN")

    def test_collision_energy_parsing(self):
        self.assertEqual(normalize_collision_energy("20"), (20.0, None))
        self.assertEqual(normalize_collision_energy("20.0"), (20.0, None))
        # HCD names a dissociation method, not a unit, so no unit is invented for it.
        self.assertEqual(normalize_collision_energy("45HCD"), (45.0, None))
        self.assertEqual(normalize_collision_energy("35 % (nominal)"), (35.0, "percent"))

    def test_missing_collision_energy_is_not_zero(self):
        for raw in ("nan", "NAN", "", "   ", None):
            with self.subTest(raw=raw):
                value, unit = normalize_collision_energy(raw)
                self.assertIsNone(value)
                self.assertIsNone(unit)
                self.assertNotEqual(value, 0.0)

    def test_collision_energy_units_are_recorded_only_when_stated(self):
        self.assertEqual(normalize_collision_energy("20 eV"), (20.0, "eV"))
        self.assertEqual(normalize_collision_energy("NCE 35"), (35.0, "NCE"))
        self.assertEqual(normalize_collision_energy("35%"), (35.0, "percent"))
        # 'level' contains 'ev'; a unit is a token, not a substring.
        self.assertEqual(normalize_collision_energy("30 (nominal level)"), (30.0, None))

    def test_collision_energy_ramp_uses_its_midpoint(self):
        self.assertEqual(normalize_collision_energy("Ramp 20-40 eV"), (30.0, "eV"))

    def test_collision_energy_bins(self):
        self.assertEqual(collision_energy_bin(20.0), "CE20-30")
        self.assertEqual(collision_energy_bin(25.4), "CE20-30")
        self.assertEqual(collision_energy_bin(30.0), "CE30-40")
        self.assertEqual(collision_energy_bin(0.0), "CE0-10")
        self.assertEqual(collision_energy_bin(45.0, 5.0), "CE45-50")
        self.assertEqual(collision_energy_bin(None), "CE_UNKNOWN")

    def test_collision_energy_bin_refuses_nonpositive_width(self):
        with self.assertRaises(ValueError):
            collision_energy_bin(20.0, 0.0)

    def test_inchikey_skeleton(self):
        self.assertEqual(inchikey_skeleton(INCHIKEY), SKELETON)
        self.assertEqual(inchikey_skeleton(INCHIKEY.lower()), SKELETON)
        self.assertEqual(inchikey_skeleton(SKELETON), SKELETON)
        self.assertIsNone(inchikey_skeleton(None))
        self.assertIsNone(inchikey_skeleton(""))
        self.assertIsNone(inchikey_skeleton("nan"))
        self.assertIsNone(inchikey_skeleton("NOTAKEY-XXXXXXXXXX-Y"))

    def test_ion_mode_normalization(self):
        self.assertEqual(normalize_ion_mode("POSITIVE"), "Positive")
        self.assertEqual(normalize_ion_mode("Positive"), "Positive")
        self.assertEqual(normalize_ion_mode("N"), "Negative")
        self.assertIsNone(normalize_ion_mode("nan"))
        self.assertIsNone(normalize_ion_mode("both"))

    def test_collision_energy_recovered_from_comment(self):
        text = _msp(_record("recovered", COLLISIONENERGY="nan", COMMENT="DB#=1|CE=32 eV|SPLASH=x"))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "library.msp"
            path.write_text(text, encoding="utf-8")
            (record,) = list(read_msp(path))
        reference = reference_record_from_msp(record)
        self.assertEqual(reference.collision_energy_raw, "32 eV")
        self.assertEqual(reference.collision_energy_value, 32.0)
        self.assertEqual(reference.collision_energy_unit, "eV")


class BuildConsensusTest(unittest.TestCase):
    def test_shared_skeleton_and_ce_bin_merge_while_other_bin_stays_separate(self):
        text = _msp(
            _record("quercetin A", COLLISIONENERGY="20"),
            _record("quercetin B", COLLISIONENERGY="25", _peaks=[(100.0, 90), (200.0, 60)]),
            _record("quercetin C", COLLISIONENERGY="45", _peaks=[(100.0, 10), (150.0, 100)]),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "library.msp"
            path.write_text(text, encoding="utf-8")
            records = [reference_record_from_msp(record) for record in read_msp(path)]
        spectra = list(build_consensus(records))
        self.assertEqual(len(spectra), 2)
        merged, separate = spectra
        self.assertEqual(merged.key.collision_energy_bin, "CE20-30")
        self.assertEqual(merged.key.instrument_class, "TOF")
        self.assertEqual(merged.key.inchikey_skeleton, SKELETON)
        self.assertEqual(merged.key.ion_mode, "Positive")
        self.assertEqual(merged.member_count, 2)
        self.assertEqual(merged.record_names, ["quercetin A", "quercetin B"])
        self.assertEqual(merged.formula, "C15H10O7")
        self.assertEqual(_mz_values(merged), [100.0, 200.0])
        self.assertAlmostEqual(merged.peaks[1][1], (0.5 + 60 / 90) / 2, places=6)
        self.assertEqual(separate.key.collision_energy_bin, "CE40-50")
        self.assertEqual(separate.member_count, 1)
        self.assertEqual(_mz_values(separate), [100.0, 150.0])

    def test_minimum_member_fraction_drops_a_peak_only_one_member_has(self):
        members = [
            _reference(0, [(100.0, 100), (150.0, 50), (200.0, 80)]),
            _reference(1, [(100.0, 100), (200.0, 80)]),
            _reference(2, [(100.0, 100), (200.0, 80)]),
        ]
        (spectrum,) = list(build_consensus(members))
        self.assertEqual(spectrum.member_count, 3)
        self.assertEqual(_mz_values(spectrum), [100.0, 200.0])

    def test_a_peak_two_of_three_members_have_is_kept(self):
        members = [
            _reference(0, [(100.0, 100), (150.0, 50)]),
            _reference(1, [(100.0, 100), (150.0, 40)]),
            _reference(2, [(100.0, 100)]),
        ]
        (spectrum,) = list(build_consensus(members))
        self.assertEqual(_mz_values(spectrum), [100.0, 150.0])

    def test_single_member_group_passes_through_unchanged(self):
        (spectrum,) = list(build_consensus([_reference(0, [(100.0, 10), (200.0, 20), (300.0, 5)])]))
        self.assertEqual(spectrum.member_count, 1)
        self.assertEqual(
            [(round(mz, 6), round(intensity, 6)) for mz, intensity in spectrum.peaks],
            [(100.0, 0.5), (200.0, 1.0), (300.0, 0.25)],
        )

    def test_records_without_a_skeleton_are_never_merged(self):
        members = [
            _reference(0, [(100.0, 100)], inchikey=None, inchikey_skeleton=None),
            _reference(1, [(100.0, 100)], inchikey=None, inchikey_skeleton=None),
        ]
        spectra = list(build_consensus(members))
        self.assertEqual([spectrum.member_count for spectrum in spectra], [1, 1])


class IngestReferenceLibraryTest(unittest.TestCase):
    def _write(self, directory: Path, text: str) -> Path:
        path = directory / "library.msp"
        path.write_text(text, encoding="utf-8")
        return path

    def test_ingest_writes_rows_normalized_conditions_and_consensus(self):
        text = _msp(
            _record("quercetin A", COLLISIONENERGY="20", INSTRUMENTTYPE="LC-ESI-ITFT"),
            _record("quercetin B", COLLISIONENERGY="25", INSTRUMENTTYPE="LC-ESI-ITFT",
                    _peaks=[(100.0, 90), (200.0, 60)]),
            _record("quercetin C", COLLISIONENERGY="nan", INSTRUMENTTYPE="LC-ESI-ITFT",
                    _peaks=[(100.0, 10), (150.0, 100)]),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            path = self._write(root, text)
            report = ingest_reference_library(
                database,
                path,
                library_name="demo",
                library_version="v1",
                source_uri="https://example.invalid/demo.msp",
                license="CC-BY-4.0",
            )
            self.assertTrue(report.valid)
            self.assertEqual(report.records_read, 3)
            self.assertEqual(report.records_skipped, 0)
            self.assertEqual(report.consensus_spectra, 2)
            connection = connect(database)
            try:
                library = connection.execute(
                    "SELECT * FROM reference_library WHERE library_id = ?", (report.library_id,)
                ).fetchone()
                self.assertEqual(library["library_kind"], "experimental_reference")
                self.assertEqual(library["record_count"], 3)
                self.assertEqual(library["byte_size"], path.stat().st_size)
                self.assertEqual(library["license"], "CC-BY-4.0")
                rows = connection.execute(
                    "SELECT * FROM reference_spectrum WHERE library_id = ? "
                    "ORDER BY library_record_index",
                    (report.library_id,),
                ).fetchall()
                self.assertEqual(len(rows), 3)
                self.assertEqual([row["instrument_class"] for row in rows], ["FT", "FT", "FT"])
                self.assertEqual([row["inchikey_skeleton"] for row in rows], [SKELETON] * 3)
                self.assertEqual([row["ion_mode"] for row in rows], ["Positive"] * 3)
                self.assertEqual([row["collision_energy_value"] for row in rows], [20.0, 25.0, None])
                self.assertEqual(rows[2]["collision_energy_raw"], "nan")
                self.assertEqual([row["peak_count"] for row in rows], [2, 2, 2])
            finally:
                connection.close()
            bins = sorted(
                spectrum.key.collision_energy_bin
                for spectrum in iter_library_consensus(database, report.library_id)
            )
            self.assertEqual(bins, ["CE20-30", "CE_UNKNOWN"])

    def test_consensus_payloads_are_stored_with_their_parameters(self):
        text = _msp(_record("quercetin A"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            report = ingest_reference_library(
                database, self._write(root, text), library_name="demo", library_version="v1"
            )
            payloads = _consensus_payloads(database)
            self.assertEqual(len(payloads), 1)
            payload = payloads[0]
            self.assertEqual(payload["library_id"], report.library_id)
            self.assertEqual(payload["library_kind"], "experimental_reference")
            self.assertEqual(payload["consensus_key"]["inchikey_skeleton"], SKELETON)
            self.assertEqual(payload["consensus_key"]["collision_energy_bin"], "CE20-30")
            self.assertEqual(payload["parameters"]["minimum_member_fraction"], 0.5)
            self.assertEqual(payload["parameters"]["mz_bin_width"], 0.01)

    def test_limit_is_honoured(self):
        text = _msp(*[_record(f"compound {index}") for index in range(5)])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            report = ingest_reference_library(
                database, self._write(root, text), library_name="demo", limit=2
            )
            self.assertEqual(report.records_read, 2)
            connection = connect(database)
            try:
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM reference_spectrum").fetchone()[0], 2
                )
                self.assertEqual(
                    connection.execute(
                        "SELECT record_count FROM reference_library WHERE library_id = ?",
                        (report.library_id,),
                    ).fetchone()[0],
                    2,
                )
            finally:
                connection.close()

    def test_precursor_mz_range_is_honoured(self):
        text = _msp(
            _record("in range", PRECURSORMZ="300.1"),
            _record("above range", PRECURSORMZ="900.4"),
            _record("below range", PRECURSORMZ="100.2"),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            report = ingest_reference_library(
                database,
                self._write(root, text),
                library_name="demo",
                precursor_mz_range=(200.0, 400.0),
            )
            self.assertEqual(report.records_read, 1)
            self.assertEqual(report.records_skipped, 0)
            connection = connect(database)
            try:
                names = [
                    row[0]
                    for row in connection.execute("SELECT record_name FROM reference_spectrum")
                ]
            finally:
                connection.close()
            self.assertEqual(names, ["in range"])

    def test_malformed_records_are_counted_rather_than_raised(self):
        text = _msp(
            _record("usable"),
            _record("bad count", _num_peaks="abc"),
            _record("no peaks", _peaks=[]),
            _record("no precursor", PRECURSORMZ=""),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            report = ingest_reference_library(database, self._write(root, text), library_name="demo")
            self.assertTrue(report.valid)
            self.assertEqual(report.records_read, 4)
            self.assertEqual(report.records_skipped, 3)
            connection = connect(database)
            try:
                names = [
                    row[0]
                    for row in connection.execute("SELECT record_name FROM reference_spectrum")
                ]
            finally:
                connection.close()
            self.assertEqual(names, ["usable"])

    def test_reingesting_the_same_file_does_not_duplicate_rows(self):
        text = _msp(_record("quercetin A"), _record("quercetin B", _peaks=[(100.0, 90), (200.0, 60)]))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            path = self._write(root, text)
            first = ingest_reference_library(database, path, library_name="demo", library_version="v1")
            second = ingest_reference_library(database, path, library_name="demo", library_version="v1")
            self.assertEqual(first.library_id, second.library_id)
            self.assertEqual(second.records_read, 2)
            self.assertEqual(second.blobs_written, 0)
            self.assertEqual(second.consensus_spectra, first.consensus_spectra)
            self.assertEqual(second.warnings, [])
            connection = connect(database)
            try:
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM reference_library").fetchone()[0], 1
                )
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM reference_spectrum").fetchone()[0], 2
                )
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM spectrum_blob").fetchone()[0],
                    first.blobs_written,
                )
            finally:
                connection.close()

    def test_byte_identical_records_share_one_blob(self):
        text = _msp(_record("quercetin"), _record("quercetin"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            report = ingest_reference_library(
                database, self._write(root, text), library_name="demo", consensus=False
            )
            self.assertEqual(report.records_read, 2)
            self.assertEqual(report.blobs_written, 1)
            connection = connect(database)
            try:
                digests = {
                    row[0]
                    for row in connection.execute("SELECT payload_sha256 FROM reference_spectrum")
                }
                self.assertEqual(len(digests), 1)
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM spectrum_blob").fetchone()[0], 1
                )
            finally:
                connection.close()

    def test_records_differing_outside_the_peak_list_keep_separate_blobs(self):
        # The payload keeps the record's verbatim fields, matching ingest._put_blob, so two records with
        # the same peaks but different names are two payloads. Pinned here because it is the visible
        # consequence of sharing one blob store with the MS-DIAL ingest path.
        text = _msp(_record("quercetin"), _record("kaempferol"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            report = ingest_reference_library(
                database, self._write(root, text), library_name="demo", consensus=False
            )
            self.assertEqual(report.blobs_written, 2)

    def test_changed_file_under_the_same_version_warns(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            path = self._write(root, _msp(_record("quercetin")))
            ingest_reference_library(database, path, library_name="demo", library_version="v1")
            path.write_text(_msp(_record("kaempferol")), encoding="utf-8")
            report = ingest_reference_library(database, path, library_name="demo", library_version="v1")
            self.assertTrue(report.valid)
            self.assertEqual(len(report.warnings), 1)
            self.assertIn("digest", report.warnings[0])

    def test_in_silico_library_kind_is_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "catalog.sqlite"
            report = ingest_reference_library(
                database,
                self._write(root, _msp(_record("predicted"))),
                library_name="predictions",
                library_kind="in_silico_predicted",
            )
            connection = connect(database)
            try:
                kind = connection.execute(
                    "SELECT library_kind FROM reference_library WHERE library_id = ?",
                    (report.library_id,),
                ).fetchone()[0]
            finally:
                connection.close()
            self.assertEqual(kind, "in_silico_predicted")

    def test_unknown_library_kind_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(ValueError):
                ingest_reference_library(
                    root / "catalog.sqlite",
                    self._write(root, _msp(_record("quercetin"))),
                    library_name="demo",
                    library_kind="guessed",
                )


def _consensus_payloads(database: Path) -> list[dict]:
    connection = connect(database)
    try:
        payloads = []
        for row in connection.execute("SELECT payload FROM spectrum_blob"):
            payload = json.loads(zlib.decompress(row[0]).decode("utf-8"))
            if payload.get("kind") == "skeleton_consensus":
                payloads.append(payload)
        return payloads
    finally:
        connection.close()


if __name__ == "__main__":
    unittest.main()
