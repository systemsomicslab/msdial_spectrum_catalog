import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msdial_spectrum_catalog.devset import build_devset, ledger_state, read_done_units


class DevsetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        folder = self.root / "campaign"
        folder.mkdir()
        self.ledger = folder / "ledger.sqlite"
        self.writer = sqlite3.connect(self.ledger)
        self.writer.executescript("""
            CREATE TABLE unit(unit_key TEXT, state TEXT, workspace TEXT,
              catalog_unit_id TEXT, parent_unit_key TEXT, role TEXT,
              repository TEXT, accession TEXT, ion_mode TEXT, instrument TEXT, handoff_path TEXT);
            CREATE TABLE transition(unit_key TEXT, to_state TEXT, seq INTEGER);
        """)
        for n, key in enumerate(("part", "excluded", "pending")):
            unit = self.root / key
            out = unit / "output"
            out.mkdir(parents=True)
            (out / "AlignResult.mdmsp").write_text(
                "NAME: Invented\nNum Peaks: 1\n50 1\n\nNAME: Empty\nNum Peaks: 0\n", encoding="utf-8")
            (out / "sample.mdpeak").touch()
            (out / "analysis_files.csv").write_text("acquisition_type\nSWATH\n", encoding="utf-8")
            state = "running" if key == "pending" else "done"
            self.writer.execute("INSERT INTO unit VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                (key, state, str(unit), "parent", "parent", "part",
                                 "synthetic", "S1", "Both", "Invented", None))
            self.writer.execute("INSERT INTO transition VALUES (?, ?, ?)", (key, state, n))
        self.writer.commit()
        self.addCleanup(self.writer.close)

    def build(self, filename="set.json", **kwargs):
        with patch("msdial_spectrum_catalog.devset.git_commit", return_value="synthetic"):
            return build_devset(self.root, ["campaign"], self.root / filename, **kwargs)

    def test_split_absent_from_manifest_exclusion_flags_counts(self):
        result = self.build(exclusions={"excluded": "synthetic exclusion"})
        self.assertEqual(result["totals"], dict(units=1, samples=1, representative_spectra_nonempty=1))
        unit = result["units"][0]
        self.assertEqual(unit["unit_key"], "part")
        self.assertEqual(unit["parent_unit_key"], "parent")
        self.assertTrue(unit["dia"])
        self.assertTrue(unit["ion_mode_both"])
        self.assertEqual(unit["alignment_mdmsp"]["records"], 2)
        self.assertEqual(result["campaigns"][0]["done_units"], 2)

    def test_content_digest_determinism_and_no_overwrite(self):
        first = self.build()
        second = self.build("again.json")
        self.assertEqual(first["content_sha256"], second["content_sha256"])
        with self.assertRaises(FileExistsError):
            self.build()
        self.assertEqual(json.loads((self.root / "set.json").read_text())["content_sha256"],
                         first["content_sha256"])

    def test_live_wal_requires_snapshot_and_snapshot_sees_done(self):
        self.writer.execute("PRAGMA journal_mode=WAL")
        self.writer.execute("INSERT INTO transition VALUES ('part', 'done', 100)")
        self.writer.commit()
        self.assertGreater(ledger_state(self.ledger)["wal_bytes"], 0)
        with self.assertRaisesRegex(ValueError, "--snapshot"):
            self.build()
        before = ledger_state(self.ledger)
        result = self.build(snapshot=True, name="stage3b-snapshot")
        self.assertEqual(result["campaigns"][0]["read_mode"], "ro-snapshot")
        self.assertEqual(result["campaigns"][0]["ledger_state"], before)
        self.assertEqual(before, ledger_state(self.ledger))
        self.assertEqual({u["unit_key"] for u in result["units"]}, {"part", "excluded"})
        self.assertNotIn("subsets", result)

    def test_state_disagreement_is_error(self):
        self.writer.execute("UPDATE unit SET state='running' WHERE unit_key='part'")
        self.writer.commit()
        with self.assertRaisesRegex(ValueError, "disagrees"):
            read_done_units(self.ledger)

    def test_changed_snapshot_is_error(self):
        before = ledger_state(self.ledger)
        with patch("msdial_spectrum_catalog.devset.ledger_state",
                   side_effect=[before, dict(before, wal_bytes=123)]):
            with self.assertRaisesRegex(ValueError, "changed"):
                read_done_units(self.ledger, snapshot=True)

    def test_snapshot_never_falls_back_to_immutable(self):
        with patch("msdial_spectrum_catalog.devset.sqlite3.connect",
                   side_effect=sqlite3.OperationalError("synthetic")) as connect:
            with self.assertRaises(sqlite3.OperationalError):
                read_done_units(self.ledger, snapshot=True)
            self.assertEqual(connect.call_count, 1)
            self.assertTrue(connect.call_args.args[0].endswith("?mode=ro"))

    def test_refuses_analysis_output(self):
        with self.assertRaisesRegex(ValueError, "analysis root"):
            with patch("msdial_spectrum_catalog.devset.git_commit", return_value="synthetic"):
                build_devset(self.root, ["campaign"],
                             self.root / "MSDIAL_Public_Reanalysis/analysis/set.json")
