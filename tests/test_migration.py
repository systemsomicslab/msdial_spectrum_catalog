import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog import storage
from msdial_spectrum_catalog.schema import MIGRATIONS, SCHEMA_VERSION

NEW_TABLES = (
    "schema_migration",
    "annotation_tool_run",
    "reference_library",
    "reference_spectrum",
    "spectrum_similarity",
    "ambiguity_class",
    "ambiguity_class_member",
    "annotation_candidate",
    "annotation_claim_component",
    "criteria_rule",
    "msdial_annotation_result",
)

NEW_INDEXES = (
    "idx_annotation_assertion_spectrum",
    "idx_annotation_evidence_assertion",
    "idx_annotation_candidate_assertion",
    "idx_reference_spectrum_library",
    "idx_reference_spectrum_skeleton",
    "idx_ambiguity_member_class",
    "idx_msdial_annotation_subject",
    "idx_criteria_rule_set",
)

EXPECTED_COLUMNS = {
    "artifact": ("source_path",),
    "sample": ("raw_file_path",),
    "criteria_set": ("vocab_version",),
    "annotation_assertion": (
        "subject_kind",
        "alignment_feature_id",
        "vocab_version",
        "claim_concept_id",
        "notation_verbatim",
        "claim_unresolved",
        "candidate_count",
        "ambiguity_class_id",
        "tool_run_id",
    ),
    "annotation_evidence": (
        "evidence_concept_id",
        "evidence_subtype",
        "metric",
        "measured_value",
        "measured_unit",
        "comparison",
        "threshold_value",
        "passed",
        "criteria_rule_id",
        "source_spectrum_id",
        "source_reference_spectrum_id",
        "source_tool_run_id",
        "out_of_distribution",
        "ood_reason",
    ),
    "msdial_annotation_result": (
        "annotation_kind",
        "candidate_name",
        "candidate_is_named",
    ),
}

LEGACY_SQL = """
CREATE TABLE catalog_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

CREATE TABLE study (
    study_id TEXT PRIMARY KEY,
    repository TEXT NOT NULL,
    accession TEXT NOT NULL,
    title TEXT,
    UNIQUE(repository, accession)
);

CREATE TABLE analysis_unit (
    analysis_unit_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    external_unit_id TEXT NOT NULL,
    separation_type TEXT,
    ion_mode TEXT,
    acquisition_type TEXT,
    UNIQUE(study_id, external_unit_id)
);

CREATE TABLE analysis_run (
    run_id TEXT PRIMARY KEY,
    analysis_unit_id TEXT NOT NULL REFERENCES analysis_unit(analysis_unit_id),
    run_fingerprint TEXT NOT NULL,
    output_directory TEXT NOT NULL,
    msdial_version TEXT,
    interactive_version TEXT,
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE artifact (
    artifact_id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    artifact_type TEXT NOT NULL,
    relative_path TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    byte_size INTEGER NOT NULL,
    UNIQUE(run_id, relative_path)
);

CREATE TABLE sample (
    sample_id TEXT PRIMARY KEY,
    analysis_unit_id TEXT NOT NULL REFERENCES analysis_unit(analysis_unit_id),
    sample_name TEXT NOT NULL,
    raw_file_name TEXT,
    repository_sample_id TEXT,
    UNIQUE(analysis_unit_id, sample_name)
);

CREATE TABLE spectrum_blob (
    payload_sha256 TEXT PRIMARY KEY,
    compression TEXT NOT NULL,
    uncompressed_bytes INTEGER NOT NULL,
    payload BLOB NOT NULL
);

CREATE TABLE spectrum (
    spectrum_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    sample_id TEXT REFERENCES sample(sample_id),
    feature_id TEXT,
    alignment_feature_id TEXT,
    spectrum_kind TEXT NOT NULL,
    source_peak_id INTEGER,
    ms1_scan_index INTEGER,
    ms2_scan_index INTEGER,
    native_id TEXT,
    usi TEXT,
    precursor_mz REAL,
    rt_min REAL,
    ion_mode TEXT,
    peak_count INTEGER NOT NULL,
    payload_sha256 TEXT NOT NULL REFERENCES spectrum_blob(payload_sha256),
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_record INTEGER
);

CREATE TABLE criteria_set (
    criteria_set_id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    version TEXT,
    description TEXT,
    rules_json TEXT NOT NULL
);

CREATE TABLE annotation_assertion (
    assertion_id TEXT PRIMARY KEY,
    spectrum_id TEXT NOT NULL REFERENCES spectrum(spectrum_id),
    annotation_level TEXT NOT NULL,
    annotation_claim TEXT,
    compound_name TEXT,
    formula TEXT,
    structure_id TEXT,
    criteria_set_id TEXT REFERENCES criteria_set(criteria_set_id),
    curation_comment TEXT,
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE annotation_evidence (
    evidence_id TEXT PRIMARY KEY,
    assertion_id TEXT NOT NULL REFERENCES annotation_assertion(assertion_id),
    evidence_tag TEXT NOT NULL,
    evidence_value_json TEXT,
    source_uri TEXT
);
"""


def table_columns(connection, table):
    return {row[1] for row in connection.execute(f"PRAGMA table_info({table})")}


def object_names(connection, kind):
    return {
        row[0]
        for row in connection.execute("SELECT name FROM sqlite_master WHERE type = ?", (kind,))
    }


def build_legacy_database(path):
    connection = sqlite3.connect(path)
    try:
        connection.executescript(LEGACY_SQL)
        connection.execute(
            "INSERT INTO catalog_meta(key, value) VALUES ('schema_version', '2')"
        )
        connection.execute(
            "INSERT INTO study(study_id, repository, accession, title) VALUES (?, ?, ?, ?)",
            ("urn:msdial:study:mb_post:MPST000007", "mb_post", "MPST000007", "Legacy study"),
        )
        connection.execute(
            "INSERT INTO analysis_unit(analysis_unit_id, study_id, external_unit_id, ion_mode) "
            "VALUES (?, ?, ?, ?)",
            (
                "urn:msdial:analysis_unit:legacy",
                "urn:msdial:study:mb_post:MPST000007",
                "lcms_neg_dda",
                "Negative",
            ),
        )
        connection.execute(
            "INSERT INTO analysis_run(run_id, analysis_unit_id, run_fingerprint, output_directory) "
            "VALUES (?, ?, ?, ?)",
            (
                "urn:msdial:analysis_run:legacy",
                "urn:msdial:analysis_unit:legacy",
                "deadbeef",
                "D:\\legacy\\output",
            ),
        )
        connection.execute(
            "INSERT INTO artifact(run_id, artifact_type, relative_path, sha256, byte_size) "
            "VALUES (?, ?, ?, ?, ?)",
            ("urn:msdial:analysis_run:legacy", "alignment", "AlignResult.mdalign", "abc123", 42),
        )
        connection.execute(
            "INSERT INTO sample(sample_id, analysis_unit_id, sample_name, raw_file_name) "
            "VALUES (?, ?, ?, ?)",
            (
                "urn:msdial:sample:legacy:sample_a",
                "urn:msdial:analysis_unit:legacy",
                "sample_a",
                "sample_a.raw",
            ),
        )
        connection.execute(
            "INSERT INTO spectrum_blob(payload_sha256, compression, uncompressed_bytes, payload) "
            "VALUES (?, ?, ?, ?)",
            ("f00d", "none", 4, b"peak"),
        )
        connection.execute(
            "INSERT INTO spectrum(spectrum_id, run_id, sample_id, spectrum_kind, peak_count, "
            "payload_sha256) VALUES (?, ?, ?, ?, ?, ?)",
            (
                "urn:msdial:spectrum:legacy:sample_a:deconvoluted:7",
                "urn:msdial:analysis_run:legacy",
                "urn:msdial:sample:legacy:sample_a",
                "deconvoluted",
                1,
                "f00d",
            ),
        )
        connection.execute(
            "INSERT INTO annotation_assertion(assertion_id, spectrum_id, annotation_level, "
            "compound_name) VALUES (?, ?, ?, ?)",
            (
                "urn:msdial:annotation_assertion:legacy:1",
                "urn:msdial:spectrum:legacy:sample_a:deconvoluted:7",
                "L2",
                "Legacy candidate",
            ),
        )
        connection.execute(
            "INSERT INTO annotation_evidence(evidence_id, assertion_id, evidence_tag) "
            "VALUES (?, ?, ?)",
            (
                "urn:msdial:annotation_evidence:legacy:1",
                "urn:msdial:annotation_assertion:legacy:1",
                "SL",
            ),
        )
        connection.commit()
    finally:
        connection.close()


class MigrationTableTests(unittest.TestCase):
    def test_migrations_are_ordered_and_reach_schema_version(self):
        versions = [migration.version for migration in MIGRATIONS]
        self.assertEqual(versions, sorted(set(versions)))
        self.assertEqual(versions[-1], SCHEMA_VERSION)


class FreshDatabaseTests(unittest.TestCase):
    def test_fresh_initialize_creates_current_tables_columns_and_indexes(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            connection = storage.connect(database)
            try:
                self.assertEqual(storage.read_schema_version(connection), SCHEMA_VERSION)
                self.assertEqual(MIGRATIONS[-1].version, SCHEMA_VERSION)
                tables = object_names(connection, "table")
                for table in NEW_TABLES:
                    self.assertIn(table, tables)
                indexes = object_names(connection, "index")
                for index in NEW_INDEXES:
                    self.assertIn(index, indexes)
                for table, columns in EXPECTED_COLUMNS.items():
                    present = table_columns(connection, table)
                    for column in columns:
                        self.assertIn(column, present, f"{table}.{column}")
                recorded = {
                    row[0] for row in connection.execute("SELECT version FROM schema_migration")
                }
                self.assertEqual(recorded, {m.version for m in MIGRATIONS})
            finally:
                connection.close()

    def test_initialize_twice_is_a_no_op(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            connection = storage.connect(database)
            try:
                before = sorted(
                    tuple(row) for row in connection.execute(
                        "SELECT version, description FROM schema_migration"
                    )
                )
            finally:
                connection.close()

            storage.initialize(database)

            connection = storage.connect(database)
            try:
                after = sorted(
                    tuple(row) for row in connection.execute(
                        "SELECT version, description FROM schema_migration"
                    )
                )
                self.assertEqual(before, after)
                self.assertEqual(storage.read_schema_version(connection), SCHEMA_VERSION)
                self.assertIn("subject_kind", table_columns(connection, "annotation_assertion"))
            finally:
                connection.close()

    def test_apply_migrations_reports_nothing_for_a_current_database(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            connection = storage.connect(database)
            try:
                self.assertEqual(storage.apply_migrations(connection, from_version=SCHEMA_VERSION), ())
                self.assertEqual(
                    storage.apply_migrations(connection, from_version=None),
                    tuple(m.version for m in MIGRATIONS),
                )
            finally:
                connection.close()


class LegacyUpgradeTests(unittest.TestCase):
    def test_version_2_database_converges_and_keeps_its_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "legacy.sqlite"
            build_legacy_database(database)

            connection = sqlite3.connect(database)
            try:
                self.assertNotIn("source_path", table_columns(connection, "artifact"))
                self.assertNotIn("raw_file_path", table_columns(connection, "sample"))
                self.assertNotIn(
                    "subject_kind", table_columns(connection, "annotation_assertion")
                )
            finally:
                connection.close()

            storage.initialize(database)

            connection = storage.connect(database)
            try:
                self.assertEqual(storage.read_schema_version(connection), SCHEMA_VERSION)
                for table, columns in EXPECTED_COLUMNS.items():
                    present = table_columns(connection, table)
                    for column in columns:
                        self.assertIn(column, present, f"{table}.{column}")
                tables = object_names(connection, "table")
                for table in NEW_TABLES:
                    self.assertIn(table, tables)
                recorded = {
                    row[0] for row in connection.execute("SELECT version FROM schema_migration")
                }
                self.assertEqual(recorded, {m.version for m in MIGRATIONS})

                self.assertEqual(
                    connection.execute("SELECT accession, title FROM study").fetchone()[:],
                    ("MPST000007", "Legacy study"),
                )
                self.assertEqual(
                    connection.execute(
                        "SELECT sample_name, raw_file_name, raw_file_path FROM sample"
                    ).fetchone()[:],
                    ("sample_a", "sample_a.raw", None),
                )
                self.assertEqual(
                    connection.execute(
                        "SELECT relative_path, sha256, byte_size, source_path FROM artifact"
                    ).fetchone()[:],
                    ("AlignResult.mdalign", "abc123", 42, ""),
                )
                self.assertEqual(
                    connection.execute(
                        "SELECT annotation_level, compound_name, subject_kind, claim_unresolved, "
                        "candidate_count FROM annotation_assertion"
                    ).fetchone()[:],
                    ("L2", "Legacy candidate", "spectrum", 0, 0),
                )
                self.assertEqual(
                    connection.execute(
                        "SELECT evidence_tag, measured_value, threshold_value, passed, "
                        "out_of_distribution FROM annotation_evidence"
                    ).fetchone()[:],
                    ("SL", None, None, None, 0),
                )
            finally:
                connection.close()


class VersionMarkerTests(unittest.TestCase):
    def test_newer_marker_raises(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            connection = sqlite3.connect(database)
            try:
                connection.execute(
                    "INSERT OR REPLACE INTO catalog_meta(key, value) VALUES ('schema_version', '99')"
                )
                connection.commit()
            finally:
                connection.close()

            with self.assertRaises(storage.SchemaVersionError):
                storage.initialize(database)

    def test_missing_marker_converges_to_current_version(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            connection = sqlite3.connect(database)
            try:
                connection.execute("DELETE FROM catalog_meta WHERE key = 'schema_version'")
                connection.commit()
            finally:
                connection.close()

            connection = storage.connect(database)
            try:
                self.assertIsNone(storage.read_schema_version(connection))
            finally:
                connection.close()

            storage.initialize(database)

            connection = storage.connect(database)
            try:
                self.assertEqual(storage.read_schema_version(connection), SCHEMA_VERSION)
                self.assertIn("vocab_version", table_columns(connection, "criteria_set"))
            finally:
                connection.close()

    def test_unparseable_marker_is_treated_as_unknown(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            connection = sqlite3.connect(database)
            try:
                connection.execute(
                    "INSERT OR REPLACE INTO catalog_meta(key, value) VALUES ('schema_version', 'v3')"
                )
                connection.commit()
            finally:
                connection.close()

            storage.initialize(database)

            connection = storage.connect(database)
            try:
                self.assertEqual(storage.read_schema_version(connection), SCHEMA_VERSION)
            finally:
                connection.close()


class ForeignKeyTests(unittest.TestCase):
    def test_annotation_candidate_requires_an_existing_assertion(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            with self.assertRaises(sqlite3.IntegrityError):
                with storage.transaction(database) as connection:
                    connection.execute(
                        "INSERT INTO annotation_candidate(candidate_id, assertion_id, rank) "
                        "VALUES (?, ?, ?)",
                        ("urn:msdial:annotation_candidate:orphan:1", "no-such-assertion", 1),
                    )

    def test_criteria_rule_requires_an_existing_criteria_set(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "catalog.sqlite"
            storage.initialize(database)
            with self.assertRaises(sqlite3.IntegrityError):
                with storage.transaction(database) as connection:
                    connection.execute(
                        "INSERT INTO criteria_rule(criteria_rule_id, criteria_set_id, "
                        "evidence_concept_id) VALUES (?, ?, ?)",
                        ("urn:msdial:criteria_rule:orphan:1", "no-such-set", "SL"),
                    )


if __name__ == "__main__":
    unittest.main()
