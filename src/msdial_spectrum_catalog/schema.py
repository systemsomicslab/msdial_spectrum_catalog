SCHEMA_VERSION = 3

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS catalog_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS study (
    study_id TEXT PRIMARY KEY,
    repository TEXT NOT NULL,
    accession TEXT NOT NULL,
    title TEXT,
    UNIQUE(repository, accession)
);

CREATE TABLE IF NOT EXISTS analysis_unit (
    analysis_unit_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    external_unit_id TEXT NOT NULL,
    separation_type TEXT,
    ion_mode TEXT,
    acquisition_type TEXT,
    UNIQUE(study_id, external_unit_id)
);

CREATE TABLE IF NOT EXISTS analysis_run (
    run_id TEXT PRIMARY KEY,
    analysis_unit_id TEXT NOT NULL REFERENCES analysis_unit(analysis_unit_id),
    run_fingerprint TEXT NOT NULL,
    output_directory TEXT NOT NULL,
    msdial_version TEXT,
    interactive_version TEXT,
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS artifact (
    artifact_id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    artifact_type TEXT NOT NULL,
    relative_path TEXT NOT NULL,
    source_path TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    byte_size INTEGER NOT NULL,
    UNIQUE(run_id, relative_path)
);

CREATE TABLE IF NOT EXISTS sample (
    sample_id TEXT PRIMARY KEY,
    analysis_unit_id TEXT NOT NULL REFERENCES analysis_unit(analysis_unit_id),
    sample_name TEXT NOT NULL,
    raw_file_name TEXT,
    raw_file_path TEXT,
    repository_sample_id TEXT,
    UNIQUE(analysis_unit_id, sample_name)
);

CREATE TABLE IF NOT EXISTS feature (
    feature_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    sample_id TEXT NOT NULL REFERENCES sample(sample_id),
    master_peak_id INTEGER NOT NULL,
    local_peak_id INTEGER,
    ms1_scan_index INTEGER,
    rt_min REAL,
    precursor_mz REAL,
    height REAL,
    area REAL,
    name TEXT,
    adduct TEXT,
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_row INTEGER,
    UNIQUE(run_id, sample_id, master_peak_id)
);

CREATE TABLE IF NOT EXISTS spectrum_blob (
    payload_sha256 TEXT PRIMARY KEY,
    compression TEXT NOT NULL,
    uncompressed_bytes INTEGER NOT NULL,
    payload BLOB NOT NULL
);

CREATE TABLE IF NOT EXISTS alignment_feature (
    alignment_feature_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    alignment_master_id INTEGER NOT NULL,
    alignment_local_id INTEGER,
    parent_alignment_id INTEGER,
    representative_sample_id TEXT REFERENCES sample(sample_id),
    representative_feature_id TEXT REFERENCES feature(feature_id),
    average_rt_min REAL,
    average_mz REAL,
    name TEXT,
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_row INTEGER,
    UNIQUE(run_id, alignment_master_id)
);

CREATE TABLE IF NOT EXISTS spectrum (
    spectrum_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    sample_id TEXT REFERENCES sample(sample_id),
    feature_id TEXT REFERENCES feature(feature_id),
    alignment_feature_id TEXT REFERENCES alignment_feature(alignment_feature_id),
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

CREATE TABLE IF NOT EXISTS alignment_member (
    alignment_member_id TEXT PRIMARY KEY,
    alignment_feature_id TEXT NOT NULL REFERENCES alignment_feature(alignment_feature_id),
    sample_id TEXT NOT NULL REFERENCES sample(sample_id),
    feature_id TEXT REFERENCES feature(feature_id),
    file_id INTEGER NOT NULL,
    is_representative INTEGER NOT NULL,
    has_source_peak INTEGER NOT NULL,
    source_master_peak_id INTEGER,
    source_local_peak_id INTEGER,
    ms1_scan_index INTEGER,
    ms2_scan_index INTEGER,
    rt_min REAL,
    mz REAL,
    height REAL,
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_row INTEGER,
    UNIQUE(alignment_feature_id, sample_id)
);

CREATE TABLE IF NOT EXISTS mztab_record (
    mztab_record_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    section TEXT NOT NULL,
    record_id TEXT NOT NULL,
    parent_refs_json TEXT,
    alignment_feature_id TEXT REFERENCES alignment_feature(alignment_feature_id),
    record_json TEXT NOT NULL,
    source_artifact_id INTEGER NOT NULL REFERENCES artifact(artifact_id),
    source_row INTEGER NOT NULL,
    UNIQUE(run_id, section, record_id)
);

CREATE TABLE IF NOT EXISTS criteria_set (
    criteria_set_id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    version TEXT,
    description TEXT,
    rules_json TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS annotation_assertion (
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

CREATE TABLE IF NOT EXISTS annotation_evidence (
    evidence_id TEXT PRIMARY KEY,
    assertion_id TEXT NOT NULL REFERENCES annotation_assertion(assertion_id),
    evidence_tag TEXT NOT NULL,
    evidence_value_json TEXT,
    source_uri TEXT
);

CREATE INDEX IF NOT EXISTS idx_feature_sample ON feature(sample_id, master_peak_id);
CREATE INDEX IF NOT EXISTS idx_spectrum_feature ON spectrum(feature_id);
CREATE INDEX IF NOT EXISTS idx_spectrum_alignment ON spectrum(alignment_feature_id);
CREATE INDEX IF NOT EXISTS idx_alignment_member_feature ON alignment_member(feature_id);
CREATE INDEX IF NOT EXISTS idx_mztab_alignment ON mztab_record(alignment_feature_id);
CREATE INDEX IF NOT EXISTS idx_artifact_sha256 ON artifact(sha256);
"""
