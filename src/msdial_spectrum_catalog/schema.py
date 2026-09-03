from __future__ import annotations

from dataclasses import dataclass

SCHEMA_VERSION = 5

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
    rules_json TEXT NOT NULL,
    vocab_version TEXT
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
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    subject_kind TEXT NOT NULL DEFAULT 'spectrum',
    alignment_feature_id TEXT,
    vocab_version TEXT,
    claim_concept_id TEXT,
    notation_verbatim TEXT,
    claim_unresolved INTEGER NOT NULL DEFAULT 0,
    candidate_count INTEGER NOT NULL DEFAULT 0,
    ambiguity_class_id TEXT,
    tool_run_id TEXT
);

CREATE TABLE IF NOT EXISTS annotation_evidence (
    evidence_id TEXT PRIMARY KEY,
    assertion_id TEXT NOT NULL REFERENCES annotation_assertion(assertion_id),
    evidence_tag TEXT NOT NULL,
    evidence_value_json TEXT,
    source_uri TEXT,
    evidence_concept_id TEXT,
    evidence_subtype TEXT,
    metric TEXT,
    measured_value REAL,
    measured_unit TEXT,
    comparison TEXT,
    threshold_value REAL,
    passed INTEGER,
    criteria_rule_id TEXT,
    source_spectrum_id TEXT,
    source_reference_spectrum_id TEXT,
    source_tool_run_id TEXT,
    out_of_distribution INTEGER NOT NULL DEFAULT 0,
    ood_reason TEXT
);

CREATE TABLE IF NOT EXISTS schema_migration (
    version INTEGER PRIMARY KEY,
    applied_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    description TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS annotation_tool_run (
    tool_run_id TEXT PRIMARY KEY,
    run_id TEXT REFERENCES analysis_run(run_id),
    tool_name TEXT NOT NULL,
    tool_version TEXT,
    tool_provenance_json TEXT,
    parameters_json TEXT,
    input_fingerprint TEXT,
    started_at TEXT,
    finished_at TEXT,
    status TEXT NOT NULL DEFAULT 'completed',
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS reference_library (
    library_id TEXT PRIMARY KEY,
    library_name TEXT NOT NULL,
    library_version TEXT,
    library_kind TEXT NOT NULL,
    source_uri TEXT,
    license TEXT,
    sha256 TEXT,
    byte_size INTEGER,
    record_count INTEGER,
    tool_run_id TEXT REFERENCES annotation_tool_run(tool_run_id),
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS reference_spectrum (
    reference_spectrum_id TEXT PRIMARY KEY,
    library_id TEXT NOT NULL REFERENCES reference_library(library_id),
    library_record_index INTEGER,
    record_name TEXT,
    inchikey TEXT,
    inchikey_skeleton TEXT,
    smiles TEXT,
    formula TEXT,
    ontology TEXT,
    precursor_mz REAL,
    precursor_type TEXT,
    ion_mode TEXT,
    instrument_type TEXT,
    instrument_class TEXT,
    collision_energy_raw TEXT,
    collision_energy_value REAL,
    collision_energy_unit TEXT,
    rt_min REAL,
    ccs REAL,
    peak_count INTEGER NOT NULL,
    payload_sha256 TEXT NOT NULL REFERENCES spectrum_blob(payload_sha256),
    UNIQUE(library_id, library_record_index)
);

CREATE TABLE IF NOT EXISTS spectrum_similarity (
    similarity_id TEXT PRIMARY KEY,
    subject_kind_a TEXT NOT NULL,
    subject_id_a TEXT NOT NULL,
    subject_kind_b TEXT NOT NULL,
    subject_id_b TEXT NOT NULL,
    method TEXT NOT NULL,
    method_version TEXT,
    score REAL NOT NULL,
    score_convention TEXT NOT NULL,
    secondary_method TEXT,
    secondary_score REAL,
    matched_peak_count INTEGER,
    mz_tolerance_da REAL,
    tool_run_id TEXT REFERENCES annotation_tool_run(tool_run_id),
    UNIQUE(subject_kind_a, subject_id_a, subject_kind_b, subject_id_b, method)
);

CREATE TABLE IF NOT EXISTS ambiguity_class (
    ambiguity_class_id TEXT PRIMARY KEY,
    definition_id TEXT NOT NULL,
    library_scope_json TEXT NOT NULL,
    blocking_key TEXT,
    condition_scope_json TEXT,
    member_count INTEGER NOT NULL,
    formula_agreement TEXT,
    skeleton_agreement TEXT,
    linkage_rule TEXT NOT NULL,
    min_pairwise_score REAL,
    score_convention TEXT,
    discriminating_mz_json TEXT,
    discriminating_evidence_needed TEXT,
    tool_run_id TEXT REFERENCES annotation_tool_run(tool_run_id),
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ambiguity_class_member (
    ambiguity_class_member_id TEXT PRIMARY KEY,
    ambiguity_class_id TEXT NOT NULL REFERENCES ambiguity_class(ambiguity_class_id),
    reference_spectrum_id TEXT REFERENCES reference_spectrum(reference_spectrum_id),
    inchikey TEXT,
    inchikey_skeleton TEXT,
    smiles TEXT,
    formula TEXT,
    record_name TEXT
);

CREATE TABLE IF NOT EXISTS annotation_candidate (
    candidate_id TEXT PRIMARY KEY,
    assertion_id TEXT NOT NULL REFERENCES annotation_assertion(assertion_id),
    rank INTEGER NOT NULL,
    rank_is_positional INTEGER NOT NULL DEFAULT 0,
    compound_name TEXT,
    formula TEXT,
    inchikey TEXT,
    inchikey_skeleton TEXT,
    smiles TEXT,
    cxsmiles TEXT,
    cxsmiles_validated INTEGER NOT NULL DEFAULT 0,
    external_db_ref TEXT,
    ontology TEXT,
    score REAL,
    score_type TEXT,
    score_gap_to_next REAL,
    tool_run_id TEXT REFERENCES annotation_tool_run(tool_run_id),
    reference_spectrum_id TEXT REFERENCES reference_spectrum(reference_spectrum_id),
    exclusion_status TEXT,
    exclusion_reason TEXT,
    UNIQUE(assertion_id, rank)
);

CREATE TABLE IF NOT EXISTS annotation_claim_component (
    claim_component_id TEXT PRIMARY KEY,
    assertion_id TEXT NOT NULL REFERENCES annotation_assertion(assertion_id),
    ordinal INTEGER NOT NULL,
    claim_concept_id TEXT NOT NULL,
    UNIQUE(assertion_id, ordinal)
);

CREATE TABLE IF NOT EXISTS criteria_rule (
    criteria_rule_id TEXT PRIMARY KEY,
    criteria_set_id TEXT NOT NULL REFERENCES criteria_set(criteria_set_id),
    evidence_concept_id TEXT NOT NULL,
    evidence_token TEXT,
    operational_criterion TEXT,
    metric TEXT NOT NULL DEFAULT '',
    comparison TEXT,
    threshold_value REAL,
    threshold_unit TEXT,
    scope_json TEXT,
    example TEXT,
    notes TEXT,
    UNIQUE(criteria_set_id, evidence_concept_id, metric)
);

CREATE TABLE IF NOT EXISTS msdial_annotation_result (
    msdial_annotation_result_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    subject_kind TEXT NOT NULL,
    subject_id TEXT NOT NULL,
    feature_id TEXT REFERENCES feature(feature_id),
    alignment_feature_id TEXT REFERENCES alignment_feature(alignment_feature_id),
    rank INTEGER NOT NULL DEFAULT 1,
    annotator_id TEXT NOT NULL DEFAULT '',
    metabolite_name TEXT,
    formula TEXT,
    ontology TEXT,
    inchikey TEXT,
    smiles TEXT,
    adduct TEXT,
    annotation_tag TEXT,
    is_rt_matched INTEGER,
    is_mz_matched INTEGER,
    is_msms_matched INTEGER,
    rt_similarity REAL,
    mz_similarity REAL,
    ccs_similarity REAL,
    simple_dot_product REAL,
    weighted_dot_product REAL,
    reverse_dot_product REAL,
    matched_peaks_count REAL,
    matched_peaks_percentage REAL,
    total_score REAL,
    score_convention TEXT,
    annotation_kind TEXT,
    candidate_name TEXT,
    candidate_is_named INTEGER,
    comment TEXT,
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_row INTEGER,
    UNIQUE(run_id, subject_kind, subject_id, rank, annotator_id)
);

CREATE INDEX IF NOT EXISTS idx_feature_sample ON feature(sample_id, master_peak_id);
CREATE INDEX IF NOT EXISTS idx_spectrum_feature ON spectrum(feature_id);
CREATE INDEX IF NOT EXISTS idx_spectrum_alignment ON spectrum(alignment_feature_id);
CREATE INDEX IF NOT EXISTS idx_alignment_member_feature ON alignment_member(feature_id);
CREATE INDEX IF NOT EXISTS idx_mztab_alignment ON mztab_record(alignment_feature_id);
CREATE INDEX IF NOT EXISTS idx_artifact_sha256 ON artifact(sha256);
CREATE INDEX IF NOT EXISTS idx_annotation_assertion_spectrum ON annotation_assertion(spectrum_id);
CREATE INDEX IF NOT EXISTS idx_annotation_evidence_assertion ON annotation_evidence(assertion_id);
CREATE INDEX IF NOT EXISTS idx_annotation_candidate_assertion ON annotation_candidate(assertion_id);
CREATE INDEX IF NOT EXISTS idx_reference_spectrum_library ON reference_spectrum(library_id);
CREATE INDEX IF NOT EXISTS idx_reference_spectrum_skeleton ON reference_spectrum(inchikey_skeleton);
CREATE INDEX IF NOT EXISTS idx_ambiguity_member_class ON ambiguity_class_member(ambiguity_class_id);
CREATE INDEX IF NOT EXISTS idx_msdial_annotation_subject ON msdial_annotation_result(subject_kind, subject_id);
CREATE INDEX IF NOT EXISTS idx_criteria_rule_set ON criteria_rule(criteria_set_id);
"""


@dataclass(frozen=True)
class AddColumn:
    """Additive ALTER TABLE ADD COLUMN step, applied only when the column is absent."""

    table: str
    column: str
    definition: str


@dataclass(frozen=True)
class RunSql:
    """Idempotent statement step, applied verbatim."""

    statement: str


@dataclass(frozen=True)
class Migration:
    version: int
    description: str
    steps: tuple[AddColumn | RunSql, ...]


MIGRATIONS: tuple[Migration, ...] = (
    Migration(
        version=3,
        description="record artifact source paths and sample raw-file paths",
        steps=(
            AddColumn("artifact", "source_path", "TEXT NOT NULL DEFAULT ''"),
            AddColumn("sample", "raw_file_path", "TEXT"),
        ),
    ),
    Migration(
        version=4,
        description="annotation vocabulary, candidate lists, reference libraries and MS-DIAL score blocks",
        steps=(
            AddColumn("annotation_assertion", "subject_kind", "TEXT NOT NULL DEFAULT 'spectrum'"),
            AddColumn("annotation_assertion", "alignment_feature_id", "TEXT"),
            AddColumn("annotation_assertion", "vocab_version", "TEXT"),
            AddColumn("annotation_assertion", "claim_concept_id", "TEXT"),
            AddColumn("annotation_assertion", "notation_verbatim", "TEXT"),
            AddColumn("annotation_assertion", "claim_unresolved", "INTEGER NOT NULL DEFAULT 0"),
            AddColumn("annotation_assertion", "candidate_count", "INTEGER NOT NULL DEFAULT 0"),
            AddColumn("annotation_assertion", "ambiguity_class_id", "TEXT"),
            AddColumn("annotation_assertion", "tool_run_id", "TEXT"),
            AddColumn("annotation_evidence", "evidence_concept_id", "TEXT"),
            AddColumn("annotation_evidence", "evidence_subtype", "TEXT"),
            AddColumn("annotation_evidence", "metric", "TEXT"),
            AddColumn("annotation_evidence", "measured_value", "REAL"),
            AddColumn("annotation_evidence", "measured_unit", "TEXT"),
            AddColumn("annotation_evidence", "comparison", "TEXT"),
            AddColumn("annotation_evidence", "threshold_value", "REAL"),
            AddColumn("annotation_evidence", "passed", "INTEGER"),
            AddColumn("annotation_evidence", "criteria_rule_id", "TEXT"),
            AddColumn("annotation_evidence", "source_spectrum_id", "TEXT"),
            AddColumn("annotation_evidence", "source_reference_spectrum_id", "TEXT"),
            AddColumn("annotation_evidence", "source_tool_run_id", "TEXT"),
            AddColumn("annotation_evidence", "out_of_distribution", "INTEGER NOT NULL DEFAULT 0"),
            AddColumn("annotation_evidence", "ood_reason", "TEXT"),
            AddColumn("criteria_set", "vocab_version", "TEXT"),
        ),
    ),
    Migration(
        version=5,
        description="separate MS-DIAL precursor-only and low-score suggestions from real MS/MS matches",
        steps=(
            AddColumn("msdial_annotation_result", "annotation_kind", "TEXT"),
            AddColumn("msdial_annotation_result", "candidate_name", "TEXT"),
            AddColumn("msdial_annotation_result", "candidate_is_named", "INTEGER"),
        ),
    ),
)
