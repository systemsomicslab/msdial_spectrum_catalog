-- v10 public schema (36f1fce), synthetic rows; no Git history needed by tests.
-- Legacy tool ID: urn:msdial:tool-run:invented:1:7d82ea6dfea4dee9
BEGIN TRANSACTION;
CREATE TABLE alignment_feature (
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
INSERT INTO "alignment_feature" VALUES('af','r',3,NULL,NULL,NULL,'f',NULL,NULL,NULL,NULL,NULL);
CREATE TABLE alignment_member (
    alignment_member_id TEXT PRIMARY KEY,
    alignment_feature_id TEXT NOT NULL REFERENCES alignment_feature(alignment_feature_id),
    sample_id TEXT NOT NULL REFERENCES sample(sample_id),
    feature_id TEXT REFERENCES feature(feature_id),
    file_id INTEGER NOT NULL,
    is_representative INTEGER NOT NULL,
    has_source_peak INTEGER NOT NULL,
    source_master_peak_id INTEGER,
    source_local_peak_id INTEGER,
    peak_origin TEXT,
    ms1_scan_index INTEGER,
    ms2_scan_index INTEGER,
    rt_min REAL,
    mz REAL,
    height REAL,
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_row INTEGER,
    UNIQUE(alignment_feature_id, sample_id)
);
CREATE TABLE ambiguity_class (
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
CREATE TABLE ambiguity_class_member (
    ambiguity_class_member_id TEXT PRIMARY KEY,
    ambiguity_class_id TEXT NOT NULL REFERENCES ambiguity_class(ambiguity_class_id),
    reference_spectrum_id TEXT REFERENCES reference_spectrum(reference_spectrum_id),
    inchikey TEXT,
    inchikey_skeleton TEXT,
    smiles TEXT,
    formula TEXT,
    record_name TEXT
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
INSERT INTO "analysis_run" VALUES('r','u','fp','.',NULL,NULL,'2026-11-13 10:54:42');
CREATE TABLE analysis_unit (
    analysis_unit_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    external_unit_id TEXT NOT NULL,
    separation_type TEXT,
    ion_mode TEXT,
    acquisition_type TEXT,
    UNIQUE(study_id, external_unit_id)
);
INSERT INTO "analysis_unit" VALUES('u','s','u',NULL,NULL,NULL);
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
CREATE TABLE annotation_candidate (
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
CREATE TABLE annotation_claim_component (
    claim_component_id TEXT PRIMARY KEY,
    assertion_id TEXT NOT NULL REFERENCES annotation_assertion(assertion_id),
    ordinal INTEGER NOT NULL,
    claim_concept_id TEXT NOT NULL,
    UNIQUE(assertion_id, ordinal)
);
CREATE TABLE annotation_evidence (
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
CREATE TABLE annotation_tool_run (
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
INSERT INTO "annotation_tool_run" VALUES('urn:msdial:tool-run:invented:1:7d82ea6dfea4dee9','r','invented','1',NULL,'{"key":4}',NULL,NULL,NULL,'completed','2026-11-13 10:54:42');
CREATE TABLE artifact (
    artifact_id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    artifact_type TEXT NOT NULL,
    relative_path TEXT NOT NULL,
    source_path TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    byte_size INTEGER NOT NULL,
    UNIQUE(run_id, relative_path)
);
CREATE TABLE candidate_set_assessment (
    candidate_set_assessment_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    subject_kind TEXT NOT NULL,
    subject_id TEXT NOT NULL,
    definition_label TEXT NOT NULL,
    rules_sha256 TEXT NOT NULL,
    state TEXT NOT NULL,
    not_assessed_reason TEXT,
    ceiling_effect TEXT NOT NULL,
    candidate_count INTEGER NOT NULL,
    contender_count INTEGER NOT NULL,
    contender_selection_rule TEXT NOT NULL,
    discrimination_margin REAL,
    max_reference_mz_spread_mda REAL,
    ambiguity_class_id TEXT REFERENCES ambiguity_class(ambiguity_class_id),
    linked_contender_count INTEGER NOT NULL DEFAULT 0,
    unlinked_contender_count INTEGER NOT NULL DEFAULT 0,
    state_reason TEXT,
    assessed_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(run_id, subject_kind, subject_id, definition_label, rules_sha256)
);
CREATE TABLE catalog_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
INSERT INTO "catalog_meta" VALUES('schema_version','10');
CREATE TABLE criteria_rule (
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
CREATE TABLE criteria_set (
    criteria_set_id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    version TEXT,
    description TEXT,
    rules_json TEXT NOT NULL,
    vocab_version TEXT
);
CREATE TABLE feature (
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
INSERT INTO "feature" VALUES('f','r','sample',7,NULL,NULL,NULL,NULL,NULL,NULL,'Invented alpha',NULL,NULL,NULL);
CREATE TABLE msdial_annotation_candidate (
    msdial_annotation_candidate_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_run(run_id),
    subject_kind TEXT NOT NULL,
    subject_id TEXT NOT NULL,
    alignment_feature_id TEXT REFERENCES alignment_feature(alignment_feature_id),
    candidate_rank INTEGER NOT NULL,
    candidate_count INTEGER NOT NULL,
    is_representative INTEGER NOT NULL,
    annotator_id TEXT,
    database_id TEXT,
    source TEXT,
    priority INTEGER,
    library_id INTEGER,
    name TEXT,
    candidate_is_named INTEGER,
    formula TEXT,
    ontology TEXT,
    inchikey TEXT,
    smiles TEXT,
    reference_mz REAL,
    reference_rt_min REAL,
    reference_adduct TEXT,
    annotation_tag TEXT,
    is_reference_matched INTEGER,
    is_annotation_suggested INTEGER,
    is_precursor_mz_match INTEGER,
    is_spectrum_match INTEGER,
    is_spectrum_comparison_performed INTEGER,
    total_score REAL,
    mz_similarity REAL,
    rt_similarity REAL,
    ri_similarity REAL,
    ccs_similarity REAL,
    isotope_similarity REAL,
    simple_dot_product REAL,
    weighted_dot_product REAL,
    reverse_dot_product REAL,
    matched_peaks_count REAL,
    matched_peaks_percentage REAL,
    score_convention TEXT,
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_row INTEGER,
    UNIQUE(run_id, subject_kind, subject_id, candidate_rank)
);
INSERT INTO "msdial_annotation_candidate" VALUES('c','r','alignment_feature','af','af',1,1,1,NULL,NULL,NULL,NULL,NULL,'Invented alpha',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
INSERT INTO "msdial_annotation_candidate" VALUES('cf','r','feature','f',NULL,1,1,1,NULL,NULL,NULL,NULL,NULL,'Invented beta',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
CREATE TABLE msdial_annotation_result (
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
    evidence_source TEXT,
    measured_terms TEXT,
    comment TEXT,
    source_artifact_id INTEGER REFERENCES artifact(artifact_id),
    source_row INTEGER,
    UNIQUE(run_id, subject_kind, subject_id, rank, annotator_id)
);
CREATE TABLE mztab_record (
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
CREATE TABLE reference_consensus_spectrum (
    reference_consensus_spectrum_id TEXT PRIMARY KEY,
    library_id TEXT NOT NULL REFERENCES reference_library(library_id),
    inchikey_skeleton TEXT,
    ion_mode TEXT,
    precursor_type TEXT,
    instrument_class TEXT,
    collision_energy_bin TEXT,
    member_count INTEGER NOT NULL,
    precursor_mz REAL,
    formula TEXT,
    ontology TEXT,
    smiles TEXT,
    peak_count INTEGER NOT NULL,
    mz_bin_width REAL NOT NULL,
    minimum_member_fraction REAL NOT NULL,
    payload_sha256 TEXT NOT NULL REFERENCES spectrum_blob(payload_sha256),
    UNIQUE(library_id, inchikey_skeleton, ion_mode, precursor_type, instrument_class, collision_energy_bin)
);
CREATE TABLE reference_library (
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
CREATE TABLE reference_spectrum (
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
CREATE TABLE sample (
    sample_id TEXT PRIMARY KEY,
    analysis_unit_id TEXT NOT NULL REFERENCES analysis_unit(analysis_unit_id),
    sample_name TEXT NOT NULL,
    raw_file_name TEXT,
    raw_file_path TEXT,
    repository_sample_id TEXT,
    UNIQUE(analysis_unit_id, sample_name)
);
INSERT INTO "sample" VALUES('sample','u','invented',NULL,NULL,NULL);
CREATE TABLE schema_migration (
    version INTEGER PRIMARY KEY,
    applied_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    description TEXT NOT NULL
);
INSERT INTO "schema_migration" VALUES(3,'2026-11-13 10:54:42','record artifact source paths and sample raw-file paths');
INSERT INTO "schema_migration" VALUES(4,'2026-11-13 10:54:42','annotation vocabulary, candidate lists, reference libraries and MS-DIAL score blocks');
INSERT INTO "schema_migration" VALUES(5,'2026-11-13 10:54:42','separate MS-DIAL precursor-only and low-score suggestions from real MS/MS matches');
INSERT INTO "schema_migration" VALUES(6,'2026-11-13 10:54:42','name why an alignment member has no source peak instead of encoding it in a sentinel');
INSERT INTO "schema_migration" VALUES(7,'2026-11-13 10:54:42','keep every annotation candidate MS-DIAL retained, not only the representative');
INSERT INTO "schema_migration" VALUES(8,'2026-11-13 10:54:42','give the library consensus spectra a row, so they stop being unreachable blobs');
INSERT INTO "schema_migration" VALUES(9,'2026-11-13 10:54:42','record why each candidate set holds more than one entry, and what that does to its claim');
INSERT INTO "schema_migration" VALUES(10,'2026-11-13 10:54:42','record the evidence MS-DIAL states in a column instead of a name prefix');
CREATE TABLE spectrum (
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
CREATE TABLE spectrum_blob (
    payload_sha256 TEXT PRIMARY KEY,
    compression TEXT NOT NULL,
    uncompressed_bytes INTEGER NOT NULL,
    payload BLOB NOT NULL
);
CREATE TABLE spectrum_similarity (
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
CREATE TABLE study (
    study_id TEXT PRIMARY KEY,
    repository TEXT NOT NULL,
    accession TEXT NOT NULL,
    title TEXT,
    UNIQUE(repository, accession)
);
INSERT INTO "study" VALUES('s','synthetic','a',NULL);
CREATE INDEX idx_feature_sample ON feature(sample_id, master_peak_id);
CREATE INDEX idx_spectrum_feature ON spectrum(feature_id);
CREATE INDEX idx_spectrum_alignment ON spectrum(alignment_feature_id);
CREATE INDEX idx_alignment_member_feature ON alignment_member(feature_id);
CREATE INDEX idx_mztab_alignment ON mztab_record(alignment_feature_id);
CREATE INDEX idx_artifact_sha256 ON artifact(sha256);
CREATE INDEX idx_annotation_assertion_spectrum ON annotation_assertion(spectrum_id);
CREATE INDEX idx_annotation_evidence_assertion ON annotation_evidence(assertion_id);
CREATE INDEX idx_annotation_candidate_assertion ON annotation_candidate(assertion_id);
CREATE INDEX idx_reference_spectrum_library ON reference_spectrum(library_id);
CREATE INDEX idx_reference_consensus_library ON reference_consensus_spectrum(library_id);
CREATE INDEX idx_reference_spectrum_skeleton ON reference_spectrum(inchikey_skeleton);
CREATE INDEX idx_ambiguity_member_class ON ambiguity_class_member(ambiguity_class_id);
CREATE INDEX idx_msdial_annotation_subject ON msdial_annotation_result(subject_kind, subject_id);
CREATE INDEX idx_msdial_candidate_subject ON msdial_annotation_candidate(subject_kind, subject_id);
CREATE INDEX idx_candidate_assessment_run ON candidate_set_assessment(run_id, state);
CREATE INDEX idx_criteria_rule_set ON criteria_rule(criteria_set_id);
COMMIT;
