from __future__ import annotations

import argparse
import json

from .ambiguity import ClassDefinition, ambiguity_class_for, compute_ambiguity_classes
from .annotation import list_assertions, load_assertion, validate_annotations
from .ingest import ingest_run
from .reference_library import LIBRARY_KINDS, ingest_reference_library
from .storage import connect, initialize
from .validate import load_spectrum, validate_run
from .vocabulary import (
    CLAIM_AXIS,
    DEFAULT_VOCABULARY,
    EVIDENCE_AXIS,
    NotationError,
    available_versions,
    emit_notation,
    load_vocabulary,
    parse_notation,
    parse_notation_any,
)


def _report_dict(report):
    return {
        "run_id": report.run_id,
        "valid": report.valid,
        "samples": report.samples,
        "features": report.features,
        "spectra": report.spectra,
        "alignments": report.alignments,
        "alignment_members": report.alignment_members,
        "errors": report.errors,
        "warnings": report.warnings,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="MS-DIAL spectrum provenance catalog")
    subparsers = parser.add_subparsers(dest="command", required=True)
    init_parser = subparsers.add_parser("init", help="Initialize a catalog database")
    init_parser.add_argument("database")

    ingest_parser = subparsers.add_parser("ingest-run", help="Register one MS-DIAL output directory")
    ingest_parser.add_argument("database")
    ingest_parser.add_argument("run_directory")
    ingest_parser.add_argument("--repository", required=True)
    ingest_parser.add_argument("--accession", required=True)
    ingest_parser.add_argument("--analysis-unit", required=True)
    ingest_parser.add_argument("--study-title")
    ingest_parser.add_argument("--separation-type")
    ingest_parser.add_argument("--ion-mode")
    ingest_parser.add_argument("--acquisition-type")
    ingest_parser.add_argument("--msdial-version")
    ingest_parser.add_argument("--interactive-version")
    ingest_parser.add_argument("--analysis-files", help="analysis_files.csv used by MS-DIAL")
    ingest_parser.add_argument("--parameter-file", help="MS-DIAL method/parameter file used by the run")

    show_parser = subparsers.add_parser("show-run", help="Show registered run counts")
    show_parser.add_argument("database")
    show_parser.add_argument("run_id")

    validate_parser = subparsers.add_parser("validate-run", help="Validate stored provenance links")
    validate_parser.add_argument("database")
    validate_parser.add_argument("run_id")

    spectrum_parser = subparsers.add_parser("show-spectrum", help="Show one spectrum and its peak payload")
    spectrum_parser.add_argument("database")
    spectrum_parser.add_argument("spectrum_id")

    annotations_parser = subparsers.add_parser(
        "show-annotations", help="Summarize MS-DIAL's own annotation results for one run"
    )
    annotations_parser.add_argument("database")
    annotations_parser.add_argument("run_id")

    assess_parser = subparsers.add_parser(
        "assess-candidates",
        help="Decide why each candidate set holds more than one entry, and what that permits",
    )
    assess_parser.add_argument("database")
    assess_parser.add_argument("run_id")
    assess_parser.add_argument("--definition-label", default="candidate-ambiguity-v1")
    assess_parser.add_argument("--mass-separation-mda", type=float, default=1.0)
    assess_parser.add_argument("--discrimination-margin", type=float, default=0.01)

    candidates_parser = subparsers.add_parser(
        "show-candidates",
        help="Report how often MS-DIAL kept more than one candidate, and list one alignment feature's set",
    )
    candidates_parser.add_argument("database")
    candidates_parser.add_argument("run_id")
    candidates_parser.add_argument(
        "--alignment-feature",
        default=None,
        help="List this alignment feature's candidates instead of summarizing the run",
    )

    validate_annotations_parser = subparsers.add_parser(
        "validate-annotations", help="Validate Level-3 annotation claims, candidates and evidence"
    )
    validate_annotations_parser.add_argument("database")
    validate_annotations_parser.add_argument("run_id")

    assertion_parser = subparsers.add_parser(
        "show-assertion", help="Show one annotation claim with its candidates and evidence"
    )
    assertion_parser.add_argument("database")
    assertion_parser.add_argument("assertion_id")

    assertions_parser = subparsers.add_parser("list-assertions", help="List annotation claims")
    assertions_parser.add_argument("database")
    assertions_parser.add_argument("--run-id")
    assertions_parser.add_argument("--alignment-feature-id")
    assertions_parser.add_argument("--spectrum-id")

    vocabulary_parser = subparsers.add_parser(
        "vocabulary", help="Show a controlled-vocabulary version, or list the available versions"
    )
    vocabulary_parser.add_argument("version", nargs="?")

    notation_parser = subparsers.add_parser(
        "notation", help="Parse an annotation notation such as L3-SC[FM,DF,MN,CO]"
    )
    notation_parser.add_argument("notation")
    notation_parser.add_argument(
        "--vocabulary",
        help=(
            "Vocabulary version. Omit to report every version under which the notation resolves, "
            "which is how a re-bound token such as CP is disambiguated."
        ),
    )
    notation_parser.add_argument("--mode", default="strict", choices=("strict", "permissive", "quarantine"))

    library_parser = subparsers.add_parser(
        "ingest-reference-library",
        help="Register a reference or in-silico MSP library and build skeleton consensus spectra",
    )
    library_parser.add_argument("database")
    library_parser.add_argument("msp_path")
    library_parser.add_argument("--library-name", required=True)
    library_parser.add_argument("--library-version")
    library_parser.add_argument("--kind", default=LIBRARY_KINDS[0], choices=LIBRARY_KINDS)
    library_parser.add_argument("--source-uri")
    library_parser.add_argument("--license")
    library_parser.add_argument("--limit", type=int, help="Cap in-scope records, for a bounded pre-test")
    library_parser.add_argument(
        "--precursor-mz-range",
        nargs=2,
        type=float,
        metavar=("BEGIN", "END"),
        help="Restrict ingestion to a precursor m/z window, for a bounded pre-test",
    )
    library_parser.add_argument("--no-consensus", action="store_true")

    ambiguity_parser = subparsers.add_parser(
        "compute-ambiguity",
        help="Group reference spectra that cannot be told apart, so an annotation can say 'A or B'",
    )
    ambiguity_parser.add_argument("database")
    ambiguity_parser.add_argument("--definition-id", default=ClassDefinition().definition_id)
    ambiguity_parser.add_argument(
        "--weighted-cosine",
        type=float,
        default=ClassDefinition().weighted_cosine_threshold,
        help="Threshold on the symmetric weighted cosine. A convention, not a measurement.",
    )
    ambiguity_parser.add_argument(
        "--entropy", type=float, default=ClassDefinition().entropy_similarity_threshold
    )
    ambiguity_parser.add_argument(
        "--tolerance", type=float, default=ClassDefinition().mz_tolerance_da
    )
    ambiguity_parser.add_argument("--library-id", action="append", dest="library_ids")
    ambiguity_parser.add_argument(
        "--allow-condition-mismatch",
        action="store_true",
        help="Compare across instrument classes and collision-energy bins. Weakens every class it "
             "produces, because an ambiguity class only holds under a stated condition.",
    )

    show_ambiguity_parser = subparsers.add_parser(
        "show-ambiguity", help="Show the ambiguity class anchored on one reference spectrum"
    )
    show_ambiguity_parser.add_argument("database")
    show_ambiguity_parser.add_argument("reference_spectrum_id")

    args = parser.parse_args(argv)
    if args.command == "init":
        initialize(args.database)
        print(args.database)
        return 0
    if args.command == "ingest-run":
        report = ingest_run(
            args.database, args.run_directory, args.repository, args.accession, args.analysis_unit,
            study_title=args.study_title, separation_type=args.separation_type, ion_mode=args.ion_mode,
            acquisition_type=args.acquisition_type, msdial_version=args.msdial_version,
            interactive_version=args.interactive_version,
            analysis_files_csv=args.analysis_files, parameter_file=args.parameter_file,
        )
        print(json.dumps(_report_dict(report), ensure_ascii=False, indent=2))
        return 0 if report.valid else 2
    if args.command == "validate-run":
        report = validate_run(args.database, args.run_id)
        print(json.dumps({
            "run_id": report.run_id,
            "valid": report.valid,
            "counts": report.counts,
            "errors": report.errors,
            "warnings": report.warnings,
        }, ensure_ascii=False, indent=2))
        return 0 if report.valid else 2
    if args.command == "show-spectrum":
        spectrum = load_spectrum(args.database, args.spectrum_id)
        if spectrum is None:
            return 1
        print(json.dumps(spectrum, ensure_ascii=False, indent=2))
        return 0
    if args.command == "ingest-reference-library":
        report = ingest_reference_library(
            args.database,
            args.msp_path,
            library_name=args.library_name,
            library_version=args.library_version,
            library_kind=args.kind,
            source_uri=args.source_uri,
            license=args.license,
            consensus=not args.no_consensus,
            limit=args.limit,
            precursor_mz_range=tuple(args.precursor_mz_range) if args.precursor_mz_range else None,
        )
        print(json.dumps({
            "library_id": report.library_id,
            "valid": report.valid,
            "records_read": report.records_read,
            "records_skipped": report.records_skipped,
            "consensus_spectra": report.consensus_spectra,
            "blobs_written": report.blobs_written,
            "errors": report.errors,
            "warnings": report.warnings,
        }, ensure_ascii=False, indent=2))
        return 0 if report.valid else 2
    if args.command == "compute-ambiguity":
        definition = ClassDefinition(
            definition_id=args.definition_id,
            weighted_cosine_threshold=args.weighted_cosine,
            entropy_similarity_threshold=args.entropy,
            mz_tolerance_da=args.tolerance,
            require_condition_match=not args.allow_condition_mismatch,
        )
        report = compute_ambiguity_classes(
            args.database, definition=definition, library_ids=args.library_ids
        )
        print(json.dumps({
            "definition_id": report.definition_id,
            "valid": report.valid,
            "definition": definition.as_rules(),
            "blocks": report.blocks,
            "pairs_compared": report.pairs_compared,
            "pairs_insufficient_evidence": report.pairs_insufficient_evidence,
            "pairs_condition_mismatch": report.pairs_condition_mismatch,
            "pairs_isobaric_not_isomeric": report.pairs_isobaric_not_isomeric,
            "edges": report.edges,
            "classes": report.classes,
            "singletons": report.singletons,
            "errors": report.errors,
            "warnings": report.warnings,
        }, ensure_ascii=False, indent=2))
        return 0 if report.valid else 2
    if args.command == "show-ambiguity":
        ambiguity_class = ambiguity_class_for(args.database, args.reference_spectrum_id)
        if ambiguity_class is None:
            return 1
        print(json.dumps(ambiguity_class, ensure_ascii=False, indent=2))
        return 0
    if args.command == "validate-annotations":
        report = validate_annotations(args.database, args.run_id)
        print(json.dumps({
            "run_id": report.run_id,
            "valid": report.valid,
            "assertions": report.assertions,
            "candidates": report.candidates,
            "evidence": report.evidence,
            "errors": report.errors,
            "warnings": report.warnings,
        }, ensure_ascii=False, indent=2))
        return 0 if report.valid else 2
    if args.command == "show-assertion":
        assertion = load_assertion(args.database, args.assertion_id)
        if assertion is None:
            return 1
        print(json.dumps(assertion, ensure_ascii=False, indent=2))
        return 0
    if args.command == "list-assertions":
        assertions = list_assertions(
            args.database,
            run_id=args.run_id,
            alignment_feature_id=args.alignment_feature_id,
            spectrum_id=args.spectrum_id,
        )
        print(json.dumps(assertions, ensure_ascii=False, indent=2))
        return 0
    if args.command == "vocabulary":
        if args.version is None:
            print(json.dumps({
                "default": DEFAULT_VOCABULARY,
                "available": list(available_versions()),
            }, ensure_ascii=False, indent=2))
            return 0
        try:
            vocabulary = load_vocabulary(args.version)
        except (KeyError, ValueError) as error:
            print(json.dumps({"error": str(error)}, ensure_ascii=False, indent=2))
            return 1
        print(json.dumps({
            "vocabulary_id": vocabulary.vocabulary_id,
            "version": vocabulary.version,
            "status": vocabulary.status,
            "provenance": vocabulary.provenance,
            "axes": {
                axis: [
                    {"token": term.token, "concept_id": term.concept_id,
                     "label": term.label, "status": term.status}
                    for term in terms
                ]
                for axis, terms in vocabulary.axes.items()
            },
            "open_issues": list(vocabulary.open_issues),
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "notation":
        if args.vocabulary:
            try:
                reading = parse_notation(args.notation, args.vocabulary, mode=args.mode)
            except NotationError as error:
                print(json.dumps({"error": str(error)}, ensure_ascii=False, indent=2))
                return 2
            readings = (reading,)
        else:
            try:
                readings = parse_notation_any(args.notation, mode=args.mode)
            except NotationError as error:
                print(json.dumps({"error": str(error)}, ensure_ascii=False, indent=2))
                return 2
        if not readings:
            print(json.dumps({
                "notation": args.notation,
                "readings": [],
                "note": "No registered vocabulary version resolves every token in this notation.",
            }, ensure_ascii=False, indent=2))
            return 2
        print(json.dumps({
            "notation": args.notation,
            "readings": [
                {
                    "vocab_version": reading.vocab_version,
                    "level": reading.level,
                    "claim_concept_ids": list(reading.claim_concept_ids),
                    "claim_labels": [
                        term.label
                        for concept in reading.claim_concept_ids
                        if (term := load_vocabulary(reading.vocab_version).by_concept(CLAIM_AXIS, concept))
                    ],
                    "evidence_concept_ids": list(reading.evidence_concept_ids),
                    "evidence_labels": [
                        term.label
                        for concept in reading.evidence_concept_ids
                        if (term := load_vocabulary(reading.vocab_version).by_concept(EVIDENCE_AXIS, concept))
                    ],
                    "unknown_tokens": list(reading.unknown_tokens),
                    "unresolved": reading.unresolved,
                    "canonical": emit_notation(reading),
                }
                for reading in readings
            ],
        }, ensure_ascii=False, indent=2))
        return 0
    if args.command == "assess-candidates":
        from .candidate_ambiguity import AssessmentDefinition, assess_run_candidates

        report = assess_run_candidates(
            args.database,
            args.run_id,
            definition=AssessmentDefinition(
                label=args.definition_label,
                mass_separation_mda=args.mass_separation_mda,
                discrimination_margin=args.discrimination_margin,
            ),
        )
        print(json.dumps({
            "run_id": args.run_id,
            "definition_label": report.definition_label,
            "valid": report.valid,
            "subjects": report.subjects,
            "states": report.states,
            "errors": report.errors,
            "warnings": report.warnings,
        }, ensure_ascii=False, indent=2))
        return 0 if report.valid else 2
    if args.command == "show-candidates":
        with connect(args.database) as connection:
            if args.alignment_feature:
                rows = connection.execute(
                    """SELECT candidate_rank, candidate_count, is_representative, name,
                        candidate_is_named, database_id, library_id, formula, inchikey,
                        is_reference_matched, is_annotation_suggested,
                        is_spectrum_comparison_performed, total_score,
                        simple_dot_product, weighted_dot_product, reverse_dot_product,
                        matched_peaks_count, matched_peaks_percentage
                        FROM msdial_annotation_candidate
                        WHERE run_id = ? AND subject_id = ?
                        ORDER BY candidate_rank""",
                    (args.run_id, args.alignment_feature),
                ).fetchall()
                if not rows:
                    return 1
                print(json.dumps({
                    "run_id": args.run_id,
                    "alignment_feature_id": args.alignment_feature,
                    "candidates": [dict(row) for row in rows],
                }, ensure_ascii=False, indent=2))
                return 0
            distribution = connection.execute(
                """SELECT candidate_count, COUNT(DISTINCT subject_id) AS alignment_features
                    FROM msdial_annotation_candidate WHERE run_id = ?
                    GROUP BY candidate_count ORDER BY candidate_count""",
                (args.run_id,),
            ).fetchall()
            if not distribution:
                return 1
            # An annotation with alternatives the search could not separate is not the same claim as one
            # without them, so the split is reported rather than left for a reader to derive.
            annotated = sum(row["alignment_features"] for row in distribution)
            ambiguous = sum(row["alignment_features"] for row in distribution if row["candidate_count"] > 1)
            print(json.dumps({
                "run_id": args.run_id,
                "annotated_alignment_features": annotated,
                "ambiguous_alignment_features": ambiguous,
                "candidate_rows": connection.execute(
                    "SELECT COUNT(*) FROM msdial_annotation_candidate WHERE run_id = ?",
                    (args.run_id,),
                ).fetchone()[0],
                "by_candidate_count": [dict(row) for row in distribution],
            }, ensure_ascii=False, indent=2))
            return 0
    if args.command == "show-annotations":
        with connect(args.database) as connection:
            rows = connection.execute(
                """SELECT subject_kind, annotation_kind,
                    COUNT(*) AS records,
                    SUM(CASE WHEN candidate_is_named = 0 THEN 1 ELSE 0 END) AS unnamed_reference,
                    SUM(CASE WHEN weighted_dot_product IS NOT NULL THEN 1 ELSE 0 END) AS with_dot_product
                    FROM msdial_annotation_result WHERE run_id = ?
                    GROUP BY subject_kind, annotation_kind
                    ORDER BY subject_kind, records DESC""",
                (args.run_id,),
            ).fetchall()
            if not rows:
                return 1
            # Only a named MS/MS match can support spectral-library evidence. A precursor-only row had no
            # product-ion spectrum at all, and an unnamed in-house record identifies no compound.
            eligible = connection.execute(
                """SELECT COUNT(*) FROM msdial_annotation_result
                    WHERE run_id = ? AND annotation_kind = 'msms_matched' AND candidate_is_named = 1""",
                (args.run_id,),
            ).fetchone()[0]
            print(json.dumps({
                "run_id": args.run_id,
                "by_kind": [dict(row) for row in rows],
                "spectral_library_eligible": eligible,
                "total": sum(row["records"] for row in rows),
            }, ensure_ascii=False, indent=2))
            return 0
    with connect(args.database) as connection:
        row = connection.execute(
            """SELECT r.run_id, r.output_directory,
                (SELECT COUNT(*) FROM feature f WHERE f.run_id = r.run_id) AS features,
                (SELECT COUNT(*) FROM spectrum s WHERE s.run_id = r.run_id) AS spectra,
                (SELECT COUNT(*) FROM alignment_feature a WHERE a.run_id = r.run_id) AS alignments
                FROM analysis_run r WHERE r.run_id = ?""",
            (args.run_id,),
        ).fetchone()
        if row is None:
            return 1
        print(json.dumps(dict(row), ensure_ascii=False, indent=2))
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
