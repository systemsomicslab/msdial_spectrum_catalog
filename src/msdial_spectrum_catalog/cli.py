from __future__ import annotations

import argparse
import json

from .annotation import list_assertions, load_assertion, validate_annotations
from .ingest import ingest_run
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
