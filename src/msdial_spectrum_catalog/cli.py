from __future__ import annotations

import argparse
import json

from .ingest import ingest_run
from .storage import connect, initialize
from .validate import load_spectrum, validate_run


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
