"""Campaign foundation command registration, kept separate for CLI integration."""
from __future__ import annotations

import json
from dataclasses import asdict


def add_ingest_options(parser, *, default="full"):
    parser.add_argument("--granularity", choices=("full", "representative"), default=default)
    parser.add_argument("--no-members", action="store_true")
    parser.add_argument("--mztab", choices=("all", "none"), default="all")


def add_parsers(subparsers):
    unit = subparsers.add_parser("ingest-unit")
    unit.add_argument("database")
    unit.add_argument("unit_directory")
    unit.add_argument("--handoff")
    unit.add_argument("--campaign-id")
    add_ingest_options(unit)
    dev = subparsers.add_parser("ingest-devset")
    dev.add_argument("database")
    dev.add_argument("devset")
    add_ingest_options(dev, default="representative")
    probe = subparsers.add_parser("size-probe")
    probe.add_argument("--unit-dir", action="append", required=True)
    probe.add_argument("--out", required=True)
    privacy = subparsers.add_parser("privacy-scan")
    privacy.add_argument("database")
    report = subparsers.add_parser("build-report")
    report.add_argument("database")
    report.add_argument("--devset")
    report.add_argument("--out")
    attrs = subparsers.add_parser("export-sample-attributes")
    attrs.add_argument("database")
    attrs.add_argument("--unit")
    attrs.add_argument("--out", required=True)
    classes = subparsers.add_parser("show-classes")
    classes.add_argument("database")
    classes.add_argument("subject_id")
    load = subparsers.add_parser("load-structure-classes")
    load.add_argument("--database", required=True)
    load.add_argument("--icd", required=True)


def dispatch(args):
    from .catalog_build import build_report, ingest_devset, privacy_scan, size_probe
    from .linkage import export_sample_attributes, ingest_unit
    from .classes import list_class_predictions
    code = 0
    if args.command == "ingest-unit":
        report = ingest_unit(args.database, args.unit_directory, granularity=args.granularity,
                             members=not args.no_members, mztab=args.mztab,
                             handoff=args.handoff, campaign_id=args.campaign_id)
        result, code = asdict(report), 0 if report.valid else 2
    elif args.command == "ingest-devset":
        report = ingest_devset(args.database, args.devset, granularity=args.granularity,
                               members=not args.no_members, mztab=args.mztab)
        result = {"units": len(report["units"]), "failures": report["failures"],
                  "db_bytes": report["db_bytes"], "table_counts": report.get("table_counts")}
        code = 0 if report["failures"] == 0 else 2
    elif args.command == "size-probe":
        report = size_probe(args.unit_dir, args.out)
        result = {"units": len(report["units"]), "variants": sum(len(u["variants"]) for u in report["units"])}
        code = 0 if all(v["valid"] for u in report["units"] for v in u["variants"]) else 2
    elif args.command == "privacy-scan":
        result = privacy_scan(args.database)
        code = 0 if result["valid"] else 2
    elif args.command == "build-report":
        result = build_report(args.database, devset=args.devset, out=args.out)
        code = 0 if result["privacy"]["valid"] and result["validation_failures"] == 0 and result.get("source_census_matches", True) else 2
    elif args.command == "export-sample-attributes":
        result = {"samples": export_sample_attributes(args.database, args.out, unit=args.unit)}
    elif args.command == "show-classes":
        result = list_class_predictions(args.database, args.subject_id)
    elif args.command == "load-structure-classes":
        from .structure_classes import load_structure_classes
        result = load_structure_classes(args.database, args.icd)
    else:
        return None
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return code
