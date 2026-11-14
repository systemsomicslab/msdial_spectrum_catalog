"""Invented campaign unit; no real annotations or library records."""
import json
from pathlib import Path


def write_unit(root, *, key="unit-a", private=False):
    unit = Path(root) / key
    out = unit / "output"
    out.mkdir(parents=True)
    (unit / "provenance").mkdir()
    project = dict(repository="synthetic", accession="SYN001", analysis_unit_id=key,
                   title="Invented experiment", description="Synthetic fixture", public_url="https://example.org/study",
                   license="CC0", source_subrecord_id="subrecord", separation="LC", acquisition_mode="DDA", ion_mode="Negative")
    libraries = [dict(name="Invented-acquired.msp", sha256="a" * 64, role="msp", library_kind="acquired",
                      distribution="private" if private else "public", private=private),
                 dict(name="Invented-predicted.msp", sha256="b" * 64, role="msp", library_kind="predicted", distribution="public")]
    (unit / "provenance/run-manifest.json").write_text(json.dumps(
        dict(schema="msdial-public-reanalysis-run.v1", project=project, ignored_path="X:/ignored")), encoding="utf-8")
    (out / "run-manifest.json").write_text(json.dumps(dict(msdial_console_version="5.synthetic",
        msdial_interactive_version="1.synthetic", analysis_type="LCMS", target_omics=dict(applied="Lipidomics"),
        annotation_pipeline_profile="invented", console=dict(binary_sha256="c" * 64, source_kind="build"),
        libraries=libraries, ignored_path="X:/ignored")), encoding="utf-8")
    (unit / "campaign-record.json").write_text(json.dumps(dict(campaign_id="synthetic-campaign", unit_key=key,
        catalog_unit_id="parent", parent_unit_key="parent", role="part", state="done", manifest_digest="d" * 64,
        libraries_as_run=[{k: v for k, v in lib.items() if k in {"name", "role", "sha256", "library_kind"}} for lib in libraries])), encoding="utf-8")
    header = "Peak ID\tName\tScan\tRT (min)\tPrecursor m/z\tHeight\tArea\tAdduct\n"
    for sample in ("sample_a", "sample_b"):
        (out / (sample + ".mdpeak")).write_text(header +
            "7\tInvented alpha\t100\t2.5\t300.1\t1000\t4000\t[M-H]-\n"
            "8\tInvented beta\t110\t3.0\t350.1\t900\t3000\t[M-H]-\n", encoding="utf-8")
        (out / (sample + ".mdmsp")).write_text(
            "NAME: Invented alpha\nPRECURSORMZ: 300.1\nCOMMENT: |PEAKID=7|MS1SCAN=100|MS2SCAN=101\nNum Peaks: 1\n100 10\n\n", encoding="utf-8")
    (out / "analysis_files.csv").write_text(
        "file_path,file_name,acquisition_type\nX:/raw/sample_a.d,sample_a,DDA\nX:/raw/sample_b.d,sample_b,DDA\n", encoding="utf-8")
    (out / "AlignResult.mdalign").write_text(
        "Alignment ID\tAverage Rt(min)\tAverage Mz\tMetabolite name\tSpectrum reference file name\n"
        "3\t2.5\t300.1\tInvented alpha\tsample_a\n4\t3.0\t350.1\tUnknown\tsample_b\n", encoding="utf-8")
    (out / "AlignResult.mdmsp").write_text(
        "NAME: Invented alpha\nPRECURSORMZ: 300.1\nCOMMENT: |PEAKID=3\nNum Peaks: 1\n100 10\n\n"
        "NAME: Empty\nPRECURSORMZ: 350.1\nCOMMENT: |PEAKID=4\nNum Peaks: 0\n\n", encoding="utf-8")
    (out / "AlignResult.mdpeakid.tsv").write_text("alignment_master_id\tsample_a\tsample_b\n3\t7\t7\n4\t8\t8\n", encoding="utf-8")
    (out / "AlignResult.mdprovenance.tsv").write_text(
        "alignment_master_id\tfile_id\tfile_name\tis_representative\thas_source_peak\tsource_master_peak_id\tmz\n"
        "3\t0\tsample_a\ttrue\ttrue\t7\t300.1\n3\t1\tsample_b\tfalse\ttrue\t7\t300.1\n"
        "4\t0\tsample_a\tfalse\ttrue\t8\t350.1\n4\t1\tsample_b\ttrue\ttrue\t8\t350.1\n", encoding="utf-8")
    (out / "AlignResult.mdcandidate.tsv").write_text(
        "alignment_master_id\tcandidate_rank\tcandidate_count\tis_representative\tannotator_id\tname\tannotation_tag_vs1\tevidence_source\tmeasured_terms\n"
        "3\t1\t1\ttrue\tmsp_high_quality\tInvented alpha\t400\tReferenceSpectrum\t mz; msms \n", encoding="utf-8")
    (out / "AlignResult.mzTab").write_text("SFH\tSMF_ID\texp_mass_to_charge\nSMF\t3\t300.1\nSMF\t4\t350.1\n", encoding="utf-8")
    (out / "SYN001_sample_metadata_reviewed.tsv").write_text(
        "sample_id\tsource_name\traw_file\traw_file_paired_by\tclass_id\tsample / Species\n"
        "sample_a\tInvented A\tsample_a\tstem\tCase\tInvented species\n"
        "sample_b.d\tInvented B\tarchive-entry\tbasename\tControl\tInvented species\n", encoding="utf-8-sig")
    (out / "msp_annotator_settings.tsv").write_bytes(
        (Path(__file__).parent / "fixtures/ws1/annotator_settings.tsv").read_bytes())
    return unit
