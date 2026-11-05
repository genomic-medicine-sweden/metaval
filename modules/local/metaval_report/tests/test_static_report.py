#!/usr/bin/env python3

import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
STATIC_REPORT_PATH = REPO_ROOT / "bin" / "static_report.py"


def load_static_report_module():
    spec = importlib.util.spec_from_file_location("static_report", STATIC_REPORT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_build_html_accepts_empty_blastx_directory(tmp_path):
    static_report = load_static_report_module()

    flagged_dir = tmp_path / "flagged"
    reads_dir = tmp_path / "reads"
    blastn_dir = tmp_path / "blastn"
    blastx_dir = tmp_path / "blastx"
    coverage_dir = tmp_path / "coverage"
    coverage_plots_dir = tmp_path / "coverage_plots"

    for directory in (
        flagged_dir,
        reads_dir,
        blastn_dir,
        blastx_dir,
        coverage_dir,
        coverage_plots_dir,
    ):
        directory.mkdir()

    rows = [
        {
            "sample": "sample1",
            "instrument_platform": "ILLUMINA",
            "na_content": "DNA",
            "is_ntc": "false",
            "sample_prep": "prep1",
        }
    ]

    html = static_report.build_html(
        ticket="1000",
        version="test",
        rows=rows,
        flagged_dir=flagged_dir,
        reads_dir=reads_dir,
        blastn_dir=blastn_dir,
        blastx_dir=blastx_dir,
        coverage_dir=coverage_dir,
        coverage_plots_dir=coverage_plots_dir,
    )

    assert "Clinical Metagenomics Report" in html
    assert "sample1" in html
    assert "Sample" in html


def test_read_taxpasta_table_summarizes_read_counts(tmp_path):
    static_report = load_static_report_module()
    taxpasta_path = tmp_path / "sample1_DNA_prep1_kraken2.tsv"

    taxpasta_path.write_text(
        "\t".join(["taxonomy_id", "name", "rank", "lineage", "sample1"]) + "\n"
        + "\t".join(["10239", "Viruses", "superkingdom", "Viruses", "1200"]) + "\n"
        + "\t".join(["9606", "Homo sapiens", "species", "Eukaryota;Homo sapiens", "300"]) + "\n"
        + "\t".join(["0", "unclassified", "no rank", "unclassified", "50"]) + "\n",
        encoding="utf-8",
    )

    table = static_report.read_taxpasta_table(taxpasta_path, "sample1")

    assert table["total_classified_reads"] == "1,500"
    assert table["total_unclassified_reads"] == "50"
    assert table["host_reads"] == "300"
