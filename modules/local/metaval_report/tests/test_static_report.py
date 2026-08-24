#!/usr/bin/env python3

import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
STATIC_REPORT_PATH = REPO_ROOT / "bin" / "static_report.py"


def load_static_report_module():
    spec = importlib.util.spec_from_file_location("static_report", STATIC_REPORT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class StaticReportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.static_report = load_static_report_module()

    def test_build_html_accepts_empty_blastx_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
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
                    "library_type": "DNA",
                    "is_ntc": "false",
                    "batch": "batch1",
                }
            ]

            html = self.static_report.build_html(
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

            self.assertIn("Clinical Metagenomics Report", html)
            self.assertIn("sample1", html)
            self.assertIn("Sample", html)


if __name__ == "__main__":
    unittest.main()
