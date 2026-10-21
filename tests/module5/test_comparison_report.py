"""Tests for Module 5 cohort-level comparison reporting."""

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ibp.module5.results import (
    FIELDS,
    write_comparison_report,
    write_standardized,
)

REFERENCE = "GCF_999900001.1"

SYRI_CALL = {
    "ref_start": "801",
    "ref_end": "1400",
    "qry_start": "801",
    "qry_end": "1400",
    "length": "600",
    "syri_id": "INV3",
}


class TestComparisonReport(unittest.TestCase):

    def make_comparison(self, root, query, outcome, count):
        folder = root / query
        folder.mkdir()

        status = {
            "caller": "syri",
            "reference_assembly": REFERENCE,
            "query_assembly": query,
            "status": outcome,
            "call_count": count,
            "reference_length_bp": 2400,
            "query_length_bp": 2400,
            "reason": "simulated_error" if outcome == "failed" else "",
        }

        (folder / "status.json").write_text(json.dumps(status))

        if outcome == "failed":
            # Simulate the placeholder expected from failure handling.
            (folder / "inversions.tsv").write_text(
                "\t".join(FIELDS) + "\n"
            )
        else:
            native = [SYRI_CALL] if outcome == "success" and count else []
            write_standardized(
                folder / "inversions.tsv", status, native
            )

        return (
            "test_species", "test_cohort", "syri",
            REFERENCE, query, folder,
        )

    def test_statuses_are_distinguished(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cases = [
                ("Q0", "success", 0),
                ("Q1", "success", 1),
                ("QU", "unresolved", 0),
                ("QP", "unsupported", 0),
                ("QF", "failed", 0),
            ]

            comparisons = [
                self.make_comparison(root, *case)
                for case in cases
            ]

            report = root / "comparisons.tsv"
            write_comparison_report(report, comparisons)

            with report.open(newline="") as handle:
                rows = {
                    row["query_assembly"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }

            self.assertEqual(len(rows), 5)
            self.assertEqual(rows["Q0"]["call_count"], "0")
            self.assertEqual(rows["Q1"]["call_count"], "1")

            for query in ("QU", "QP", "QF"):
                self.assertEqual(rows[query]["call_count"], "")

            self.assertEqual(rows["QF"]["status"], "failed")
            self.assertEqual(rows["QF"]["reason"], "simulated_error")

    def test_incorrect_call_count_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            comparison = self.make_comparison(
                root, "Q1", "success", 1
            )

            status_path = root / "Q1" / "status.json"
            status = json.loads(status_path.read_text())
            status["call_count"] = 2
            status_path.write_text(json.dumps(status))

            with self.assertRaises(ValueError):
                write_comparison_report(
                    root / "comparisons.tsv", [comparison]
                )



class TestWorkerFailureBoundary(unittest.TestCase):

    def test_missing_fastas_record_failure(self):
        project = Path(__file__).resolve().parents[2]

        with tempfile.TemporaryDirectory() as temp:
            for caller in ("leen", "syri"):
                with self.subTest(caller=caller):
                    output = Path(temp) / caller

                    command = [
                        sys.executable,
                        "-m", f"ibp.module5.{caller}",
                        "--reference", str(Path(temp) / "missing_ref.fna"),
                        "--query", str(Path(temp) / "missing_query.fna"),
                        "--reference-id", "REF",
                        "--query-id", "QUERY",
                        "--reference-topology", "circular",
                        "--query-topology", "circular",
                        "--outdir", str(output),
                    ]

                    env = os.environ.copy()
                    env["PYTHONPATH"] = str(project / "src")

                    process = subprocess.run(
                        command,
                        cwd=project,
                        env=env,
                        capture_output=True,
                        text=True,
                        timeout=30,
                        check=False,
                    )

                    self.assertEqual(
                        process.returncode, 0, process.stderr
                    )

                    status = json.loads(
                        (output / "status.json").read_text()
                    )

                    self.assertEqual(status["status"], "failed")
                    self.assertIsNone(status["call_count"])
                    self.assertIn("FileNotFoundError", status["error"])

                    self.assertTrue(
                        (output / "calls.tsv").is_file()
                    )
                    self.assertTrue(
                        (output / "run.log").is_file()
                    )

                    with (output / "inversions.tsv").open() as handle:
                        rows = list(csv.DictReader(
                            handle, delimiter="\t"
                        ))

                    self.assertEqual(rows, [])

if __name__ == "__main__":
    unittest.main()
