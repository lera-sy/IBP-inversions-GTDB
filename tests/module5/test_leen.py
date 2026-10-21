"""Integration tests for the Leen single-comparison worker."""

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENOMES = (
    ROOT / "tests/module5/fixtures/prepared"
    / "escherichia_coli_562/test_3_v1/genomes"
)

REFERENCE = "GCF_999900001.1"


@unittest.skipUnless(shutil.which("minimap2"), "minimap2 required")
class TestLeenWorker(unittest.TestCase):

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.output = Path(temp.name)

    def run_worker(self, query_id, topology="circular"):
        env = os.environ.copy()
        env["PYTHONPATH"] = str(ROOT / "src")

        subprocess.run(
            [
                sys.executable, "-m", "ibp.module5.leen",
                "--reference", str(GENOMES / f"{REFERENCE}.fna"),
                "--query", str(GENOMES / f"{query_id}.fna"),
                "--reference-id", REFERENCE,
                "--query-id", query_id,
                "--reference-topology", "circular",
                "--query-topology", topology,
                "--outdir", str(self.output),
                "--window", "100",
                "--min-size", "200",
                "--max-jump", "200",
                "--threads", "2",
            ],
            cwd=ROOT,
            env=env,
            check=True,
            capture_output=True,
            text=True,
            timeout=90,
        )

        status = json.loads(
            (self.output / "status.json").read_text()
        )

        with (self.output / "calls.tsv").open() as handle:
            calls = list(csv.DictReader(handle, delimiter="\t"))

        return status, calls

    def test_identical_genome(self):
        status, calls = self.run_worker("GCF_999900002.1")

        self.assertEqual(status["status"], "success")
        self.assertEqual(status["call_count"], 0)
        self.assertEqual(calls, [])

    def test_known_inversion(self):
        status, calls = self.run_worker("GCF_999900003.1")

        self.assertEqual(status["status"], "success")
        self.assertEqual(status["call_count"], 1)
        self.assertEqual(len(calls), 1)
        self.assertEqual(int(calls[0]["length"]), 600)

    def test_unsupported_topology(self):
        status, calls = self.run_worker(
            "GCF_999900002.1", topology="linear"
        )

        self.assertEqual(status["status"], "unsupported")
        self.assertEqual(calls, [])


if __name__ == "__main__":
    unittest.main()
