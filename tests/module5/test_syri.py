"""Tests for the SyRI single-comparison worker."""

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ibp.module5.syri import parse_syri

ROOT = Path(__file__).resolve().parents[2]
GENOMES = (
    ROOT / "tests/module5/fixtures/prepared"
    / "escherichia_coli_562/test_3_v1/genomes"
)
REFERENCE = "GCF_999900001.1"


class TestSyriParsing(unittest.TestCase):

    def test_inversion_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "syri.out"
            path.write_text(
                "ref\t801\t1400\t-\t-\tqry\t801\t1400"
                "\tINV3\t-\tINV\t-\n"
            )

            calls, records = parse_syri(path, 2400, True)

            self.assertEqual(records, 1)
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0]["length"], 600)
            self.assertEqual(calls[0]["ref_start"], 801)
            self.assertEqual(calls[0]["ref_end"], 1400)


@unittest.skipUnless(
    os.getenv("IBP_TEST_SYRI") == "1",
    "Real SyRI integration tests require opt-in",
)
class TestSyriWorker(unittest.TestCase):

    def run_worker(self, query_id):
        with tempfile.TemporaryDirectory() as temp:
            env = os.environ.copy()
            env["PYTHONPATH"] = str(ROOT / "src")

            subprocess.run(
                [
                    sys.executable, "-m", "ibp.module5.syri",
                    "--reference",
                    str(GENOMES / f"{REFERENCE}.fna"),
                    "--query",
                    str(GENOMES / f"{query_id}.fna"),
                    "--reference-id", REFERENCE,
                    "--query-id", query_id,
                    "--reference-topology", "circular",
                    "--query-topology", "circular",
                    "--outdir", temp,
                    "--syri-env", "ibp-syri",
                    "--threads", "2",
                ],
                cwd=ROOT,
                env=env,
                check=True,
                capture_output=True,
                text=True,
                timeout=90,
            )

            output = Path(temp)
            status = json.loads(
                (output / "status.json").read_text()
            )

            with (output / "calls.tsv").open() as handle:
                calls = list(csv.DictReader(
                    handle, delimiter="\t"
                ))

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

        self.assertEqual(int(calls[0]["ref_start"]), 801)
        self.assertEqual(int(calls[0]["ref_end"]), 1400)
        self.assertEqual(int(calls[0]["length"]), 600)


if __name__ == "__main__":
    unittest.main()
