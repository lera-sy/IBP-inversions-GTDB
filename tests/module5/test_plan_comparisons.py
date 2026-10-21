"""Tests for Module 5 comparison planning."""

import unittest
from pathlib import Path

from ibp.module5.plan_comparisons import plan_jobs
from ibp.module5.validate_input import load_package

MANIFEST = (
    Path(__file__).resolve().parent
    / "fixtures/prepared/escherichia_coli_562"
    / "test_3_v1/manifest.json"
)


class TestPlanning(unittest.TestCase):

    def test_four_jobs(self):
        package = load_package(MANIFEST)
        jobs = plan_jobs(package)

        self.assertEqual(len(jobs), 4)

        combinations = {
            (job.caller, job.query_assembly)
            for job in jobs
        }

        self.assertEqual(combinations, {
            ("leen", "GCF_999900002.1"),
            ("syri", "GCF_999900002.1"),
            ("leen", "GCF_999900003.1"),
            ("syri", "GCF_999900003.1"),
        })

    def test_same_reference_for_all_jobs(self):
        package = load_package(MANIFEST)
        jobs = plan_jobs(package)

        self.assertTrue(all(
            job.reference_assembly == "GCF_999900001.1"
            for job in jobs
        ))

    def test_job_order_is_deterministic(self):
        package = load_package(MANIFEST)
        self.assertEqual(
            plan_jobs(package),
            plan_jobs(package),
        )


if __name__ == "__main__":
    unittest.main()
