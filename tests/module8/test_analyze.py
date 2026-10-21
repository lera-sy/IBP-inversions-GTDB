"""Module 8 smoke tests: empty-family cohorts and nonempty figures."""
import tempfile
import unittest
from pathlib import Path

from ibp.module7.aggregate import aggregate, read_table
from ibp.module8.agreement import AGREEMENT_FIELDS, CALL_FIELDS, FAMILY_SUPPORT_FIELDS
from ibp.module8.analyze import ANALYSIS_FIELDS, SPECIES_METRICS, analyze
from tests.module7.test_aggregate import fixture


class Module8Tests(unittest.TestCase):
    def test_single_species_with_family(self):
        with tempfile.TemporaryDirectory() as temp:
            report = fixture(temp, "species_a")
            out7 = Path(temp) / "module7"
            out8 = Path(temp) / "module8"
            aggregate([report], out7, "test_v1")
            result = analyze(out7, out8)
            self.assertEqual(int(result[0]["total_candidate_families"]), 1)
            self.assertEqual(len(read_table(out8 / "species_metrics.tsv", SPECIES_METRICS)), 2)
            self.assertEqual(len(read_table(out8 / "analysis_summary.tsv", ANALYSIS_FIELDS)), 2)
            rows = read_table(out8 / "caller_agreement.tsv", AGREEMENT_FIELDS)
            self.assertEqual(rows[0]["unique_candidate_pairs"], "1")
            self.assertEqual(len(read_table(out8 / "caller_agreement_calls.tsv", CALL_FIELDS)), 2)
            self.assertEqual(len(read_table(out8 / "caller_agreement_families.tsv", FAMILY_SUPPORT_FIELDS)), 2)
            for name in ("family_lengths", "genome_vs_family_sizes", "families_per_callable_genome",
                         "caller_agreement", "agreement_family_lengths"):
                self.assertGreater((out8 / "figures" / (name + ".png")).stat().st_size, 5000)

    def test_no_families_writes_valid_figures(self):
        with tempfile.TemporaryDirectory() as temp:
            report = fixture(temp, "species_a", family=False)
            out7 = Path(temp) / "module7"
            out8 = Path(temp) / "module8"
            aggregate([report], out7, "test_v1")
            result = analyze(out7, out8)
            self.assertEqual(result[0]["median_candidate_length_bp"], "")
            self.assertEqual(int(result[1]["species_with_size_observations"]), 0)
