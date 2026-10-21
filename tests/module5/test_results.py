"""Tests for the shared inversion output schema."""

import csv
import tempfile
import unittest
from pathlib import Path

from ibp.module5.results import (
    standardize_calls,
    write_standardized,
)

STATUS = {
    "caller": "leen",
    "status": "success",
    "reference_assembly": "GCF_999900001.1",
    "query_assembly": "GCF_999900003.1",
    "reference_length_bp": 2400,
    "query_length_bp": 2400,
    "call_count": 1,
}

LEEN = {
    "ref_intervals": "801-1391",
    "length": "600",
    "call_type": "clean",
}

SYRI = {
    "ref_start": "801",
    "ref_end": "1400",
    "qry_start": "801",
    "qry_end": "1400",
    "length": "600",
    "syri_id": "INV3",
}


class TestResults(unittest.TestCase):

    def test_leen_estimated_interval(self):
        row = standardize_calls(STATUS, [LEEN])[0]
        self.assertEqual(row["ref_start_0"], 801)
        self.assertEqual(row["ref_end_0"], 1391)
        self.assertEqual(row["reported_length_bp"], 600)
        self.assertEqual(row["interval_span_bp"], 590)

    def test_syri_coordinate_conversion(self):
        status = {**STATUS, "caller": "syri"}
        row = standardize_calls(status, [SYRI])[0]
        self.assertEqual(row["ref_start_0"], 800)
        self.assertEqual(row["ref_end_0"], 1400)
        self.assertEqual(row["qry_start_0"], 800)
        self.assertEqual(row["reported_length_bp"], 600)
        self.assertEqual(row["interval_span_bp"], 600)

    def test_circular_origin_crossing(self):
        call = {
            **LEEN,
            "ref_intervals": "2300-100",
            "length": "200",
        }
        row = standardize_calls(STATUS, [call])[0]
        self.assertEqual(row["interval_span_bp"], 200)
        self.assertEqual(row["spans_reference_origin"], "true")

    def test_multiple_intervals(self):
        call = {
            **LEEN,
            "ref_intervals": "100-300;500-650",
            "length": "350",
        }
        row = standardize_calls(STATUS, [call])[0]
        self.assertEqual(row["ref_start_0"], "")
        self.assertEqual(row["interval_span_bp"], 350)

    def test_unresolved_is_not_zero_inversions(self):
        status = {
            **STATUS, "status": "unresolved", "call_count": 0
        }
        with self.assertRaises(ValueError):
            standardize_calls(status, [])

    def test_successful_empty_output(self):
        status = {**STATUS, "call_count": 0}
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "inversions.tsv"
            write_standardized(path, status, [])

            with path.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(rows, [])


if __name__ == "__main__":
    unittest.main()
