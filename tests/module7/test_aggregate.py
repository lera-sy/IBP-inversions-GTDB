"""Module 7 regression tests: aggregation and fail-closed behavior."""
import tempfile
import unittest
from pathlib import Path

from ibp.module6.catalog import (
    ASSOCIATION_FIELDS,
    FAMILY_FIELDS,
    MAP_FIELDS,
    SIZE_FIELDS,
    STATE_FIELDS,
    SUMMARY_FIELDS,
)
from ibp.module7.aggregate import aggregate, read_table, write_table


def fixture(root, species, cohort="test_v1", family=True):
    folder = Path(root) / species / cohort
    folder.mkdir(parents=True)
    genome_ids = ["REF", "QUERY"]
    size_rows = [{"species_id": species, "cohort_id": cohort, "assembly_id": g,
                      "chromosome_length_bp": 10000, "genome_length_bp": 11000,
                      "is_reference": str(g == "REF").lower(),
                      "leen_callable": "true", "syri_callable": "true"} for g in genome_ids]
    families = []
    states = []
    mapping = []
    for caller, fid in (("leen", "L0001"), ("syri", "S0001")):
        if family:
            families.append({"species_id": species, "cohort_id": cohort, "caller": caller,
                                 "family_id": fid, "reference_assembly": "REF", "ref_start_0": "1000",
                                 "ref_end_0": "2000", "median_reported_length_bp": "1000",
                                 "observed_carriers": "1", "member_calls": "1",
                                 "duplicate_calls_within_genomes": "0", "geometry_status": "simple",
                                 "interpretation": "candidate_family_not_independent_event"})
            states.extend({"caller": caller, "family_id": fid, "assembly_id": g,
                               "state": "0" if g == "REF" else "1"} for g in genome_ids)
            mapping.append({"caller": caller, "family_id": fid, "query_assembly": "QUERY",
                                "call_id": f"{caller}_1", "ref_intervals_0": "1000-2000",
                                "reported_length_bp": "1000", "geometry_status": "simple"})
    associations = [{"query_assembly": "QUERY", "leen_call_id": "leen_1", "syri_call_id": "syri_1",
                         "leen_family_id": "L0001", "syri_family_id": "S0001",
                         "left_breakpoint_distance_bp": "0", "right_breakpoint_distance_bp": "0",
                         "association_status": "unique_candidate",
                         "interpretation": "provisional_breakpoint_proximity_not_consensus"}] if family else []
    summary = [{"species_id": species, "cohort_id": cohort, "caller": caller,
                    "total_genomes": "2", "callable_genomes": "2", "unknown_genomes": "0",
                    "candidate_families": "1" if family else "0",
                    "families_with_multiple_carriers": "0", "calls": "1" if family else "0",
                    "complex_geometry_calls": "0",
                    "candidate_cross_caller_pairs": "1" if family else "0",
                    "unique_cross_caller_pairs": "1" if family else "0",
                    "ambiguous_cross_caller_pairs": "0", "interpretation": "tree_free_candidate_families_not_independent_events"}
               for caller in ("leen", "syri")]
    for name, fields, rows in (
        ("summary.tsv", SUMMARY_FIELDS, summary),
        ("genome_sizes.tsv", SIZE_FIELDS, size_rows),
        ("candidate_families.tsv", FAMILY_FIELDS, families),
        ("call_family_map.tsv", MAP_FIELDS, mapping),
        ("tip_states.tsv", STATE_FIELDS, states),
        ("cross_caller_associations.tsv", ASSOCIATION_FIELDS, associations),
    ):
        write_table(folder / name, fields, rows)
    return folder / "summary.tsv"


class Module7Tests(unittest.TestCase):
    def test_aggregate_multiple_species_keeps_caller_distinction(self):
        with tempfile.TemporaryDirectory() as temp:
            a = fixture(temp, "species_a")
            b = fixture(temp, "species_b")
            output = Path(temp) / "aggregation"
            result = aggregate([b, a], output, "test_v1")
            self.assertEqual(len(result), 2)
            self.assertEqual([int(r["candidate_families_total"]) for r in result], [2, 2])
            self.assertEqual([int(r["provisional_cross_caller_pairs_total"]) for r in result], [2, 2])
            rows = read_table(output / "cross_caller_associations.tsv",
                              ["species_id", "cohort_id"] + ASSOCIATION_FIELDS)
            self.assertEqual({r["species_id"] for r in rows}, {"species_a", "species_b"})

    def test_zero_families_is_valid(self):
        with tempfile.TemporaryDirectory() as temp:
            report = fixture(temp, "species_a", family=False)
            result = aggregate([report], Path(temp) / "out", "test_v1")
            self.assertEqual(int(result[0]["candidate_families_total"]), 0)

    def test_missing_species_file_fails_without_silent_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            report = fixture(temp, "species_a")
            (report.parent / "tip_states.tsv").unlink()
            with self.assertRaises(FileNotFoundError):
                aggregate([report], Path(temp) / "out", "test_v1")

    def test_duplicate_species_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            report = fixture(temp, "species_a")
            with self.assertRaises(ValueError):
                aggregate([report, report], Path(temp) / "out", "test_v1")
