"""Regression tests for provisional Leen/SyRI support, missingness and plots."""

import tempfile
import unittest
from pathlib import Path

from ibp.module7.aggregate import OUTPUTS, aggregate, read_table, write_table
from ibp.module8.agreement import (
    AGREEMENT_FIELDS,
    CALL_FIELDS,
    FAMILY_SUPPORT_FIELDS,
    build_agreement,
)
from ibp.module8.analyze import analyze
from tests.module7.test_aggregate import fixture


def synthetic_aggregate(folder):
    """Create one species with unique, ambiguous, caller-only, and unknown cases."""
    folder = Path(folder)
    species, cohort = "species_a", "study"
    genomes = [
        {
            "species_id": species, "cohort_id": cohort, "assembly_id": g,
            "chromosome_length_bp": "10000", "genome_length_bp": "11000",
            "is_reference": str(g == "REF").lower(),
            "leen_callable": "true",
            "syri_callable": str(g != "G3").lower(),
        } for g in ("REF", "G1", "G2", "G3")
    ]
    groups = [("leen", "L1", [("G1", "l1"), ("G2", "l2")]),
              ("leen", "L2", [("G1", "l3"), ("G1", "l4")]),
              ("leen", "L3", [("G3", "l5")]),
              ("syri", "S1", [("G1", "s1")]),
              ("syri", "S2", [("G1", "s2")]),
              ("syri", "S3", [("G2", "s3")])]
    families, mappings = [], []
    for caller, fid, calls in groups:
        families.append({
            "species_id": species, "cohort_id": cohort, "caller": caller,
            "family_id": fid, "reference_assembly": "REF", "ref_start_0": "1000",
            "ref_end_0": "2000", "median_reported_length_bp": "1000",
            "observed_carriers": str(len({g for g, _ in calls})),
            "member_calls": str(len(calls)),
            "duplicate_calls_within_genomes": str(len(calls) - len({g for g, _ in calls})),
            "geometry_status": "simple", "interpretation": "candidate_family_not_independent_event",
        })
        for g, call_id in calls:
            mappings.append({
                "species_id": species, "cohort_id": cohort, "caller": caller,
                "family_id": fid, "query_assembly": g, "call_id": call_id,
                "ref_intervals_0": "1000-2000", "reported_length_bp": "1000",
                "geometry_status": "simple",
            })
    edges = []
    for lc, sc, status in [("l1", "s1", "unique_candidate"),
                           ("l3", "s2", "ambiguous_candidate"),
                           ("l4", "s2", "ambiguous_candidate")]:
        edges.append({
            "species_id": species, "cohort_id": cohort, "query_assembly": "G1",
            "leen_call_id": lc, "syri_call_id": sc,
            "leen_family_id": "L1" if lc == "l1" else "L2",
            "syri_family_id": "S1" if sc == "s1" else "S2",
            "left_breakpoint_distance_bp": "0", "right_breakpoint_distance_bp": "0",
            "association_status": status,
            "interpretation": "provisional_breakpoint_proximity_not_consensus",
        })
    summaries = [{
        "species_id": species, "cohort_id": cohort, "caller": caller,
        "total_genomes": "4", "callable_genomes": "4" if caller == "leen" else "3",
        "unknown_genomes": "0" if caller == "leen" else "1",
        "candidate_families": "3", "families_with_multiple_carriers": "1" if caller == "leen" else "0",
        "calls": "5" if caller == "leen" else "3", "complex_geometry_calls": "0",
        "candidate_cross_caller_pairs": "3", "unique_cross_caller_pairs": "1",
        "ambiguous_cross_caller_pairs": "2",
        "interpretation": "tree_free_candidate_families_not_independent_events",
    } for caller in ("leen", "syri")]
    for name, rows in (
        ("species_summary.tsv", summaries),
        ("genome_sizes.tsv", genomes),
        ("candidate_families.tsv", families),
        ("call_family_map.tsv", mappings),
        ("cross_caller_associations.tsv", edges),
        ("tip_states.tsv", []),
    ):
        write_table(folder / name, OUTPUTS[name], rows)
    return folder


class CallerAgreementTests(unittest.TestCase):
    def test_partition_all_call_categories_without_double_counting_ambiguous_edges(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = synthetic_aggregate(tmp)
            rows, calls, families = build_agreement(folder)
            row = rows[0]
            self.assertEqual(int(row["joint_callable_query_genomes"]), 2)
            self.assertEqual(int(row["unique_candidate_pairs"]), 1)
            self.assertEqual({k: int(row[k]) for k in (
                "leen_only_calls", "syri_only_calls", "ambiguous_leen_calls",
                "ambiguous_syri_calls", "unassessable_leen_calls", "unassessable_syri_calls"
            )}, {
                "leen_only_calls": 1, "syri_only_calls": 1,
                "ambiguous_leen_calls": 2, "ambiguous_syri_calls": 1,
                "unassessable_leen_calls": 1, "unassessable_syri_calls": 0,
            })
            mapping = {(r["caller"], r["call_id"]): r for r in calls}
            self.assertEqual(mapping["leen", "l1"]["category"], "unique_candidate")
            self.assertEqual(mapping["leen", "l2"]["category"], "caller_only")
            self.assertEqual(mapping["leen", "l3"]["category"], "ambiguous_candidate")
            self.assertEqual(mapping["leen", "l4"]["possible_partner_call_ids"], "s2")
            self.assertEqual(mapping["leen", "l5"]["category"], "unassessable")
            self.assertEqual(mapping["syri", "s2"]["possible_partner_call_ids"], "l3;l4")
            groups = {(r["caller"], r["family_id"]): r for r in families}
            self.assertEqual(groups["leen", "L1"]["support_class"], "has_unique_candidate")
            self.assertEqual(groups["leen", "L1"]["caller_only_calls"], 1)
            self.assertEqual(groups["leen", "L2"]["support_class"], "ambiguous_only")
            self.assertEqual(groups["leen", "L3"]["support_class"], "not_evaluable")
            self.assertEqual(groups["syri", "S3"]["support_class"],
                             "caller_only_with_evaluable_calls")
            for caller in ("leen", "syri"):
                self.assertEqual(sum(r["caller"] == caller for r in calls),
                                 int(row[f"{caller}_calls"]))

    def test_outputs_include_all_tables_and_five_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = synthetic_aggregate(Path(tmp) / "input")
            out = Path(tmp) / "out"
            analyze(folder, out)
            self.assertEqual(len(read_table(out / "caller_agreement.tsv", AGREEMENT_FIELDS)), 1)
            self.assertEqual(len(read_table(out / "caller_agreement_calls.tsv", CALL_FIELDS)), 8)
            self.assertEqual(len(read_table(out / "caller_agreement_families.tsv",
                                            FAMILY_SUPPORT_FIELDS)), 6)
            for name in ("caller_agreement", "agreement_family_lengths"):
                self.assertGreater((out / "figures" / f"{name}.png").stat().st_size, 5000)

    def test_zero_family_case_and_isolated_source_cohort(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = fixture(tmp, "species_a", cohort="study__species_a", family=False)
            folder = Path(tmp) / "agg"
            aggregate([source], folder, "study", source_cohorts={
                "species_a": "study__species_a",
            })
            rows, calls, families = build_agreement(folder)
            self.assertEqual(rows[0]["unique_candidate_pairs"], 0)
            self.assertEqual(calls, [])
            self.assertEqual(families, [])
            out = Path(tmp) / "out"
            analyze(folder, out)
            self.assertGreater((out / "figures" / "caller_agreement.png").stat().st_size, 5000)

    def test_inconsistent_pair_status_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = synthetic_aggregate(tmp)
            edges = read_table(folder / "cross_caller_associations.tsv",
                               OUTPUTS["cross_caller_associations.tsv"])
            edges[1]["association_status"] = "unique_candidate"
            write_table(folder / "cross_caller_associations.tsv",
                        OUTPUTS["cross_caller_associations.tsv"], edges)
            with self.assertRaisesRegex(ValueError, "cardinality"):
                build_agreement(folder)

    def test_missing_other_caller_not_relabelled_as_caller_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = synthetic_aggregate(tmp)
            rows, calls, _ = build_agreement(folder)
            self.assertEqual(rows[0]["unassessable_leen_calls"], 1)
            call = next(r for r in calls if r["call_id"] == "l5")
            self.assertEqual(call["other_caller_callable"], "false")
            self.assertEqual(call["category"], "unassessable")

    def test_orphan_candidate_pair_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = synthetic_aggregate(tmp)
            edges = read_table(folder / "cross_caller_associations.tsv",
                               OUTPUTS["cross_caller_associations.tsv"])
            edges[0]["syri_call_id"] = "missing"
            write_table(folder / "cross_caller_associations.tsv",
                        OUTPUTS["cross_caller_associations.tsv"], edges)
            with self.assertRaisesRegex(ValueError, "orphan"):
                build_agreement(folder)


if __name__ == "__main__":
    unittest.main()
