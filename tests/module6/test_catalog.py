"""Regression tests for caller-separated families and explicit unknown states."""

import json
import tempfile
import unittest
from pathlib import Path

from ibp.module5.results import COMPARISON_FIELDS, FIELDS
from ibp.module6.catalog import build_catalogue, endpoints, read_table, write_table

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/module5/fixtures/prepared/escherichia_coli_562/test_3_v1"


class Module6Test(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        manifest = json.loads((FIXTURE / "manifest.json").read_text())
        self.reference = manifest["reference_assembly_id"]
        self.genomes = ["GCF_999900002.1", "GCF_999900003.1"]
        self.rows = {}

    def add(self, caller, genome, records=None, status="success"):
        records = [] if records is None else records
        folder = self.root / caller / genome
        folder.mkdir(parents=True)
        data = []
        for i, (start, end) in enumerate(records, 1):
            data.append({
                "caller": caller, "reference_assembly": self.reference,
                "query_assembly": genome, "call_id": f"{caller}_{genome}_{i}",
                "ref_start_0": start if isinstance(start, int) else "",
                "ref_end_0": end if isinstance(end, int) else "",
                "ref_intervals_0": f"{start}-{end}" if isinstance(start, int) else start,
                "qry_start_0": "", "qry_end_0": "",
                "reported_length_bp": 600, "interval_span_bp": 600,
                "call_type": "INV", "spans_reference_origin": "false",
                "coordinate_method": "syri_native" if caller == "syri" else "window_alignment_estimate",
            })
        write_table(folder / "inversions.tsv", FIELDS, data)
        self.rows[(caller, genome)] = {
            "species_id": "escherichia_coli_562", "cohort_id": "test_3_v1",
            "caller": caller, "reference_assembly": self.reference,
            "query_assembly": genome, "status": status,
            "call_count": str(len(data)) if status == "success" else "",
            "reason": "not_resolved" if status != "success" else "",
            "inversions_path": str(folder / "inversions.tsv"),
            "status_path": str(folder / "status.json"),
        }

    def run_catalogue(self):
        report = self.root / "comparisons.tsv"
        write_table(report, COMPARISON_FIELDS, list(self.rows.values()))
        output = self.root / "output"
        build_catalogue(FIXTURE / "manifest.json", report, output, tolerance=50)
        return output

    def test_family_states_and_caller_associations(self):
        self.add("leen", self.genomes[0])
        self.add("syri", self.genomes[0])
        self.add("leen", self.genomes[1], [(600, 1200)])
        self.add("syri", self.genomes[1], [(620, 1225)])
        folder = self.run_catalogue()
        _, families = read_table(folder / "candidate_families.tsv")
        _, pairings = read_table(folder / "cross_caller_associations.tsv")
        _, states = read_table(folder / "tip_states.tsv")
        self.assertEqual(len(families), 2)
        self.assertEqual(len(pairings), 1)
        self.assertEqual(pairings[0]["association_status"], "unique_candidate")
        self.assertEqual(pairings[0]["left_breakpoint_distance_bp"], "20")
        self.assertEqual(len(states), 6)  # 2 families times 3 genomes
        self.assertEqual({s["state"] for s in states}, {"0", "1"})
        self.assertTrue(all("not_independent_event" in f["interpretation"] for f in families))

    def test_failed_comparison_is_unknown_not_absence(self):
        self.add("leen", self.genomes[0], [(650, 1250)])
        self.add("syri", self.genomes[0], status="failed")
        self.add("leen", self.genomes[1], [(655, 1255)])
        self.add("syri", self.genomes[1], [(660, 1260)])
        folder = self.run_catalogue()
        _, states = read_table(folder / "tip_states.tsv")
        leen = [s for s in states if s["caller"] == "leen"]
        self.assertEqual({s["family_id"] for s in leen}, {"L0001"})
        self.assertEqual({s["state"] for s in leen if s["assembly_id"] in self.genomes}, {"1"})
        syri = [s for s in states if s["caller"] == "syri"]
        self.assertEqual([s["state"] for s in syri if s["assembly_id"] == self.genomes[0]], ["?"])
        _, summary = read_table(folder / "summary.tsv")
        self.assertEqual(next(x for x in summary if x["caller"] == "syri")["unknown_genomes"], "1")

    def test_complex_intervals_not_force_merged(self):
        self.assertEqual(endpoints("2200-2400;0-400", 2400, True), (2200, 400))
        self.assertIsNone(endpoints("100-200;500-600", 2400, True))
        self.add("leen", self.genomes[0], [("100-200;500-600", None)])
        self.add("leen", self.genomes[1], [("100-200;500-600", None)])
        self.add("syri", self.genomes[0])
        self.add("syri", self.genomes[1])
        folder = self.run_catalogue()
        _, families = read_table(folder / "candidate_families.tsv")
        self.assertEqual(len(families), 2)
        self.assertTrue(all(f["geometry_status"] == "complex_unclustered" for f in families))

    def test_one_to_many_cross_caller_is_ambiguous(self):
        self.add("leen", self.genomes[0])
        self.add("syri", self.genomes[0])
        self.add("leen", self.genomes[1], [(600, 1200), (625, 1225)])
        self.add("syri", self.genomes[1], [(610, 1210)])
        folder = self.run_catalogue()
        _, pairings = read_table(folder / "cross_caller_associations.tsv")
        self.assertEqual(len(pairings), 2)
        self.assertEqual({p["association_status"] for p in pairings}, {"ambiguous_candidate"})

    def test_no_calls_still_writes_complete_empty_catalogue(self):
        for caller in ("leen", "syri"):
            for query in self.genomes:
                self.add(caller, query)
        folder = self.run_catalogue()
        _, families = read_table(folder / "candidate_families.tsv")
        _, summary = read_table(folder / "summary.tsv")
        self.assertEqual(families, [])
        self.assertEqual(len(summary), 2)
        self.assertEqual({r["candidate_families"] for r in summary}, {"0"})

    def test_reject_missing_comparison(self):
        self.add("leen", self.genomes[0])
        self.add("syri", self.genomes[0])
        self.add("leen", self.genomes[1])
        with self.assertRaisesRegex(ValueError, "Missing Module 5 comparisons"):
            self.run_catalogue()


if __name__ == "__main__":
    unittest.main()
