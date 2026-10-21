import csv
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from ibp.loci.build import (build_families, build_loci, group_overlapping, prepare_calls,
                            process_species, same_breakpoints)
from ibp.module5.results import FIELDS

L = 1_000_000


def leen(q, iv, cid, t="clean", ln=None):
    a, b = map(int, iv.split(";")[0].split("-"))
    return {"caller": "leen", "query_assembly": q, "call_id": cid, "ref_intervals_0": iv,
            "ref_start_0": "", "ref_end_0": "", "call_type": t,
            "reported_length_bp": ln or (b - a if b > a else L - a + b)}


def syri(q, a, b, cid):
    return {"caller": "syri", "query_assembly": q, "call_id": cid, "ref_intervals_0": f"{a}-{b}",
            "ref_start_0": a, "ref_end_0": b, "call_type": "INV", "reported_length_bp": b - a}


class GeometryTest(unittest.TestCase):
    def test_syri_longer_arc_is_folded(self):
        calls = prepare_calls(pd.DataFrame([syri("q", 100_000, 900_000, "s1")]), L, 6000)
        r = calls.iloc[0]
        self.assertEqual((r.start, r.end, r.size_bp), (900_000, 100_000, 200_000))

    def test_floor(self):
        calls = prepare_calls(pd.DataFrame([syri("q", 100_000, 104_000, "s1")]), L, 6000)
        self.assertFalse(calls.iloc[0].kept)

    def test_group_and_wrap(self):
        g = group_overlapping(["950000-20000", "10000-15000", "400000-500000"], L)
        self.assertEqual([(a, b, m) for a, b, m in g], [(400000, 500000, [2]), (950000, 20000, [0, 1])])

    def test_breakpoints_either_order(self):
        self.assertTrue(same_breakpoints(900_000, 100_000, 101_000, 899_000, L, 5000))


class WideLocusTest(unittest.TestCase):
    def test_chain_around_circle_uses_largest_call(self):
        raw = pd.DataFrame([leen("q", "0-300000", "a"), leen("q", "250000-550000", "b"),
                            leen("q", "500000-800000", "c"), leen("q", "750000-1000000", "d")])
        calls = prepare_calls(raw, L, 6000)
        loci = build_loci(calls[calls.kept], L, 5000)
        self.assertEqual(len(loci), 1)
        self.assertTrue(loci.iloc[0]["wide_locus"])
        self.assertEqual(loci.iloc[0]["size_bp"], 300000)


class LociTiersTest(unittest.TestCase):
    def test_tiers_and_families(self):
        raw = pd.DataFrame([
            leen("q1", "100000-300000", "l1"), leen("q1", "150000-160000", "l2", "in_place"),
            syri("q1", 101000, 299000, "s1"),                       # agrees -> high
            leen("q2", "100500-300500", "l3"),                       # same family, Leen only -> low
            leen("q2", "600000-700000", "l4"), syri("q2", 620000, 690000, "s2"),  # overlap -> medium
            syri("q3", 100000, 900000, "s3"),                        # folded to 900000-100000 -> low
        ])
        calls = prepare_calls(raw, L, 6000)
        loci = build_loci(calls[calls.kept], L, 5000)
        tiers = dict(zip(loci["locus"], loci["tier"]))
        self.assertEqual(tiers["q1:1"], "high")
        self.assertEqual(loci.set_index("locus").loc["q1:1", "n_calls"], 3)
        self.assertIn("medium", tiers.values())
        fams, loci = build_families(loci, L, 5000)
        big = fams[fams["n_genomes"] == 2].iloc[0]
        self.assertEqual(big["tier"], "high")       # best tier of its loci
        self.assertEqual(len(fams), 3)


class ProcessSpeciesTest(unittest.TestCase):
    def test_end_to_end_and_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pkg = tmp / "prepared/sp_1/c1"
            pkg.mkdir(parents=True)
            rows = [("ref", "circular", "true"), ("q1", "circular", "false"), ("q2", "linear", "false")]
            pd.DataFrame([{"assembly_id": a, "chromosome_length_bp": L, "genome_length_bp": L + 5000,
                           "topology": t, "is_reference": r} for a, t, r in rows]).to_csv(pkg / "genomes.tsv", sep="\t", index=False)

            def write(caller, q, status, calls):
                d = tmp / "m5/sp_1/c1" / caller / q
                d.mkdir(parents=True)
                (d / "status.json").write_text(json.dumps({"query_assembly": q, "status": status}))
                with open(d / "inversions.tsv", "w", newline="") as h:
                    w = csv.DictWriter(h, fieldnames=FIELDS, delimiter="\t")
                    w.writeheader()
                    for c in calls:
                        w.writerow({k: c.get(k, "") for k in FIELDS})
            write("leen", "q1", "success", [leen("q1", "100000-300000", "l1")])
            write("syri", "q1", "success", [syri("q1", 101000, 299000, "s1")])
            write("leen", "q2", "unsupported", [])
            write("syri", "q2", "success", [syri("q2", 101000, 299000, "s2")])
            summary, comps, calls, loci, fams = process_species(tmp / "m5", tmp / "prepared", "sp_1", "c1",
                                                                6000, 5000, 0.2)
            self.assertEqual(summary["excluded"], "")
            self.assertEqual(summary["genomes_both_callers"], 1)     # linear q2 excluded
            self.assertEqual(list(fams["tier"]), ["high"])


if __name__ == "__main__":
    unittest.main()
