import tempfile
import unittest
from pathlib import Path

import pandas as pd

from ibp.module1.gtdb_select import select, species_key, slug
from ibp.module1.gtdb_species_table import load_and_filter

TAX = "d__Bacteria;p__Pseudomonadota;c__Gamma;o__Ent;f__Ent;g__Escherichia;s__{}"


def gtdb_rows():
    rows = []
    # Species A: 6 genomes, 2 share a BioSample, 1 fails CheckM2
    for i in range(6):
        rows.append({"accession": f"RS_GCF_00000{i:04d}.1", "gtdb_taxonomy": TAX.format("Escherichia coli"),
                     "gtdb_genome_representative": "RS_GCF_000005845.2",
                     "checkm2_completeness": "99", "checkm2_contamination": "0.5",
                     "genome_size": str(5_000_000 + i), "ncbi_assembly_level": "Complete Genome",
                     "ncbi_genome_category": "none", "ncbi_biosample": "SAMA1" if i < 2 else f"SAMA{i}"})
    rows[5]["checkm2_contamination"] = "9"
    # Species B: 3 genomes (too few at min 4), one is a GenBank-only genome, one is a scaffold
    for i in range(3):
        rows.append({"accession": f"RS_GCF_00010{i:04d}.1", "gtdb_taxonomy": TAX.format("Escherichia_A albertii"),
                     "gtdb_genome_representative": "RS_GCF_000100000.1",
                     "checkm2_completeness": "98", "checkm2_contamination": "1",
                     "genome_size": "4600000", "ncbi_assembly_level": "Complete Genome",
                     "ncbi_genome_category": "none", "ncbi_biosample": f"SAMB{i}"})
    rows.append({**rows[-1], "accession": "GB_GCA_000200000.1"})
    rows.append({**rows[-1], "accession": "RS_GCF_000300000.1", "ncbi_assembly_level": "Scaffold"})
    return pd.DataFrame(rows)


def candidates(accessions):
    return pd.DataFrame({"Assembly Accession": accessions, "Organism Taxonomic ID": "562",
                         "Organism Name": "Escherichia coli K-12",
                         "Assembly BioSample Accession": "", "Sequencing Category": "long-read-only"})


class GtdbSelectTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        path = Path(self.tmp.name) / "meta.tsv"
        gtdb_rows().to_csv(path, sep="\t", index=False)
        self.gtdb = load_and_filter(path, 95.0, 5.0)

    def tearDown(self):
        self.tmp.cleanup()

    def test_key_and_slug(self):
        self.assertEqual(species_key("RS_GCF_000005845.2"), 5845)
        self.assertEqual(slug("Escherichia_A albertii"), "escherichia_a_albertii")
        with self.assertRaises(ValueError):
            species_key("not_an_accession")

    def test_filters_dedup_and_threshold(self):
        # GTDB filter keeps 5 of A (one contaminated) and 3 of B (GenBank and scaffold dropped)
        self.assertEqual(len(self.gtdb), 8)
        # Only 4 of A and all of B are long-read
        cand = candidates([f"GCF_00000{i:04d}.1" for i in range(5) if i != 3]
                          + [f"GCF_00010{i:04d}.1" for i in range(3)])
        selected, species, summary, steps, thresholds = select(self.gtdb, cand, 3, 2, seed=1)
        # A: 4 long-read, BioSample duplicate removed -> 3 >= 3, sample 2; B: 3 -> sample 2
        self.assertEqual(sorted(species["species_id"]), ["escherichia_a_albertii_100000", "escherichia_coli_5845"])
        self.assertEqual(len(selected), 4)
        self.assertEqual(set(selected["Organism Taxonomic ID"]), {"5845", "100000"})
        self.assertIn("NCBI Organism Taxonomic ID", selected.columns)
        self.assertEqual(list(steps["genomes"][:3]), [8, 7, 6])

    def test_version_suffix_does_not_block_join(self):
        cand = candidates(["GCF_000000000.2", "GCF_000000002.1", "GCF_000000004.1"])
        selected = select(self.gtdb, cand, 2, 50, seed=1)[0]
        self.assertIn("GCF_000000000.2", set(selected["Assembly Accession"]))
        self.assertTrue(selected.set_index("Assembly Accession").loc["GCF_000000000.2", "version_changed"])

    def test_min_genomes_drops_species(self):
        cand = candidates([f"GCF_00000{i:04d}.1" for i in range(5)] + [f"GCF_00010{i:04d}.1" for i in range(3)])
        selected, species, summary, _, thresholds = select(self.gtdb, cand, 4, 50, seed=1)
        self.assertEqual(list(species["species_id"]), ["escherichia_coli_5845"])
        self.assertEqual(len(selected), 4)  # 5 minus the BioSample duplicate
        self.assertFalse(summary.set_index("gtdb_species").loc["Escherichia_A albertii", "eligible"])
        self.assertEqual(int(thresholds.set_index("min_genomes").loc[5, "species"]), 0)

    def test_same_seed_same_selection(self):
        cand = candidates([f"GCF_00000{i:04d}.1" for i in range(5)])
        a = select(self.gtdb, cand, 2, 2, seed=42)[0]["Assembly Accession"].tolist()
        b = select(self.gtdb, cand, 2, 2, seed=42)[0]["Assembly Accession"].tolist()
        self.assertEqual(a, b)


class GtdbTreeTest(unittest.TestCase):
    def test_prune_and_rename(self):
        try:
            import dendropy  # noqa: F401
        except ImportError:
            self.skipTest("dendropy not installed")
        from ibp.module1.gtdb_tree import prune
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "t.tree").write_text("((RS_GCF_000000001.1:1,RS_GCF_000000002.1:1):1,"
                                        "(RS_GCF_000000003.1:1,RS_GCF_000000004.1:1):1);")
            pd.DataFrame({"gtdb_representative": ["RS_GCF_000000001.1", "RS_GCF_000000003.1",
                                                  "RS_GCF_000000004.1", "RS_GCF_000000009.1"],
                          "species_id": ["a_1", "c_3", "d_4", "x_9"]}).to_csv(tmp / "s.tsv", sep="\t", index=False)
            prune(tmp / "t.tree", tmp / "s.tsv", tmp / "out.tree", tmp / "missing.tsv")
            text = (tmp / "out.tree").read_text()
            for label in ["a_1", "c_3", "d_4"]:
                self.assertIn(label, text)
            self.assertNotIn("RS_GCF", text)
            self.assertEqual(pd.read_csv(tmp / "missing.tsv", sep="\t")["species_id"].tolist(), ["x_9"])


if __name__ == "__main__":
    unittest.main()
