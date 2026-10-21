"""Select a deterministic Mash medoid using only accepted chromosomes."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from ibp.module4 import select_reference as reference


class ReferenceTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.distances = self.root / "distances.tsv"

    def write_distances(self, matrix):
        rows = [f"{first}.fna\t{second}.fna\t{value}\t0\t100/100\n"
                for (first, second), value in matrix.items()]
        self.distances.write_text("".join(rows))

    def test_medoid_and_small_cohorts(self):
        matrix = {("a", "a"): 0, ("b", "b"): 0, ("c", "c"): 0,
                  ("a", "b"): 0.1, ("b", "a"): 0.1,
                  ("a", "c"): 0.3, ("c", "a"): 0.3,
                  ("b", "c"): 0.2, ("c", "b"): 0.2}
        self.write_distances(matrix)
        ranking = reference.rank_distances(self.distances, ["c", "b", "a"])
        self.assertEqual(ranking["assembly_id"].tolist(), ["b", "a", "c"])
        self.assertAlmostEqual(ranking.iloc[0]["mean_mash_distance"], 0.15)
        self.assertEqual(ranking["is_reference"].tolist(), ["true", "false", "false"])
        for assemblies in [["b", "a"], ["a"]]:
            self.write_distances({pair: value for pair, value in matrix.items() if set(pair) <= set(assemblies)})
            ranking = reference.rank_distances(self.distances, assemblies)
            self.assertEqual(ranking.iloc[0]["assembly_id"], "a")
            if len(assemblies) == 1:
                self.assertTrue(pd.isna(ranking.iloc[0]["mean_mash_distance"]))

    def test_reference_must_be_circular(self):
        matrix = {("a", "a"): 0, ("b", "b"): 0, ("c", "c"): 0,
                  ("a", "b"): 0.1, ("b", "a"): 0.1,
                  ("a", "c"): 0.3, ("c", "a"): 0.3,
                  ("b", "c"): 0.2, ("c", "b"): 0.2}
        self.write_distances(matrix)
        # b is the medoid but linear; a is the best circular one
        ranking = reference.rank_distances(self.distances, ["a", "b", "c"], eligible=["a", "c"])
        self.assertEqual(ranking.set_index("assembly_id").loc["a", "is_reference"], "true")
        self.assertEqual((ranking["is_reference"] == "true").sum(), 1)
        # no circular chromosome at all: fall back to the medoid
        ranking = reference.rank_distances(self.distances, ["a", "b", "c"], eligible=[])
        self.assertEqual(ranking.set_index("assembly_id").loc["b", "is_reference"], "true")

    def test_invalid_distance_tables(self):
        good = {("a", "a"): 0, ("a", "b"): 0.1, ("b", "a"): 0.1, ("b", "b"): 0}
        self.write_distances(good)
        original = self.distances.read_text()
        cases = [original + original.splitlines(keepends=True)[0],
                 "\n".join(original.splitlines()[:-1]), original.replace("a.fna", "foreign.fna"),
                 original.replace("0.1", "nan"), original.replace("0.1", "inf"),
                 original.replace("0.1", "-0.1"), original.replace("0.1", "1.1"),
                 original.replace("0.1", "0.2", 1), original.replace("\t0\t0\t", "\t0.1\t0\t", 1),
                 "malformed\n"]
        for text in cases:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.distances.write_text(text)
                reference.rank_distances(self.distances, ["a", "b"])

    def test_script_accepts_linear_and_circular_but_rejects_invalid_inputs(self):
        fastas = [self.root / "GCF_a.1.fna", self.root / "GCF_b.1.fna"]
        metadata = pd.DataFrame({"Assembly Accession": [path.stem for path in fastas],
                                 "Organism Taxonomic ID": ["1", "1"]})
        report = metadata.copy()
        report["Normalization Status"] = "PASS"
        report["Topology"] = ["linear", "circular"]
        report["Chromosome Length"] = "7"
        report["NCBI Length"] = "7"
        report["Normalized FASTA"] = [str(path) for path in fastas]
        for path in fastas:
            path.write_text(f">{path.stem}\nAACCGTA\n")
        metadata.to_csv(self.root / "metadata.tsv", sep="\t", index=False)
        report.to_csv(self.root / "report.tsv", sep="\t", index=False)
        command = ["select_reference", "--metadata", str(self.root / "metadata.tsv"),
                   "--normalization", str(self.root / "report.tsv"), "--taxid", "1",
                   "--output-dir", str(self.root / "output"), "--fastas", *map(str, fastas)]

        def mash(command, **kwargs):
            self.assertEqual(command[command.index("-p") + 1], "1")
            if command[1] == "sketch":
                self.assertEqual(command[command.index("-k") + 1], "21")
                self.assertEqual(command[command.index("-s") + 1], "10000")
                Path(command[command.index("-o") + 1] + ".msh").touch()
            else:
                for first in fastas:
                    for second in fastas:
                        kwargs["stdout"].write(f"{first}\t{second}\t0\t0\t100/100\n")

        with patch.object(sys, "argv", command), patch.object(reference.subprocess, "run", side_effect=mash):
            reference.main()
        # Tie on distance: the circular chromosome wins (a linear reference cannot be used by the callers).
        self.assertEqual((self.root / "output/reference.txt").read_text(), "GCF_b.1\n")

        for column, value in [("Normalization Status", "REVIEW"), ("Topology", "unknown"),
                              ("Organism Taxonomic ID", "2"), ("Chromosome Length", "8"),
                              ("NCBI Length", "8"), ("Normalized FASTA", str(self.root / "other.fna"))]:
            changed = report.copy()
            changed.loc[0, column] = value
            changed.to_csv(self.root / "report.tsv", sep="\t", index=False)
            with (
                self.subTest(column=column), patch.object(sys, "argv", command),
                patch.object(reference.subprocess, "run") as run, self.assertRaises(ValueError),
            ):
                reference.main()
            run.assert_not_called()
        report.to_csv(self.root / "report.tsv", sep="\t", index=False)
        fastas[0].write_text(">wrong_header\nAACCGTA\n")
        with patch.object(sys, "argv", command), self.assertRaisesRegex(ValueError, "FASTA ID"):
            reference.main()
        fastas[0].unlink()
        with patch.object(sys, "argv", command), self.assertRaises(FileNotFoundError):
            reference.main()


if __name__ == "__main__":
    unittest.main()
