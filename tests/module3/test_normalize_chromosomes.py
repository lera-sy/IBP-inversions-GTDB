"""Preserve linear DNA unchanged and require a dnaA anchor for circular DNA."""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pandas as pd

from ibp.module3.normalize_chromosomes import (
    check_transformation,
    normalize_one,
    read_chromosome,
)


class NormalizeTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.args = SimpleNamespace(chromosomes=root / "chromosomes", normalized=root / "normalized", logs=root / "logs")
        for directory in vars(self.args).values():
            directory.mkdir()
        self.source = self.args.chromosomes / "GCF_test.1.fna"
        self.source.write_text(">CHR.1 original chromosome\nAACCGTA\n")
        self.row = {"Assembly Accession": "GCF_test.1", "Chromosome Accession": "CHR.1",
                    "Topology Status": "PASS", "Topology Notes": "", "Topology": "circular",
                    "Chromosome Length": "7", "NCBI Length": "7"}

    def test_exact_transformations(self):
        self.assertEqual(check_transformation("AACCGTA", "AACCGTA"), ("unchanged", 0))
        self.assertEqual(check_transformation("AACCGTA", "CCGTAAA"), ("rotation", 2))
        self.assertEqual(check_transformation("AACCGTA", "GGTTTAC"), ("reverse_complement_and_rotation", 3))
        with self.assertRaisesRegex(ValueError, "unexpected_sequence_change"):
            check_transformation("AACCGTA", "AACCGTT")
        with self.assertRaisesRegex(ValueError, "chromosome_length_changed"):
            check_transformation("AACCGTA", "AACCGT")

    def test_anchor_and_sequence_checks_before_publishing(self):
        cases = [
            ("dnaA", "CCGTAAA", "PASS", ""),
            ("Contig_already_reoriented", "AACCGTA", "PASS", ""),
            ("No_MMseqs2_hits", "AACCGTA", "REVIEW", "dnaA_anchor_not_confirmed"),
            ("repA", "CCGTAAA", "REVIEW", "dnaA_anchor_not_confirmed"),
            ("dnaA", "AACCGTT", "REVIEW", "unexpected_sequence_change"),
        ]
        for anchor, sequence, status, note in cases:
            with self.subTest(anchor=anchor, sequence=sequence):
                target = self.args.normalized / "GCF_test.1.fna"
                target.unlink(missing_ok=True)

                def dnaapler(command, sequence=sequence, anchor=anchor, **kwargs):
                    self.assertEqual(command[command.index("--db") + 1], "dnaa")
                    self.assertEqual(command[command.index("--autocomplete") + 1], "none")
                    self.assertEqual(command[command.index("--threads") + 1], "1")
                    work = Path(command[command.index("--output") + 1])
                    work.mkdir()
                    (work / "GCF_test.1_reoriented.fasta").write_text(f">CHR.1 rotated=True\n{sequence}\n")
                    pd.DataFrame([{"Contig": "CHR.1 original chromosome", "Gene_Reoriented": anchor,
                                   "Top_Hit": "protein_DNAA"}]).to_csv(
                        work / "GCF_test.1_all_reorientation_summary.tsv", sep="\t", index=False)
                    return Mock(returncode=0)

                with patch("ibp.module3.normalize_chromosomes.subprocess.run", side_effect=dnaapler):
                    result = normalize_one(self.row, self.args)
                self.assertEqual(result["Normalization Status"], status)
                self.assertEqual(result["Normalization Notes"], note)
                self.assertEqual(target.exists(), status == "PASS")
                self.assertEqual(self.source.read_text(), ">CHR.1 original chromosome\nAACCGTA\n")
                if status == "PASS":
                    self.assertEqual(read_chromosome(target), ("GCF_test.1", sequence))
                self.assertEqual(list(self.args.logs.glob("tmp*")), [])

    def test_linear_chromosome_is_retained_unchanged_without_dnaapler(self):
        self.row["Topology"] = "linear"
        with patch("ibp.module3.normalize_chromosomes.subprocess.run") as run:
            result = normalize_one(self.row, self.args)
        run.assert_not_called()
        self.assertEqual(result["Normalization Status"], "PASS")
        self.assertEqual(result["Normalization Notes"], "")
        self.assertEqual(result["Transformation"], "unchanged_linear")
        self.assertEqual(result["Rotation Offset"], 0)
        self.assertEqual(result["Anchor Result"], "not_applicable")
        self.assertEqual(result["Dnaapler Log"], "")
        self.assertEqual(result["Dnaapler Version"], "")
        self.assertEqual(read_chromosome(Path(result["Normalized FASTA"])), ("GCF_test.1", "AACCGTA"))
        self.assertEqual(self.source.read_text(), ">CHR.1 original chromosome\nAACCGTA\n")
        self.assertEqual(list(self.args.logs.iterdir()), [])

    def test_invalid_linear_and_unknown_chromosomes_are_not_published(self):
        cases = [
            {"Topology": "linear", "NCBI Length": "8"},
            {"Topology": "linear", "Chromosome Length": "8"},
            {"Topology": "linear", "Chromosome Accession": "OTHER.1"},
            {"Topology": "linear", "Topology Status": "REVIEW", "Topology Notes": "chromosome_length_mismatch"},
            {"Topology": "unknown", "Topology Status": "REVIEW", "Topology Notes": "unknown_topology"},
            {"Topology": "unknown"},
        ]
        for changes in cases:
            with self.subTest(changes=changes), patch("ibp.module3.normalize_chromosomes.subprocess.run") as run:
                result = normalize_one({**self.row, **changes}, self.args)
                run.assert_not_called()
                self.assertEqual(result["Normalization Status"], "REVIEW")
                self.assertTrue(result["Normalization Notes"])
                self.assertEqual(list(self.args.normalized.glob("*.fna")), [])

    def test_tool_failure_stops_the_job(self):
        with (
            patch("ibp.module3.normalize_chromosomes.subprocess.run", return_value=Mock(returncode=7)),
            self.assertRaisesRegex(RuntimeError, "code 7"),
        ):
            normalize_one(self.row, self.args)
        self.assertEqual(list(self.args.normalized.glob("*.fna")), [])


if __name__ == "__main__":
    unittest.main()


class NormalizeCacheTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.root = root
        for d in ["chromosomes", "normalized", "logs", "old_fna", "old_logs", "cache_out"]:
            (root / d).mkdir()
        (root / "chromosomes" / "GCF_test.1.fna").write_text(">CHR.1 original\nAACCGTA\n")
        self.row = {"Assembly Accession": "GCF_test.1", "Chromosome Accession": "CHR.1",
                    "Topology Status": "PASS", "Topology Notes": "", "Topology": "circular",
                    "Chromosome Length": "7", "NCBI Length": "7"}
        pd.DataFrame([{"Contig": "CHR.1 original", "Gene_Reoriented": "dnaA", "Top_Hit": "DNAA"}]).to_csv(
            root / "old_logs" / "GCF_test.1.tsv", sep="\t", index=False)
        self.args = SimpleNamespace(chromosomes=root / "chromosomes", normalized=root / "normalized",
                                    logs=root / "logs", cache=[f"{root / 'old_fna'}:{root / 'old_logs'}"],
                                    cache_out=root / "cache_out")

    def test_valid_cache_is_reused_without_dnaapler(self):
        (self.root / "old_fna" / "GCF_test.1.fna").write_text(">GCF_test.1\nCCGTAAA\n")
        with patch("ibp.module3.normalize_chromosomes.subprocess.run") as run:
            result = normalize_one(self.row, self.args)
        run.assert_not_called()
        self.assertEqual(result["Normalization Status"], "PASS")
        self.assertEqual(result["Normalization Notes"], "reused_from_cache")
        self.assertEqual((result["Transformation"], result["Rotation Offset"]), ("rotation", 2))
        self.assertTrue((self.root / "cache_out" / "GCF_test.1.fna").is_file())

    def test_cache_for_a_different_sequence_is_rejected(self):
        (self.root / "old_fna" / "GCF_test.1.fna").write_text(">GCF_test.1\nAACCGTT\n")
        with patch("ibp.module3.normalize_chromosomes.subprocess.run") as run:
            result = normalize_one(self.row, self.args)
        run.assert_not_called()
        self.assertEqual(result["Normalization Status"], "REVIEW")
        self.assertEqual(result["Normalization Notes"], "unexpected_sequence_change")
