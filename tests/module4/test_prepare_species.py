"""Module 5 packages must be complete, portable and published only after validation."""

import errno
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from ibp.module4 import prepare_species as preparation


class PreparationTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.ids = [f"GCF_{number:09d}.1" for number in range(1, 5)]
        selected = pd.DataFrame({"Assembly Accession": self.ids, "Organism Taxonomic ID": "1"})
        selected.to_csv(self.root / "selected.tsv", sep="\t", index=False)
        selected.iloc[:2].to_csv(self.root / "accepted.tsv", sep="\t", index=False)
        qc = selected.copy()
        qc["QC Status"] = ["PASS", "PASS", "PASS", "REVIEW"]
        qc["QC Notes"] = ["", "", "", "missing_length"]
        qc["Chromosome Accession"] = [f"NC_{number:06d}.1" for number in range(1, 5)]
        qc["Chromosome Length"] = "7"
        qc["Assembly Length"] = "12"
        qc["Chromosome FASTA"] = "/scratch/original.fna"
        qc.to_csv(self.root / "qc.tsv", sep="\t", index=False)
        report = qc.iloc[:3].copy()
        report["Normalization Status"] = ["PASS", "PASS", "REVIEW"]
        report["Normalization Notes"] = ["", "", "unknown_topology"]
        report["Topology"] = ["linear", "circular", "unknown"]
        report["NCBI Length"] = "7"
        report["Transformation"] = ["unchanged_linear", "rotation", ""]
        report["Rotation Offset"] = ["0", "2", ""]
        report["Dnaapler Version"] = ["", "1.4.0", ""]
        self.fastas = [self.root / f"{accession}.fna" for accession in self.ids[:2]]
        report["Normalized FASTA"] = [*map(str, self.fastas), ""]
        report.to_csv(self.root / "normalization.tsv", sep="\t", index=False)
        for path in self.fastas:
            path.write_text(f">{path.stem}\nAACCGTA\n")
        (self.root / "reference.txt").write_text(self.ids[1] + "\n")
        pd.DataFrame([{"taxid": "1", "scientific_name": "Example species", "species_id": "example_species_1"}]).to_csv(
            self.root / "species.tsv", sep="\t", index=False)
        self.package = self.root / "prepared/example_species_1/selected_50_v1"
        self.command = ["prepare_species", "--species", "example_species_1", "--cohort", "selected_50_v1",
                        "--species-table", str(self.root / "species.tsv"), "--reference", str(self.root / "reference.txt"),
                        "--output-dir", str(self.package), "--fastas", *map(str, self.fastas)]
        for name in ["selected", "accepted", "qc", "normalization"]:
            self.command.extend(["--" + name, str(self.root / f"{name}.tsv")])

    def publish(self):
        with patch.object(sys, "argv", self.command):
            preparation.main()

    def test_contract_counts_provenance_and_portability(self):
        self.publish()
        preparation.validate_package(self.package)
        manifest = json.loads((self.package / "manifest.json").read_text())
        self.assertEqual((manifest["selected_count"], manifest["accepted_count"]), (4, 2))
        genomes = pd.read_csv(self.package / "genomes.tsv", sep="\t")
        self.assertEqual(genomes["assembly_id"].tolist(), self.ids[:2])
        self.assertEqual(genomes["genome_length_bp"].tolist(), [12, 12])
        self.assertEqual(genomes["chromosome_length_bp"].tolist(), [7, 7])
        self.assertEqual(genomes["topology"].tolist(), ["linear", "circular"])
        self.assertEqual(genomes["is_reference"].tolist(), [False, True])
        self.assertEqual(self.fastas[0].stat().st_ino, (self.package / genomes.iloc[0]["fasta_path"]).stat().st_ino)
        provenance = pd.read_csv(self.package / "preparation.tsv", sep="\t", keep_default_na=False)
        self.assertEqual(provenance["Preparation Status"].tolist(), ["ACCEPTED"] * 2 + ["EXCLUDED"] * 2)
        self.assertEqual(provenance["Transformation"].tolist()[:2], ["unchanged_linear", "rotation"])
        self.assertEqual(provenance["Rotation Offset"].tolist()[:2], ["0", "2"])
        self.assertEqual(provenance.iloc[2]["Normalization Notes"], "unknown_topology")
        self.assertEqual(provenance.iloc[3]["QC Notes"], "missing_length")
        self.assertNotIn("Chromosome FASTA", provenance.columns)
        self.assertNotIn("Normalized FASTA", provenance.columns)
        transferred = self.root / "transferred"
        shutil.copytree(self.package, transferred)
        for path in self.fastas:
            path.unlink()
        preparation.validate_package(transferred)
        preparation.validate_package(self.package)

    def test_existing_package_is_not_overwritten(self):
        self.publish()
        original = (self.package / "manifest.json").read_bytes()
        with self.assertRaisesRegex(FileExistsError, "not overwritten"):
            self.publish()
        self.assertEqual((self.package / "manifest.json").read_bytes(), original)

    def test_cross_filesystem_falls_back_to_copy(self):
        with patch.object(preparation.os, "link", side_effect=OSError(errno.EXDEV, "different filesystem")):
            self.publish()
        preparation.validate_package(self.package)
        self.assertNotEqual(self.fastas[0].stat().st_ino, (self.package / "genomes" / self.fastas[0].name).stat().st_ino)

    def test_failed_validation_does_not_publish_a_partial_package(self):
        with patch.object(preparation, "validate_package", side_effect=ValueError("bad package")), self.assertRaises(ValueError):
            self.publish()
        self.assertFalse(self.package.exists())
        self.assertEqual(list(self.package.parent.glob(".prepare-*")), [])

    def test_invalid_membership_and_metadata_do_not_publish(self):
        path = self.root / "normalization.tsv"
        original = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        for column, value in [("Normalization Status", "REVIEW"), ("Topology", "unknown"),
                              ("Chromosome Length", "8"), ("NCBI Length", "8"),
                              ("Chromosome Accession", "OTHER.1"), ("Normalized FASTA", "/wrong.fna")]:
            modified = original.copy()
            modified.loc[0, column] = value
            modified.to_csv(path, sep="\t", index=False)
            with self.subTest(column=column), self.assertRaises(ValueError):
                self.publish()
            self.assertFalse(self.package.exists())
        original.to_csv(path, sep="\t", index=False)
        (self.root / "reference.txt").write_text(self.ids[2] + "\n")
        with self.assertRaisesRegex(ValueError, "not an accepted"):
            self.publish()

    def test_validator_rejects_bad_paths_flags_lengths_and_checksums(self):
        self.publish()
        table_path = self.package / "genomes.tsv"
        original = pd.read_csv(table_path, sep="\t", dtype=str, keep_default_na=False)
        for column, value in [("fasta_path", "../outside.fna"), ("is_reference", "true"),
                              ("chromosome_length_bp", "8"), ("genome_length_bp", "6"),
                              ("sha256", "0" * 64), ("sequence_id", "wrong"), ("topology", "unknown")]:
            modified = original.copy()
            modified.loc[0, column] = value
            modified.to_csv(table_path, sep="\t", index=False)
            with self.subTest(column=column), self.assertRaises(ValueError):
                preparation.validate_package(self.package)
        original.to_csv(table_path, sep="\t", index=False)
        (self.package / "genomes" / self.fastas[0].name).write_text(f">{self.ids[0]}\nAACCGTT\n")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            preparation.validate_package(self.package)

    def test_validator_rejects_inconsistent_manifest_and_duplicate_assemblies(self):
        self.publish()
        path = self.package / "manifest.json"
        original = json.loads(path.read_text())
        for field, value in [("reference_assembly_id", self.ids[2]), ("accepted_count", 3),
                             ("selected_count", 1), ("schema_version", "other")]:
            path.write_text(json.dumps({**original, field: value}))
            with self.subTest(field=field), self.assertRaises(ValueError):
                preparation.validate_package(self.package)
        path.write_text(json.dumps(original))
        path = self.package / "genomes.tsv"
        genomes = pd.read_csv(path, sep="\t", dtype=str)
        genomes.loc[1, "assembly_id"] = genomes.loc[0, "assembly_id"]
        genomes.to_csv(path, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "Invalid genome table"):
            preparation.validate_package(self.package)


if __name__ == "__main__":
    unittest.main()
