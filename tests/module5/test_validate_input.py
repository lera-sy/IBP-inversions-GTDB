"""Tests for the Module 4 -> Module 5 input contract."""

import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from ibp.module5.validate_input import load_package

FIXTURE = (
    Path(__file__).resolve().parent
    / "fixtures/prepared/escherichia_coli_562/test_3_v1"
)


class TestInputContract(unittest.TestCase):

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)

        self.root = Path(temp.name) / "species"
        shutil.copytree(FIXTURE, self.root)

        self.manifest = self.root / "manifest.json"
        self.inventory = self.root / "genomes.tsv"

    def change_manifest(self, key, value):
        data = json.loads(self.manifest.read_text())
        data[key] = value
        self.manifest.write_text(json.dumps(data))

    def change_inventory(self, row_index, key, value):
        with self.inventory.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames
            rows = list(reader)

        rows[row_index][key] = value

        with self.inventory.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_valid_package(self):
        package = load_package(self.manifest)

        self.assertEqual(package.species_id, "escherichia_coli_562")
        self.assertEqual(package.reference.assembly_id,
                         "GCF_999900001.1")
        self.assertEqual(len(package.queries), 2)

    def test_missing_fasta(self):
        fasta = (
            self.root / "genomes/GCF_999900002.1.fna"
        )
        fasta.unlink()

        with self.assertRaisesRegex(ValueError, "Missing input"):
            load_package(self.manifest)

    def test_wrong_checksum(self):
        self.change_inventory(1, "sha256", "0" * 64)

        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            load_package(self.manifest)

    def test_duplicate_assembly(self):
        self.change_inventory(
            1, "assembly_id", "GCF_999900001.1"
        )

        with self.assertRaisesRegex(ValueError, "duplicate"):
            load_package(self.manifest)

    def test_wrong_reference(self):
        self.change_manifest(
            "reference_assembly_id", "GCF_999900003.1"
        )

        with self.assertRaisesRegex(ValueError, "reference"):
            load_package(self.manifest)

    def test_wrong_genome_count(self):
        self.change_manifest("accepted_count", 2)

        with self.assertRaisesRegex(ValueError, "accepted_count"):
            load_package(self.manifest)

    def test_invalid_topology(self):
        self.change_inventory(1, "topology", "invalid")

        with self.assertRaisesRegex(ValueError, "topology"):
            load_package(self.manifest)

    def test_linear_chromosome_accepted(self):
        self.change_inventory(1, "topology", "linear")

        package = load_package(self.manifest)
        self.assertEqual(package.queries[0].topology, "linear")

    def test_path_traversal_rejected(self):
        self.change_inventory(
            1, "fasta_path", "../../outside.fna"
        )

        with self.assertRaisesRegex(ValueError, "escapes"):
            load_package(self.manifest)


if __name__ == "__main__":
    unittest.main()


class PackageFolderTest(unittest.TestCase):
    def test_folder_or_manifest_both_load(self):
        from pathlib import Path

        from ibp.module5.validate_input import load_package
        folder = Path("tests/module5/fixtures/prepared/escherichia_coli_562/test_3_v1")
        self.assertEqual(load_package(folder).species_id, load_package(folder / "manifest.json").species_id)
