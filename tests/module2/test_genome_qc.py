"""Small examples covering chromosome selection and unsafe input cases."""

import json
import tempfile
import unittest
from pathlib import Path

from ibp.module2.genome_qc import qc_assembly, read_fasta


class GenomeQCTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.folder = self.root / "ncbi_dataset/data/GCF_test.1"
        self.folder.mkdir(parents=True)
        self.chromosomes = self.root / "chromosomes"
        self.chromosomes.mkdir()

    def run_qc(self, sequences, reports):
        fasta = "".join(f">{header}\n{sequence}\n" for header, sequence in sequences)
        (self.folder / "test_genomic.fna").write_text(fasta)
        (self.folder / "sequence_report.jsonl").write_text(
            "\n".join(json.dumps(report) for report in reports)
        )
        return qc_assembly("GCF_test.1", self.root, self.root, self.chromosomes)

    def report(self, sequence_id, length, molecule="Chromosome"):
        return {"assemblyAccession": "GCF_test.1", "refseqAccession": sequence_id,
                "length": length, "assignedMoleculeLocationType": molecule}

    def test_largest_annotated_chromosome_with_larger_plasmid(self):
        row = self.run_qc(
            [("plasmid", "A" * 20), ("chr2", "AC"), ("chr1 description", "ACGTN")],
            [self.report("chr1", 5), self.report("plasmid", 20, "Plasmid"), self.report("chr2", 2)],
        )
        self.assertEqual(row["QC Status"], "PASS")
        self.assertEqual(row["Chromosome Accession"], "chr1")
        self.assertEqual(row["Excluded Chromosomes"], "chr2")
        self.assertEqual(row["Excluded Plasmids"], "plasmid")
        self.assertEqual(read_fasta(Path(row["Chromosome FASTA"])), [("chr1 description", "ACGTN")])
        self.assertEqual(row["N Percent"], 20)
        self.assertEqual(row["GC Percent"], 40)

    def test_genbank_alias(self):
        report = self.report("refseq", 4)
        report["genbankAccession"] = "genbank"
        row = self.run_qc([("genbank", "ACGT")], [report])
        self.assertEqual(row["QC Status"], "PASS")

    def test_review_cases_do_not_extract_a_chromosome(self):
        missing_length = self.report("chr1", 4)
        del missing_length["length"]
        cases = [
            ([("chr1", "ACGT")], [missing_length], "missing_or_invalid_report_length"),
            ([("chr1", "ACGT"), ("chr2", "AAAA")],
             [self.report("chr1", 4), self.report("chr2", 4)], "largest_chromosome_length_tie"),
            # The assembly totals match, but individual sequence lengths do not.
            ([("chr1", "ACGT"), ("plasmid", "AA")],
             [self.report("chr1", 3), self.report("plasmid", 3, "Plasmid")], "sequence_length_mismatch"),
            ([("chr1", "ACGT"), ("extra", "AA")],
             [self.report("chr1", 4), self.report("extra", 2, "Unknown")], "unclassified_sequence"),
            ([("other", "ACGT")], [self.report("chr1", 4)], "fasta_sequence_id_missing_from_report"),
            ([("chr1", "ACGT"), ("chr1", "ACGT")],
             [self.report("chr1", 4), self.report("chr2", 4)], "duplicate_fasta_sequence"),
            ([("chr1", "ACGR")], [self.report("chr1", 4)], "non_ACGTN_chromosome_bases"),
        ]
        for sequences, reports, reason in cases:
            with self.subTest(reason=reason):
                row = self.run_qc(sequences, reports)
                self.assertEqual(row["QC Status"], "REVIEW")
                self.assertEqual(row["QC Notes"], reason)
                self.assertEqual(list(self.chromosomes.glob("*.fna")), [])

    def test_sequence_before_header_is_rejected(self):
        path = self.root / "bad.fna"
        path.write_text("ACGT\n>chr1\nACGT\n")
        with self.assertRaisesRegex(ValueError, "sequence_before_fasta_header"):
            read_fasta(path)


if __name__ == "__main__":
    unittest.main()
