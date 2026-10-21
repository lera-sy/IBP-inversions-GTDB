"""Create a synthetic Module 4 package for testing Module 5."""

import csv
import hashlib
import json
import random
from pathlib import Path

ROOT = Path(
    "tests/module5/fixtures/prepared/"
    "escherichia_coli_562/test_3_v1"
)

GENOMES = ROOT / "genomes"
GENOMES.mkdir(parents=True, exist_ok=True)

rng = random.Random(42)
reference = "".join(rng.choices("ACGT", k=2400))

def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

inverted = (
    reference[:800]
    + reverse_complement(reference[800:1400])
    + reference[1400:]
)

assemblies = [
    ("GCF_999900001.1", reference, True),
    ("GCF_999900002.1", reference, False),
    ("GCF_999900003.1", inverted, False),
]

columns = [
    "assembly_id", "fasta_path", "sequence_id",
    "chromosome_accession", "chromosome_length_bp",
    "genome_length_bp", "topology", "is_reference",
    "sha256",
]

with (ROOT / "genomes.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()

    for assembly_id, sequence, is_reference in assemblies:
        sequence_id = assembly_id + "_chromosome"
        fasta = GENOMES / (assembly_id + ".fna")

        lines = [
            sequence[i:i+80]
            for i in range(0, len(sequence), 80)
        ]
        fasta.write_text(">" + sequence_id + "\n"
                         + "\n".join(lines) + "\n")

        writer.writerow({
            "assembly_id": assembly_id,
            "fasta_path": "genomes/" + fasta.name,
            "sequence_id": sequence_id,
            "chromosome_accession": "MOCK_" + assembly_id,
            "chromosome_length_bp": len(sequence),
            "genome_length_bp": len(sequence),
            "topology": "circular",
            "is_reference": str(is_reference).lower(),
            "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
        })

manifest = {
    "schema_version": "1.0",
    "species_id": "escherichia_coli_562",
    "taxid": 562,
    "scientific_name": "Escherichia coli",
    "cohort_id": "test_3_v1",
    "genomes_file": "genomes.tsv",
    "reference_assembly_id": "GCF_999900001.1",
    "analysis_scope": "primary_chromosome",
    "selected_count": 3,
    "accepted_count": 3,
}

(ROOT / "manifest.json").write_text(
    json.dumps(manifest, indent=2) + "\n"
)

print("Mock package created:", ROOT)
print("Genomes: 3")
print("Queries: 2")
print("Known inversion: 600 bp")
