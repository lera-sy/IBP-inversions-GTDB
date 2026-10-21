"""Read and validate a Module 4 species/cohort package."""

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

COLUMNS = {
    "assembly_id", "fasta_path", "sequence_id",
    "chromosome_accession", "chromosome_length_bp",
    "genome_length_bp", "topology", "is_reference", "sha256",
}

ALPHABET = set("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class Genome:
    assembly_id: str
    fasta_path: Path
    sequence_id: str
    chromosome_length_bp: int
    genome_length_bp: int
    topology: str
    chromosome_accession: str


@dataclass(frozen=True)
class SpeciesPackage:
    species_id: str
    cohort_id: str
    reference: Genome
    queries: tuple[Genome, ...]


def _package_path(root: Path, relative: str) -> Path:
    path = Path(relative)

    if path.is_absolute() or not relative:
        raise ValueError(f"Expected relative file path: {relative!r}")

    full = (root / path).resolve()

    if not full.is_relative_to(root):
        raise ValueError(f"Path escapes species package: {relative}")

    if not full.is_file():
        raise ValueError(f"Missing input file: {full}")

    return full


def _check_fasta(
    path: Path,
    expected_id: str,
    expected_length: int,
    checksum: str,
) -> None:
    digest = hashlib.sha256()
    header = None
    length = 0

    with path.open("rb") as handle:
        for raw in handle:
            digest.update(raw)
            line = raw.strip().decode("ascii")

            if not line:
                continue

            if line.startswith(">"):
                if header is not None or len(line) == 1:
                    raise ValueError(
                        f"Expected exactly one FASTA record: {path}"
                    )
                header = line[1:].split()[0]
            else:
                if header is None or not set(line.upper()) <= ALPHABET:
                    raise ValueError(f"Invalid FASTA sequence: {path}")
                length += len(line)

    if header != expected_id:
        raise ValueError(
            f"FASTA ID mismatch: {header!r} != {expected_id!r}"
        )

    if length != expected_length:
        raise ValueError(
            f"FASTA length mismatch: {length} != {expected_length}"
        )

    if digest.hexdigest() != checksum.lower():
        raise ValueError(f"FASTA checksum mismatch: {path}")


def load_package(manifest_path: str | Path) -> SpeciesPackage:
    """Load and validate one species/cohort package."""

    manifest_path = Path(manifest_path).resolve()
    # Snakemake can pass the checkpoint package folder instead of the
    # manifest inside it; accept both.
    if manifest_path.is_dir():
        manifest_path = manifest_path / "manifest.json"
    root = manifest_path.parent

    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )

    required = {
        "schema_version", "species_id", "taxid",
        "scientific_name", "cohort_id", "genomes_file",
        "reference_assembly_id", "analysis_scope",
        "selected_count", "accepted_count",
    }

    missing = required - manifest.keys()
    if missing:
        raise ValueError(
            f"Manifest missing fields: {sorted(missing)}"
        )

    if manifest["schema_version"] != "1.0":
        raise ValueError("Unsupported Module 4 schema version")

    inventory = _package_path(root, manifest["genomes_file"])

    with inventory.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        missing = COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Inventory missing columns: {sorted(missing)}"
            )

        rows = list(reader)

    if len(rows) != manifest["accepted_count"]:
        raise ValueError(
            "accepted_count does not match genomes.tsv"
        )

    if not rows or manifest["selected_count"] < len(rows):
        raise ValueError("Invalid cohort counts")

    genomes = []
    references = []
    seen = set()

    for row in rows:
        accession = row["assembly_id"]

        if not accession or accession in seen:
            raise ValueError(
                f"Missing or duplicate assembly ID: {accession!r}"
            )

        seen.add(accession)

        if row["topology"] not in {"circular", "linear", "unknown"}:
            raise ValueError(f"Invalid topology for {accession}")

        if row["is_reference"] not in {"true", "false"}:
            raise ValueError(f"Invalid is_reference for {accession}")

        chromosome_length = int(row["chromosome_length_bp"])
        genome_length = int(row["genome_length_bp"])

        if chromosome_length <= 0 or genome_length < chromosome_length:
            raise ValueError(f"Invalid lengths for {accession}")

        fasta = _package_path(root, row["fasta_path"])

        _check_fasta(
            fasta,
            row["sequence_id"],
            chromosome_length,
            row["sha256"],
        )

        genome = Genome(
            assembly_id=accession,
            fasta_path=fasta,
            sequence_id=row["sequence_id"],
            chromosome_length_bp=chromosome_length,
            genome_length_bp=genome_length,
            topology=row["topology"],
            chromosome_accession=row["chromosome_accession"],
        )

        genomes.append(genome)

        if row["is_reference"] == "true":
            references.append(genome)

    if (
        len(references) != 1
        or references[0].assembly_id
        != manifest["reference_assembly_id"]
    ):
        raise ValueError(
            "Exactly one reference must match reference_assembly_id"
        )

    return SpeciesPackage(
        species_id=manifest["species_id"],
        cohort_id=manifest["cohort_id"],
        reference=references[0],
        queries=tuple(
            genome for genome in genomes
            if genome != references[0]
        ),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate Module 4 input package"
    )
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()

    package = load_package(args.manifest)

    print(f"VALID: {package.species_id} / {package.cohort_id}")
    print(f"Reference: {package.reference.assembly_id}")
    print(f"Queries: {len(package.queries)}")

    for query in package.queries:
        print(f"  {query.assembly_id} ({query.topology})")
