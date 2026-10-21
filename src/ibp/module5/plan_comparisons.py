"""Plan independent inversion-calling comparisons."""

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

from ibp.module5.validate_input import SpeciesPackage, load_package

CALLERS = ("leen", "syri")


@dataclass(frozen=True)
class CallerJob:
    species_id: str
    cohort_id: str
    caller: str
    reference_assembly: str
    query_assembly: str
    reference_fasta: Path
    query_fasta: Path
    reference_topology: str
    query_topology: str


def plan_jobs(package: SpeciesPackage) -> tuple[CallerJob, ...]:
    """Build one job per query genome and caller."""

    reference = package.reference

    return tuple(
        CallerJob(
            species_id=package.species_id,
            cohort_id=package.cohort_id,
            caller=caller,
            reference_assembly=reference.assembly_id,
            query_assembly=query.assembly_id,
            reference_fasta=reference.fasta_path,
            query_fasta=query.fasta_path,
            reference_topology=reference.topology,
            query_topology=query.topology,
        )
        for query in sorted(
            package.queries, key=lambda g: g.assembly_id
        )
        for caller in CALLERS
    )


def main():
    parser = argparse.ArgumentParser(
        description="Show planned inversion caller jobs"
    )
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()

    package = load_package(args.manifest)
    jobs = plan_jobs(package)

    writer = csv.writer(
        sys.stdout, delimiter="\t", lineterminator="\n"
    )

    writer.writerow([
        "species_id", "cohort_id", "caller",
        "reference_assembly", "query_assembly",
        "reference_fasta", "query_fasta",
        "reference_topology", "query_topology",
    ])

    for job in jobs:
        writer.writerow([
            job.species_id,
            job.cohort_id,
            job.caller,
            job.reference_assembly,
            job.query_assembly,
            job.reference_fasta,
            job.query_fasta,
            job.reference_topology,
            job.query_topology,
        ])


if __name__ == "__main__":
    main()
