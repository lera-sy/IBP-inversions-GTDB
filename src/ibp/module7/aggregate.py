"""Aggregate Module 6 inversion-only catalogues across species.

A caller-specific candidate family is not an independently evolved inversion.
Missing comparisons within a species remain unknown in Module 6 summaries.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

from ibp.module6.catalog import (
    ASSOCIATION_FIELDS,
    FAMILY_FIELDS,
    MAP_FIELDS,
    SIZE_FIELDS,
    STATE_FIELDS,
    SUMMARY_FIELDS,
)

DATASETS = {
    "candidate_families.tsv": FAMILY_FIELDS,
    "call_family_map.tsv": MAP_FIELDS,
    "tip_states.tsv": STATE_FIELDS,
    "cross_caller_associations.tsv": ASSOCIATION_FIELDS,
    "genome_sizes.tsv": SIZE_FIELDS,
    "summary.tsv": SUMMARY_FIELDS,
}
EXTRA_ID = ["species_id", "cohort_id"]
OUTPUTS = {
    "candidate_families.tsv": FAMILY_FIELDS,
    "call_family_map.tsv": EXTRA_ID + MAP_FIELDS,
    "tip_states.tsv": EXTRA_ID + STATE_FIELDS,
    "cross_caller_associations.tsv": EXTRA_ID + ASSOCIATION_FIELDS,
    "genome_sizes.tsv": SIZE_FIELDS,
    "species_summary.tsv": SUMMARY_FIELDS,
}
RUN_FIELDS = [
    "cohort_id", "caller", "species_completed", "genomes_total_across_species",
    "callable_genomes_total", "unknown_genomes_total", "candidate_families_total",
    "families_with_multiple_carriers_total", "calls_total",
    "complex_geometry_calls_total", "provisional_cross_caller_pairs_total",
    "unique_cross_caller_pairs_total", "ambiguous_cross_caller_pairs_total",
    "interpretation",
]


def read_table(path, fields):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != fields:
            raise ValueError(f"Schema mismatch in {path}: {reader.fieldnames}")
        return list(reader)


def write_table(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def aggregate(summaries, outdir, cohort, source_cohorts=None):
    """Validate and combine catalogues, preserving caller and missingness semantics.

    ``source_cohorts`` is for isolated production runs: a mapping from species
    IDs to their original per-species cohort IDs. The published aggregate uses
    ``cohort`` but the original source-cohort identities remain available in
    the production run's species_status.tsv. The default remains strict and
    backwards-compatible with Module 7's Snakemake rule.
    """
    if not summaries:
        raise ValueError("Module 7 requires at least one species catalogue")
    by_name = {name: [] for name in OUTPUTS}
    seen_species = set()
    counts = Counter()
    for summary_file in sorted(map(Path, summaries), key=lambda p: str(p)):
        if summary_file.name != "summary.tsv":
            raise ValueError(f"Expected Module 6 summary.tsv, found {summary_file}")
        folder = summary_file.parent
        tables = {name: read_table(folder / name, fields)
                  for name, fields in DATASETS.items()}
        summary = tables["summary.tsv"]
        if {r["caller"] for r in summary} != {"leen", "syri"} or len(summary) != 2:
            raise ValueError(f"Expected one Leen and one SyRI row in {summary_file}")
        species_values = {r["species_id"] for r in summary}
        cohort_values = {r["cohort_id"] for r in summary}
        if len(species_values) != 1:
            raise ValueError(f"Species mismatch: {summary_file}")
        species = species_values.pop()
        if source_cohorts is None:
            source_cohort = cohort
        else:
            if species not in source_cohorts:
                raise ValueError(f"Unregistered species in study: {species}")
            source_cohort = source_cohorts[species]
        if cohort_values != {source_cohort}:
            raise ValueError(f"Species or cohort mismatch: {summary_file}")
        if not species or species in seen_species or folder.name != source_cohort or folder.parent.name != species:
            raise ValueError(f"Duplicate or inconsistent species catalogue: {folder}")
        seen_species.add(species)

        families = tables["candidate_families.tsv"]
        ids = {(r["caller"], r["family_id"]) for r in families}
        if len(ids) != len(families):
            raise ValueError(f"Duplicate candidate family ID: {folder}")
        if any((r["species_id"], r["cohort_id"]) != (species, source_cohort)
               for r in families + tables["genome_sizes.tsv"]):
            raise ValueError(f"Invalid row identity: {folder}")
        for name in ("call_family_map.tsv", "tip_states.tsv"):
            if any((r["caller"], r["family_id"]) not in ids for r in tables[name]):
                raise ValueError(f"Orphan family in {name}: {folder}")
        genome_ids = {r["assembly_id"] for r in tables["genome_sizes.tsv"]}
        if len(genome_ids) != len(tables["genome_sizes.tsv"]):
            raise ValueError(f"Duplicate genome: {folder}")
        if any(r["assembly_id"] not in genome_ids or r["state"] not in {"0", "1", "?"}
               for r in tables["tip_states.tsv"]):
            raise ValueError(f"Invalid tip states: {folder}")
        if len(tables["tip_states.tsv"]) != len(ids) * len(genome_ids):
            raise ValueError(f"Missing tip states: {folder}")
        states = {(r["caller"], r["family_id"], r["assembly_id"])
                  for r in tables["tip_states.tsv"]}
        if len(states) != len(tables["tip_states.tsv"]):
            raise ValueError(f"Duplicate tip states: {folder}")
        for row in tables["cross_caller_associations.tsv"]:
            if ("leen", row["leen_family_id"]) not in ids or (
                "syri", row["syri_family_id"]) not in ids:
                raise ValueError(f"Orphan cross-caller association: {folder}")
        for row in summary:
            caller = row["caller"]
            if int(row["total_genomes"]) != len(genome_ids):
                raise ValueError(f"Genome count mismatch: {folder}/{caller}")
            if int(row["candidate_families"]) != sum(
                r["caller"] == caller for r in families
            ):
                raise ValueError(f"Family count mismatch: {folder}/{caller}")
            for field in (
                "total_genomes", "callable_genomes", "unknown_genomes",
                "candidate_families", "families_with_multiple_carriers", "calls",
                "complex_geometry_calls",
            ):
                counts[(caller, field)] += int(row[field])
            counts[(caller, "species")] += 1
        counts[("both", "associations")] += len(tables["cross_caller_associations.tsv"])
        counts[("both", "unique")] += sum(
            r["association_status"] == "unique_candidate"
            for r in tables["cross_caller_associations.tsv"]
        )
        counts[("both", "ambiguous")] += sum(
            r["association_status"] == "ambiguous_candidate"
            for r in tables["cross_caller_associations.tsv"]
        )

        for name, rows in tables.items():
            output_name = "species_summary.tsv" if name == "summary.tsv" else name
            for row in rows:
                result = dict(row)
                if name in ("call_family_map.tsv", "tip_states.tsv", "cross_caller_associations.tsv"):
                    result.update(species_id=species, cohort_id=cohort)
                elif source_cohorts is not None:
                    result["cohort_id"] = cohort
                by_name[output_name].append(result)

    total_pairs = counts[("both", "associations")]
    if total_pairs != counts[("both", "unique")] + counts[("both", "ambiguous")]:
        raise ValueError("Cross-caller association statuses are inconsistent")
    overall = []
    for caller in ("leen", "syri"):
        overall.append({
            "cohort_id": cohort,
            "caller": caller,
            "species_completed": counts[(caller, "species")],
            "genomes_total_across_species": counts[(caller, "total_genomes")],
            "callable_genomes_total": counts[(caller, "callable_genomes")],
            "unknown_genomes_total": counts[(caller, "unknown_genomes")],
            "candidate_families_total": counts[(caller, "candidate_families")],
            "families_with_multiple_carriers_total": counts[(caller, "families_with_multiple_carriers")],
            "calls_total": counts[(caller, "calls")],
            "complex_geometry_calls_total": counts[(caller, "complex_geometry_calls")],
            "provisional_cross_caller_pairs_total": total_pairs,
            "unique_cross_caller_pairs_total": counts[("both", "unique")],
            "ambiguous_cross_caller_pairs_total": counts[("both", "ambiguous")],
            "interpretation": "caller_specific_candidate_families_not_independent_events",
        })

    # Validate all data before writing aggregated output files.
    outdir = Path(outdir)
    for name, fields in OUTPUTS.items():
        write_table(outdir / name, fields, by_name[name])
    write_table(outdir / "run_summary.tsv", RUN_FIELDS, overall)
    print(f"Module 7: {len(seen_species)} species, "
          f"{len(by_name['candidate_families.tsv'])} caller-specific candidate families, "
          f"{total_pairs} provisional cross-caller pairs")
    return overall


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summaries", nargs="+", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cohort", required=True)
    args = parser.parse_args()
    aggregate(args.summaries, args.outdir, args.cohort)


if __name__ == "__main__":
    main()
