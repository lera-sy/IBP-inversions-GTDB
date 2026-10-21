"""Descriptive inversion-only analyses and figures for Module 7 outputs.

No candidate family or cross-caller association is an inferred evolutionary event.
"""

import argparse
from collections import defaultdict
from pathlib import Path
from statistics import median

from ibp.module7.aggregate import OUTPUTS, read_table, write_table
from ibp.module8.agreement import (
    AGREEMENT_FIELDS,
    CALL_FIELDS,
    FAMILY_SUPPORT_FIELDS,
    build_agreement,
    make_agreement_figures,
)

SPECIES_METRICS = [
    "species_id", "cohort_id", "caller", "genomes", "callable_genomes",
    "unknown_genomes", "median_chromosome_size_bp", "median_genome_size_bp",
    "candidate_families", "families_per_callable_genome",
    "median_candidate_length_bp", "median_family_to_chromosome_ratio",
    "interpretation",
]
ANALYSIS_FIELDS = [
    "caller", "species_with_results", "total_candidate_families",
    "median_candidate_length_bp", "median_chromosome_size_bp_across_species",
    "species_with_size_observations", "interpretation",
]


def build_statistics(indir):
    indir = Path(indir)
    summaries = read_table(indir / "species_summary.tsv", OUTPUTS["species_summary.tsv"])
    families = read_table(indir / "candidate_families.tsv", OUTPUTS["candidate_families.tsv"])
    genomes = read_table(indir / "genome_sizes.tsv", OUTPUTS["genome_sizes.tsv"])
    if not summaries:
        raise ValueError("No species summaries to analyze")
    grouped_genomes = defaultdict(list)
    grouped_families = defaultdict(list)
    for row in genomes:
        grouped_genomes[(row["species_id"], row["cohort_id"])].append(row)
    for row in families:
        grouped_families[(row["species_id"], row["cohort_id"], row["caller"])].append(row)
    metrics = []
    for summary in summaries:
        species, cohort, caller = (summary[k] for k in ("species_id", "cohort_id", "caller"))
        ids = (species, cohort)
        gs = grouped_genomes[ids]
        chrom = [int(g["chromosome_length_bp"]) for g in gs]
        genomes_bp = [int(g["genome_length_bp"]) for g in gs]
        if not chrom or any(x <= 0 for x in chrom + genomes_bp):
            raise ValueError(f"Invalid/missing genome lengths for {ids}")
        med_chrom = median(chrom)
        med_genome = median(genomes_bp)
        fs = grouped_families[(species, cohort, caller)]
        lengths = [float(f["median_reported_length_bp"]) for f in fs]
        if any(x <= 0 for x in lengths):
            raise ValueError(f"Invalid family length in {species}/{caller}")
        observed_n = int(summary["candidate_families"])
        callable_n = int(summary["callable_genomes"])
        if observed_n != len(lengths) or not 0 < callable_n <= len(gs):
            raise ValueError(f"Summary inconsistency in {species}/{caller}")
        metrics.append({
            "species_id": species, "cohort_id": cohort, "caller": caller,
            "genomes": len(gs), "callable_genomes": callable_n,
            "unknown_genomes": int(summary["unknown_genomes"]),
            "median_chromosome_size_bp": med_chrom,
            "median_genome_size_bp": med_genome,
            "candidate_families": observed_n,
            "families_per_callable_genome": f"{observed_n / callable_n:.9f}",
            "median_candidate_length_bp": median(lengths) if lengths else "",
            "median_family_to_chromosome_ratio": (
                f"{median(lengths) / med_chrom:.9f}" if lengths else ""
            ),
            "interpretation": "descriptive_candidate_families_not_independent_events",
        })
    analyses = []
    for caller in ("leen", "syri"):
        caller_rows = [r for r in metrics if r["caller"] == caller]
        lengths = [float(f["median_reported_length_bp"]) for f in families if f["caller"] == caller]
        analyses.append({
            "caller": caller,
            "species_with_results": len(caller_rows),
            "total_candidate_families": len(lengths),
            "median_candidate_length_bp": median(lengths) if lengths else "",
            "median_chromosome_size_bp_across_species": median(
                float(x["median_chromosome_size_bp"]) for x in caller_rows
            ) if caller_rows else "",
            "species_with_size_observations": sum(
                x["median_candidate_length_bp"] != "" for x in caller_rows
            ),
            "interpretation": "tree_free_descriptive_analysis_only",
        })
    return metrics, analyses, families


def make_figures(metrics, families, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = Path(outdir) / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    labels = ("leen", "syri")
    # Reuse the legacy inversion-family size comparisons, but from standardized data.
    fig, ax = plt.subplots(figsize=(8, 5))
    for caller in labels:
        values = [float(r["median_reported_length_bp"]) / 1000
                  for r in families if r["caller"] == caller]
        if values:
            ax.hist(values, bins=min(18, max(1, len(values))), alpha=0.5,
                    label=f"{caller} (n={len(values)})")
    ax.set(xlabel="Candidate inversion-family length (kb)", ylabel="Number of families",
           title="Caller-specific candidate inversion lengths")
    ax.legend(loc="best") if ax.has_data() else None
    fig.tight_layout()
    fig.savefig(figures / "family_lengths.png", dpi=175)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for caller in labels:
        rows = [r for r in metrics if r["caller"] == caller and r["median_candidate_length_bp"] != ""]
        if rows:
            ax.scatter([float(r["median_genome_size_bp"]) / 1_000_000 for r in rows],
                       [float(r["median_candidate_length_bp"]) / 1000 for r in rows],
                       label=f"{caller} (n={len(rows)})", alpha=0.8)
    ax.set(xlabel="Median genome length (Mb; accepted assemblies)",
           ylabel="Median candidate-family length (kb)",
           title="Genome size vs candidate inversion-family length")
    ax.legend(loc="best") if ax.has_data() else None
    fig.tight_layout()
    fig.savefig(figures / "genome_vs_family_sizes.png", dpi=175)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for caller in labels:
        rows = [r for r in metrics if r["caller"] == caller]
        if rows:
            ax.scatter([float(r["median_genome_size_bp"]) / 1_000_000 for r in rows],
                       [float(r["families_per_callable_genome"]) for r in rows],
                       label=f"{caller} (n={len(rows)})", alpha=0.8)
    ax.set(xlabel="Median genome length (Mb; accepted assemblies)",
           ylabel="Candidate families / callable genomes",
           title="Candidate inversion-family burden across species")
    ax.legend(loc="best") if ax.has_data() else None
    fig.tight_layout()
    fig.savefig(figures / "families_per_callable_genome.png", dpi=175)
    plt.close(fig)


def analyze(indir, outdir):
    metrics, analyses, families = build_statistics(indir)
    outdir = Path(outdir)
    make_figures(metrics, families, outdir)
    agreement, agreement_calls, agreement_families = build_agreement(indir)
    make_agreement_figures(agreement, agreement_families, outdir)
    write_table(outdir / "caller_agreement.tsv", AGREEMENT_FIELDS, agreement)
    write_table(outdir / "caller_agreement_calls.tsv", CALL_FIELDS, agreement_calls)
    write_table(outdir / "caller_agreement_families.tsv", FAMILY_SUPPORT_FIELDS, agreement_families)
    write_table(outdir / "species_metrics.tsv", SPECIES_METRICS, metrics)
    write_table(outdir / "analysis_summary.tsv", ANALYSIS_FIELDS, analyses)
    print(f"Module 8: {len(metrics)//2} species; five descriptive figures written to {outdir/'figures'}")
    return analyses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aggregate-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.aggregate_dir, args.outdir)


if __name__ == "__main__":
    main()
