"""Module 8: descriptive statistics and publication-oriented figure outputs."""
from pathlib import Path

MODULE8_RESULTS = str(config.get(
    "module8_results_root", Path(DOWNLOAD_DIR) / "module8"
)).rstrip("/")
MODULE8_OUTDIR = f"{MODULE8_RESULTS}/{COHORT}"

rule module8_complete:
    input:
        f"{MODULE8_OUTDIR}/analysis_summary.tsv",
        f"{MODULE8_OUTDIR}/species_metrics.tsv",
        f"{MODULE8_OUTDIR}/caller_agreement.tsv",
        f"{MODULE8_OUTDIR}/caller_agreement_calls.tsv",
        f"{MODULE8_OUTDIR}/caller_agreement_families.tsv",
        f"{MODULE8_OUTDIR}/figures/family_lengths.png",
        f"{MODULE8_OUTDIR}/figures/genome_vs_family_sizes.png",
        f"{MODULE8_OUTDIR}/figures/families_per_callable_genome.png",
        f"{MODULE8_OUTDIR}/figures/caller_agreement.png",
        f"{MODULE8_OUTDIR}/figures/agreement_family_lengths.png",

rule module8_analyze:
    input:
        summary=f"{MODULE7_OUTDIR}/run_summary.tsv",
        species=f"{MODULE7_OUTDIR}/species_summary.tsv",
        families=f"{MODULE7_OUTDIR}/candidate_families.tsv",
        genomes=f"{MODULE7_OUTDIR}/genome_sizes.tsv",
        mapping=f"{MODULE7_OUTDIR}/call_family_map.tsv",
        associations=f"{MODULE7_OUTDIR}/cross_caller_associations.tsv",
        agreement_code="src/ibp/module8/agreement.py",
        code="src/ibp/module8/analyze.py",
    output:
        summary=f"{MODULE8_OUTDIR}/analysis_summary.tsv",
        metrics=f"{MODULE8_OUTDIR}/species_metrics.tsv",
        agree=f"{MODULE8_OUTDIR}/caller_agreement.tsv",
        calls=f"{MODULE8_OUTDIR}/caller_agreement_calls.tsv",
        support=f"{MODULE8_OUTDIR}/caller_agreement_families.tsv",
        fig1=f"{MODULE8_OUTDIR}/figures/family_lengths.png",
        fig2=f"{MODULE8_OUTDIR}/figures/genome_vs_family_sizes.png",
        fig3=f"{MODULE8_OUTDIR}/figures/families_per_callable_genome.png",
        fig4=f"{MODULE8_OUTDIR}/figures/caller_agreement.png",
        fig5=f"{MODULE8_OUTDIR}/figures/agreement_family_lengths.png",
    params:
        root=MODULE7_OUTDIR,
        outdir=MODULE8_OUTDIR,
    shell:
        """
        PYTHONPATH=src python -m ibp.module8.analyze \
           --aggregate-dir {params.root:q} --outdir {params.outdir:q}
        """
