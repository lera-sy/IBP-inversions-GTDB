"""Module 7: strict, species-aware aggregation of Module 6 catalogues."""
from pathlib import Path

MODULE7_RESULTS = str(config.get(
    "module7_results_root", Path(DOWNLOAD_DIR) / "module7"
)).rstrip("/")
MODULE7_OUTDIR = f"{MODULE7_RESULTS}/{COHORT}"

rule module7_complete:
    input:
        f"{MODULE7_OUTDIR}/run_summary.tsv"

rule module7_aggregate:
    input:
        reports=module6_complete_targets,
        code="src/ibp/module7/aggregate.py",
    output:
        summary=f"{MODULE7_OUTDIR}/run_summary.tsv",
        species=f"{MODULE7_OUTDIR}/species_summary.tsv",
        families=f"{MODULE7_OUTDIR}/candidate_families.tsv",
        mapping=f"{MODULE7_OUTDIR}/call_family_map.tsv",
        states=f"{MODULE7_OUTDIR}/tip_states.tsv",
        associations=f"{MODULE7_OUTDIR}/cross_caller_associations.tsv",
        genomes=f"{MODULE7_OUTDIR}/genome_sizes.tsv",
    params:
        outdir=MODULE7_OUTDIR,
    shell:
        """
        PYTHONPATH=src python -m ibp.module7.aggregate \
          --summaries {input.reports:q} \
          --cohort {COHORT:q} --outdir {params.outdir:q}
        """
