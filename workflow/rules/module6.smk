"""Module 6: tree-free inversion families and provisional caller associations.

Independent evolutionary event inference is NOT performed here.
Module 5 comparison statuses are never converted to zero calls on failure.
"""

from pathlib import Path

if MODULE5_INTEGRATED:
    MODULE6_RESULTS = str(config.get("module6_results_root", Path(DOWNLOAD_DIR) / "module6")).rstrip("/")
else:
    MODULE6_RESULTS = str(config.get("module6_results_root", Path(RESULTS).parent / "module6")).rstrip("/")


def module6_complete_targets(wc):
    # Reuse Module 5's checkpoint-aware species/cohort discovery.
    reports = module5_complete_targets(wc)
    return [
        str(Path(MODULE6_RESULTS) / Path(str(report)).relative_to(RESULTS).parent / "summary.tsv")
        for report in reports
    ]


rule module6_complete:
    input:
        module6_complete_targets


def module6_package_file(wc, name):
    """Resolve published Module 4 files only after its checkpoint completes.

    A literal manifest path is not sufficient to build the DAG: the manifest is
    inside a checkpoint-produced package directory, so it has no direct rule.
    Module 5 already supplies the checkpoint-aware package resolver.
    """
    return str(Path(module5_package(wc)) / name)


rule module6_catalogue:
    input:
        report=f"{RESULTS}/{{species_id}}/{{cohort_id}}/comparisons.tsv",
        package=module5_package,
        manifest=lambda wc: module6_package_file(wc, "manifest.json"),
        inventory=lambda wc: module6_package_file(wc, "genomes.tsv"),
        code="src/ibp/module6/catalog.py",
    output:
        families=f"{MODULE6_RESULTS}/{{species_id}}/{{cohort_id}}/candidate_families.tsv",
        mapping=f"{MODULE6_RESULTS}/{{species_id}}/{{cohort_id}}/call_family_map.tsv",
        states=f"{MODULE6_RESULTS}/{{species_id}}/{{cohort_id}}/tip_states.tsv",
        associations=f"{MODULE6_RESULTS}/{{species_id}}/{{cohort_id}}/cross_caller_associations.tsv",
        genomes=f"{MODULE6_RESULTS}/{{species_id}}/{{cohort_id}}/genome_sizes.tsv",
        summary=f"{MODULE6_RESULTS}/{{species_id}}/{{cohort_id}}/summary.tsv",
    params:
        outdir=lambda wc: f"{MODULE6_RESULTS}/{wc.species_id}/{wc.cohort_id}",
        tolerance=lambda wc: int(config.get("module6", {}).get("breakpoint_tolerance_bp", 5000)),
    shell:
        """
        PYTHONPATH=src python -m ibp.module6.catalog \
            --manifest {input.manifest:q} \
            --report {input.report:q} \
            --outdir {params.outdir:q} \
            --breakpoint-tolerance {params.tolerance}
        """
