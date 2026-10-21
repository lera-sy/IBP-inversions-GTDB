"""Module 5: inversion calling in standalone or upstream-integrated workflows.

Standalone mode uses a config with `packages` (e.g. the mock fixture).
Integrated mode uses the normalization and species-preparation checkpoints.
"""

import shlex
import sys
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path("src").resolve()))

from ibp.module5.validate_input import load_package
from ibp.module5.plan_comparisons import plan_jobs

MODULE5_INTEGRATED = "packages" not in config

if MODULE5_INTEGRATED:
    # These variables/functions are supplied by the root Snakefile.
    PREPARED = Path(DOWNLOAD_DIR) / "prepared"
    RESULTS = str(config.get("module5_results_root", Path(DOWNLOAD_DIR) / "module5")).rstrip("/")
else:
    PREPARED = Path(config["prepared_root"])
    RESULTS = str(config["results_root"]).rstrip("/")


@lru_cache(maxsize=None)
def module5_jobs(species, cohort):
    """Validate each published package once and plan both callers."""
    manifest = PREPARED / species / cohort / "manifest.json"
    package = load_package(manifest)
    if package.species_id != species or package.cohort_id != cohort:
        raise ValueError(f"Package identity mismatch: {manifest}")
    return plan_jobs(package)


def module5_package(wc):
    """Wait for the Module 4 package in integrated mode."""
    if not MODULE5_INTEGRATED:
        return str(PREPARED / wc.species_id / wc.cohort_id)
    if wc.cohort_id != COHORT:
        raise ValueError(f"Unexpected cohort {wc.cohort_id}; expected {COHORT}")
    return str(checkpoints.prepare_species.get(species=wc.species_id).output.package)


def module5_job(wc):
    if MODULE5_INTEGRATED:
        module5_package(wc)  # Trigger checkpoint reevaluation if necessary.
    matches = [
        job for job in module5_jobs(wc.species_id, wc.cohort_id)
        if job.caller == wc.caller and job.query_assembly == wc.query_assembly
    ]
    if len(matches) != 1:
        raise ValueError(f"Unknown or duplicate Module 5 comparison: {wc}")
    return matches[0]


def module5_options(wc):
    if wc.caller == "leen":
        options = config.get("leen", {})
        return (
            f"--window {int(options.get('window', 2000))} "
            f"--min-size {int(options.get('min_size', 6000))} "
            f"--max-jump {int(options.get('max_jump', 20000))}"
        )
    if wc.caller == "syri":
        environment = config.get("syri_env", "ibp-syri")
        return "--syri-env " + shlex.quote(str(environment)) if environment else ""
    raise ValueError(f"Unknown caller: {wc.caller}")


def module5_report_inputs(wc):
    """Discover comparison files only after Module 4 publishes the package."""
    module5_package(wc)
    return [
        f"{RESULTS}/{wc.species_id}/{wc.cohort_id}/{job.caller}/"
        f"{job.query_assembly}/{filename}"
        for job in module5_jobs(wc.species_id, wc.cohort_id)
        for filename in ("status.json", "calls.tsv", "inversions.tsv")
    ]


def module5_complete_targets(wc):
    if MODULE5_INTEGRATED:
        # selected_species() depends on the normalization checkpoint.
        taxids = selected_species()
        names = pd.read_csv(SPECIES_TABLE, sep="\t", dtype=str)
        by_taxid = names.set_index("taxid")
        if by_taxid.index.duplicated().any() or set(taxids) - set(by_taxid.index):
            raise ValueError("Missing or duplicate TaxIDs in data/processed/species.tsv")
        species_ids = by_taxid.loc[taxids, "species_id"].tolist()
        return [f"{RESULTS}/{species}/{COHORT}/comparisons.tsv" for species in species_ids]
    return [
        f"{RESULTS}/{item['species_id']}/{item['cohort_id']}/comparisons.tsv"
        for item in config["packages"]
    ]


rule module5_complete:
    input:
        module5_complete_targets


rule inversion_pair:
    input:
        # Directory is an explicit dependency on the upstream package.
        package=module5_package,
        manifest=lambda wc: str(PREPARED / wc.species_id / wc.cohort_id / "manifest.json"),
        inventory=lambda wc: str(PREPARED / wc.species_id / wc.cohort_id / "genomes.tsv"),
        reference=lambda wc: str(module5_job(wc).reference_fasta),
        query=lambda wc: str(module5_job(wc).query_fasta),
    output:
        status=f"{RESULTS}/{{species_id}}/{{cohort_id}}/{{caller}}/{{query_assembly}}/status.json",
        calls=f"{RESULTS}/{{species_id}}/{{cohort_id}}/{{caller}}/{{query_assembly}}/calls.tsv",
        standardized=f"{RESULTS}/{{species_id}}/{{cohort_id}}/{{caller}}/{{query_assembly}}/inversions.tsv",
    params:
        outdir=lambda wc: f"{RESULTS}/{wc.species_id}/{wc.cohort_id}/{wc.caller}/{wc.query_assembly}",
        reference_id=lambda wc: module5_job(wc).reference_assembly,
        reference_topology=lambda wc: module5_job(wc).reference_topology,
        query_topology=lambda wc: module5_job(wc).query_topology,
        options=module5_options,
    threads: 2
    shell:
        """
        PYTHONPATH=src python -m ibp.module5.{wildcards.caller} \
            --reference {input.reference:q} \
            --query {input.query:q} \
            --reference-id {params.reference_id:q} \
            --query-id {wildcards.query_assembly:q} \
            --reference-topology {params.reference_topology:q} \
            --query-topology {params.query_topology:q} \
            --outdir {params.outdir:q} \
            --threads {threads} {params.options}
        """


rule module5_cohort_report:
    input:
        module5_report_inputs
    output:
        f"{RESULTS}/{{species_id}}/{{cohort_id}}/comparisons.tsv"
    run:
        from ibp.module5.results import write_comparison_report

        comparisons = [
            (
                wildcards.species_id, wildcards.cohort_id, job.caller,
                job.reference_assembly, job.query_assembly,
                Path(RESULTS) / wildcards.species_id / wildcards.cohort_id
                / job.caller / job.query_assembly,
            )
            for job in module5_jobs(wildcards.species_id, wildcards.cohort_id)
        ]
        write_comparison_report(output[0], comparisons)
