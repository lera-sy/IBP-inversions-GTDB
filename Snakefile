import os

import pandas as pd

# All definitions and settings live in config/config.yaml.
# Run-specific overrides: --configfile config/runs/<run>.yaml
configfile: "config/config.yaml"

envvars:
    "VSC_SCRATCH"

# Isolated runs keep Modules 2-8 separate while sharing Module 1's metadata.
# This mode is optional so all previously validated cohorts keep their paths.
COHORT = str(config.get("cohort", "selected_50_v1"))
if not COHORT or "/" in COHORT or "\\" in COHORT or COHORT in {".", ".."}:
    raise ValueError("cohort must be a single directory name.")
ISOLATE_RUN = str(config.get("isolate_run", False)).lower() in {"true", "1", "yes"}
ROOT_DOWNLOAD_DIR = os.path.join(os.environ["VSC_SCRATCH"], "IBP-main")
DOWNLOAD_DIR = (
    os.path.join(ROOT_DOWNLOAD_DIR, "runs", COHORT)
    if ISOLATE_RUN else ROOT_DOWNLOAD_DIR
)
RUN_PROCESSED_DIR = (
    f"data/processed/runs/{COHORT}" if ISOLATE_RUN else "data/processed"
)
REFERENCE_DIR = DOWNLOAD_DIR + "/reference_selection/{species}/" + COHORT
PREPARED_DIR = DOWNLOAD_DIR + "/prepared/{species}/" + COHORT
# Written by Module 1 (GTDB selection); taxid = GTDB species key, not an NCBI TaxID.
SPECIES_TABLE = "data/processed/species.tsv"
# Normal runs use Module 1's full multi-species assembly selection.
# Only an explicitly configured local integration test filters downstream assemblies.
LOCAL_SELECTION_TEST = "run_assembly_limit" in config
RUN_ASSEMBLIES = (
    f"{RUN_PROCESSED_DIR}/run_selected_assemblies.tsv"
    if LOCAL_SELECTION_TEST
    else "data/processed/selected_assemblies.tsv"
)


def selected_species():
    result = checkpoints.normalize_chromosomes.get()
    metadata = pd.read_csv(result.output.passed, sep="\t", dtype=str)
    available = sorted(metadata["Organism Taxonomic ID"].unique())
    if "species" in config and "species_file" in config:
        raise ValueError("Use species or species_file, not both.")
    species = config.get("species", available)
    if "species_file" in config:
        species = pd.read_csv(config["species_file"], sep="\t", dtype=str)["species_id"].tolist()
    if not isinstance(species, list) or not species:
        raise ValueError("species must be a nonempty list of species IDs or TaxIDs.")
    names = pd.read_csv(SPECIES_TABLE, sep="\t", dtype=str).set_index("species_id")["taxid"].to_dict()
    species = [names.get(str(identifier), str(identifier)) for identifier in species]
    if set(species) - set(available):
        raise ValueError("Requested species has no accepted chromosomes.")
    return species


def reference_targets(wildcards):
    return expand(REFERENCE_DIR + "/reference.txt", species=selected_species())


def prepared_targets(wildcards):
    species = selected_species()
    names = pd.read_csv(SPECIES_TABLE, sep="\t", dtype=str).set_index("taxid")
    if names.index.duplicated().any() or set(species) - set(names.index):
        raise ValueError("Each requested TaxID needs one entry in config/species.tsv.")
    return expand(PREPARED_DIR, species=names.loc[species, "species_id"].tolist())


def reference_fastas(wildcards):
    result = checkpoints.normalize_chromosomes.get()
    metadata = pd.read_csv(result.output.passed, sep="\t", dtype=str)
    taxid = wildcards.species.rsplit("_", 1)[-1]
    selected = metadata.loc[metadata["Organism Taxonomic ID"] == taxid, "Assembly Accession"]
    if selected.empty:
        raise ValueError(f"No accepted chromosomes for species {wildcards.species}.")
    return [os.path.join(result.output.normalized, accession + ".fna") for accession in sorted(selected)]

rule all:
    input:
        "data/processed/selected_assemblies.tsv",
        "data/processed/selection_summary.tsv",
        "data/processed/species.tree",
        DOWNLOAD_DIR + "/genomes",
        DOWNLOAD_DIR + "/sequence_reports",
        f"{RUN_PROCESSED_DIR}/genome_qc.tsv",
        f"{RUN_PROCESSED_DIR}/qc_passed_assemblies.tsv",
        DOWNLOAD_DIR + "/chromosomes",
        f"{RUN_PROCESSED_DIR}/chromosome_topology.tsv",
        f"{RUN_PROCESSED_DIR}/normalization.tsv",
        f"{RUN_PROCESSED_DIR}/normalized_assemblies.tsv",
        DOWNLOAD_DIR + "/normalized_chromosomes",
        prepared_targets,
        lambda wildcards: module5_complete_targets(wildcards)

# Module 1 is kept separate from downstream download and preparation rules.
include: "workflow/rules/module1.smk"

# The small-cohort selector belongs to the integration test, not production.
if LOCAL_SELECTION_TEST:
    include: "workflow/testing/local_species_selection.smk"


# Modules 2–4: execution and preparation stages.
include: "workflow/rules/module2.smk"
include: "workflow/rules/module3.smk"
include: "workflow/rules/module4.smk"

# Module 5 extends the original workflow without reorganizing Modules 1-4.
# Its input functions resolve accepted species and genomes after checkpoints.
include: "workflow/rules/module5.smk"
include: "workflow/rules/module6.smk"
include: "workflow/rules/module7.smk"
include: "workflow/rules/module8.smk"
