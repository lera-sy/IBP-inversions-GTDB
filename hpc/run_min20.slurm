#!/bin/bash
# Min-20 run, start to finish in one job: GTDB selection (>= 20 genomes),
# download, QC, topology, normalization (reusing the min-50 chromosomes),
# references, packages, both callers, module 6. Own folders; the min-50
# results are not touched. Resubmitting resumes where it stopped.
#
# Submit from the repo folder:
#   sbatch --account=lp_edu_large_omics hpc/run_min20.slurm
#
#SBATCH --job-name=ibp-min20
#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/min20_%j.log

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$VSC_DATA/miniforge3/bin/activate"
conda activate ibp-pipeline

MEM_MB=$(( SLURM_MEM_PER_NODE * 85 / 100 ))
echo "Start: $(date)  cores: $SLURM_CPUS_PER_TASK  memory budget: ${MEM_MB} MB"
snakemake module6_complete data/processed/runs/selected_min20_v1/species.tree \
    --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb="$MEM_MB" \
    --retries 1 --rerun-incomplete --keep-going \
    --configfile config/runs/min20.yaml
echo "End: $(date)"
