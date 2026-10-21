#!/bin/bash
# Phase 2 -> 4 preparation as one batch job: downloads (reused if done), QC,
# topology, dnaA normalization, Mash reference and species packages.
# Normalization runs one single-thread Dnaapler per core, so cores = speed.
#
# Submit from the repo folder:
#   sbatch --account=<your lp_ account> hpc/prepare_chromosomes.slurm
#
#SBATCH --job-name=ibp-prepare
#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=72G
#SBATCH --time=04:00:00
#SBATCH --output=logs/prepare_%j.log

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
source "$VSC_DATA/miniforge3/bin/activate"
conda activate ibp-pipeline

echo "Start: $(date)  cores: $SLURM_CPUS_PER_TASK"
snakemake --cores "$SLURM_CPUS_PER_TASK" --rerun-incomplete --keep-going prepare_packages
echo "End: $(date)"
