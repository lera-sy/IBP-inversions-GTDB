#!/bin/bash
# Phase 4 pilot: Leen + minimap2/SyRI on the pilot species, then the
# caller-specific catalogue (module 6). Uses the packages made by
# hpc/prepare_chromosomes.slurm; nothing upstream is rerun.
#
# Submit from the repo folder:
#   sbatch --account=lp_edu_large_omics hpc/pilot_callers.slurm
#
#SBATCH --job-name=ibp-pilot
#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/pilot_%j.log

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$VSC_DATA/miniforge3/bin/activate"
conda activate ibp-pipeline

echo "Start: $(date)  cores: $SLURM_CPUS_PER_TASK"
snakemake --cores "$SLURM_CPUS_PER_TASK" --configfile config/runs/pilot.yaml \
    --rerun-incomplete --keep-going module6_complete
echo "End: $(date)"
