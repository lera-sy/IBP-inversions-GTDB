#!/bin/bash
# Phase 5: Leen + minimap2/SyRI on ALL selected species, then the
# caller-specific catalogue (module 6). Uses the packages made by
# hpc/prepare_chromosomes.slurm; nothing upstream is rerun.
#
# Submit from the repo folder:
#   sbatch --account=lp_edu_large_omics hpc/run_callers.slurm
#
#SBATCH --job-name=ibp-callers
#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=72G
#SBATCH --time=12:00:00
#SBATCH --output=logs/callers_%j.log

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source "$VSC_DATA/miniforge3/bin/activate"
conda activate ibp-pipeline

echo "Start: $(date)  cores: $SLURM_CPUS_PER_TASK"
snakemake --cores "$SLURM_CPUS_PER_TASK" \
    --rerun-incomplete --keep-going module6_complete
echo "End: $(date)"
