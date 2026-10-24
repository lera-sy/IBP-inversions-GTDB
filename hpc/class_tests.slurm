#!/bin/bash
# Spanning vs local, formal tests (interaction test, empirical p vs both nulls with the
# same model, spanning share with clustered CIs, symmetry). Too heavy for the login node
# (CPU time limit). Needs results/robustness/<cohort>/termini_gc_skew.tsv.
# Submit from the repo folder:
#   sbatch --account=lp_edu_large_omics hpc/class_tests.slurm
#   COHORT=selected_min20_v1 sbatch --account=lp_edu_large_omics --export=ALL hpc/class_tests.slurm
#
#SBATCH --job-name=ibp-classtests
#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/class_tests_%j.log

source "$VSC_DATA/miniforge3/bin/activate"
conda activate ibp-pipeline
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
COHORT=${COHORT:-selected_50_v1}
if [ "$COHORT" = selected_50_v1 ]; then SEL=data/processed; else SEL=data/processed/runs/$COHORT; fi
mkdir -p logs
PYTHONPATH=src python -u -m ibp.stats.class_tests --loci "results/loci/$COHORT" --tree "$SEL/species.tree" \
  --taxonomy "$SEL/species.tsv" --termini "results/robustness/$COHORT/termini_gc_skew.tsv" \
  --replicates 200 --workers "${SLURM_CPUS_PER_TASK:-1}" --out "results/class_tests/$COHORT"
