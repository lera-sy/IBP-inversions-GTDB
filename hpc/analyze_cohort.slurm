#!/bin/bash
# Full analysis of one cohort after its callers have finished, in one job:
# loci/tiers -> models + room-alone null -> spanning/local split -> local checks -> figures
# -> robustness -> class tests -> distribution + ori/ter checks -> poster figures.
# Every step prints a header; the log is logs/analysis_<cohort>_<jobid>.log.
#   COHORT=selected_min20_v1 sbatch --account=lp_edu_large_omics --export=ALL hpc/analyze_cohort.slurm
#
#SBATCH --job-name=ibp-analysis
#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/analysis_%j.log

source "$VSC_DATA/miniforge3/bin/activate"
conda activate ibp-pipeline
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src
export COHORT=${COHORT:-selected_min20_v1}
if [ "$COHORT" = selected_50_v1 ]; then SEL=data/processed; else SEL=data/processed/runs/$COHORT; fi
step() { echo; echo "########## $1   ($(date '+%H:%M'))"; }
echo "Cohort: $COHORT"

step "1. Loci, families, tiers";            bash hpc/build_loci.sh
step "2. Room-alone null and models";       bash hpc/run_model.sh
step "3. Spanning vs local, asked separately"; bash hpc/split_classes.sh
step "4. Local-inversion checks";           bash hpc/check_local.sh | head -n 30
step "5. Figures";                          bash hpc/plot_run.sh
step "6. Robustness (terminus, nulls, ancestry, percentiles)"; bash hpc/robustness.sh
step "7. Class tests (interaction, empirical nulls, enrichment, symmetry)"
python -u -m ibp.stats.class_tests --loci "results/loci/$COHORT" --tree "$SEL/species.tree" \
  --taxonomy "$SEL/species.tsv" --termini "results/robustness/$COHORT/termini_gc_skew.tsv" \
  --replicates 200 --workers "${SLURM_CPUS_PER_TASK:-1}" --out "results/class_tests/$COHORT"
step "8. Distribution and ori vs ter";      bash hpc/extra_checks.sh
step "9. Poster figures"
python -m ibp.stats.poster_figures --loci "results/loci/$COHORT" --model "results/model/$COHORT" \
  --split "results/figures/$COHORT/split" --class-tests "results/class_tests/$COHORT" --out "results/poster/$COHORT"
step "Done"
