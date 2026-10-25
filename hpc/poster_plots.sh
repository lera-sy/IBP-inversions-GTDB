#!/bin/bash
# Poster figures: question.png, split_boxes.png, split_scatter.png (seconds; login node is fine).
# Needs run_model.sh, split_classes.sh and class_tests.slurm to have run for the cohort.
#   bash hpc/poster_plots.sh
#   COHORT=selected_min20_v1 bash hpc/poster_plots.sh
if [ "${CONDA_DEFAULT_ENV:-}" != ibp-pipeline ]; then source "$VSC_DATA/miniforge3/bin/activate"; conda activate ibp-pipeline; fi
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
PYTHONPATH=src python -m ibp.stats.poster_plots --loci "results/loci/$COHORT" --model "results/model/$COHORT" \
  --split "results/figures/$COHORT/split" --class-tests "results/class_tests/$COHORT" --out "results/poster/$COHORT"
