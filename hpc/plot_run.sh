#!/bin/bash
# Figures for one run only (no team comparison). ibp-pipeline active:
#   bash hpc/plot_run.sh                       # min-50 run
#   COHORT=selected_min20_v1 bash hpc/plot_run.sh
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
PYTHONPATH=src python -m ibp.stats.plot_run --loci "results/loci/$COHORT" --model "results/model/$COHORT" --out "results/figures/$COHORT"
