#!/bin/bash
# Distribution/normality checks and the ori vs ter split (about 5 min; login node is fine).
#   bash hpc/extra_checks.sh
#   COHORT=selected_min20_v1 bash hpc/extra_checks.sh
if [ "${CONDA_DEFAULT_ENV:-}" != ibp-pipeline ]; then source "$VSC_DATA/miniforge3/bin/activate"; conda activate ibp-pipeline; fi
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
if [ "$COHORT" = selected_50_v1 ]; then SEL=data/processed; else SEL=data/processed/runs/$COHORT; fi
export PYTHONPATH=src OMP_NUM_THREADS=1
python -m ibp.stats.distribution_checks --loci "results/loci/$COHORT" --tree "$SEL/species.tree" --boot 1000 --out "results/extra/$COHORT"
python -m ibp.stats.ori_ter --loci "results/loci/$COHORT" --tree "$SEL/species.tree" \
  --termini "results/robustness/$COHORT/termini_gc_skew.tsv" --replicates 200 --out "results/extra/$COHORT" | grep -v "null "
