#!/bin/bash
# Checks on the local-inversion slope. Run after hpc/run_model.sh, ibp-pipeline active:
#   bash hpc/check_local.sh                     # min-50 run
#   COHORT=selected_min20_v1 bash hpc/check_local.sh
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
if [ "$COHORT" = selected_50_v1 ]; then SEL=data/processed; else SEL=data/processed/runs/$COHORT; fi
PYTHONPATH=src python -m ibp.stats.check_local --loci "results/loci/$COHORT" --tree "$SEL/species.tree" --null "results/model/$COHORT/null_slopes.tsv" --out "results/model/$COHORT/local_checks"
