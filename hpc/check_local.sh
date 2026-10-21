#!/bin/bash
# Checks on the local-inversion slope. Run after hpc/run_model.sh, ibp-pipeline active:
#   bash hpc/check_local.sh
set -euo pipefail
COHORT=selected_50_v1
PYTHONPATH=src python -m ibp.stats.check_local --loci results/loci/$COHORT --tree data/processed/species.tree --null results/model/$COHORT/null_slopes.tsv --out results/model/$COHORT/local_checks
