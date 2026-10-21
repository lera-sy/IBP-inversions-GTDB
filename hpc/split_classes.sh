#!/bin/bash
# The research question asked separately for spanning and local inversions. ibp-pipeline active:
#   bash hpc/split_classes.sh
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
if [ "$COHORT" = selected_50_v1 ]; then SEL=data/processed; else SEL=data/processed/runs/$COHORT; fi
PYTHONPATH=src python -m ibp.stats.split_classes --loci "results/loci/$COHORT" --tree "$SEL/species.tree" --null "results/model/$COHORT/null_slopes.tsv" --out "results/figures/$COHORT/split"
