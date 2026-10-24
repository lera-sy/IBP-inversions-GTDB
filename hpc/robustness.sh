#!/bin/bash
# Robustness of the spanning vs local result: terminus position (shift, margin, GC skew),
# random-breakpoint null, one species per genus/family, whole distribution.
# A few minutes; login node is fine. ibp-pipeline active:
#   bash hpc/robustness.sh
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
if [ "$COHORT" = selected_50_v1 ]; then SEL=data/processed; ROOT=$VSC_SCRATCH/IBP-main
else SEL=data/processed/runs/$COHORT; ROOT=$VSC_SCRATCH/IBP-main/runs/$COHORT; fi
PYTHONPATH=src python -m ibp.stats.robustness --loci "results/loci/$COHORT" --tree "$SEL/species.tree" --taxonomy "$SEL/species.tsv" --prepared-root "$ROOT/prepared" --cohort "$COHORT" --replicates 100 --picks 200 --out "results/robustness/$COHORT"
