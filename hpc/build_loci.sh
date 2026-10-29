#!/bin/bash
# Calls -> loci -> families for all species. Light (a few minutes); login node is fine.
# Run from the repo folder with ibp-pipeline active:
#   bash hpc/build_loci.sh                      # min-50 run
#   COHORT=selected_min20_v1 bash hpc/build_loci.sh
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
if [ "$COHORT" = selected_50_v1 ]; then ROOT=$VSC_SCRATCH/IBP-main; SEL=data/processed
else ROOT=$VSC_SCRATCH/IBP-main/runs/$COHORT; SEL=data/processed/runs/$COHORT; fi
PYTHONPATH=src python -m ibp.loci.build --module5-root "$ROOT/module5" --prepared-root "$ROOT/prepared" --species-table "$SEL/species.tsv" --cohort "$COHORT" --floor-bp 6000 --tolerance-bp 5000 --max-failed 0.20 --workers "${SLURM_CPUS_PER_TASK:-1}" --out "results/loci/$COHORT"
