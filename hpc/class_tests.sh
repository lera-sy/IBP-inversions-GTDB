#!/bin/bash
# Spanning vs local, formal tests: phylogenetic interaction test, empirical p-values against
# both nulls (same model and procedure), spanning share with clustered CIs, symmetry.
# Needs results/robustness/<cohort>/termini_gc_skew.tsv (run hpc/robustness.sh first).
# About 10-20 minutes; login node is fine. Run from the repo folder:
#   source $VSC_DATA/miniforge3/bin/activate && conda activate ibp-pipeline
#   bash hpc/class_tests.sh
# activates ibp-pipeline itself if it is not active
if [ "${CONDA_DEFAULT_ENV:-}" != ibp-pipeline ]; then source "$VSC_DATA/miniforge3/bin/activate"; conda activate ibp-pipeline; fi
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
if [ "$COHORT" = selected_50_v1 ]; then SEL=data/processed; else SEL=data/processed/runs/$COHORT; fi
PYTHONPATH=src python -u -m ibp.stats.class_tests --loci "results/loci/$COHORT" --tree "$SEL/species.tree" --taxonomy "$SEL/species.tsv" --termini "results/robustness/$COHORT/termini_gc_skew.tsv" --replicates 200 --out "results/class_tests/$COHORT"
