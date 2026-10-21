#!/bin/bash
# Spanning vs local: significance per class and the class difference. ibp-pipeline active:
#   bash hpc/plot_classes.sh
set -euo pipefail
COHORT=${COHORT:-selected_50_v1}
PYTHONPATH=src python -m ibp.stats.plot_classes --model "results/model/$COHORT" --out "results/figures/$COHORT"
