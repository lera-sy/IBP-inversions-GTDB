#!/bin/bash
# Finish the min-50 run: loci with the raw-alignment tiers -> S0 + models ->
# local checks -> comparison figures with the team run. A few minutes; login
# node is fine. Run from the repo folder with ibp-pipeline active:
#   bash hpc/finish_50.sh
set -euo pipefail
export COHORT=selected_50_v1
echo "########## 1. Loci and tiers (Leen + SyRI + raw alignment)"
bash hpc/build_loci.sh
echo; echo "########## 2. S0 and models"
bash hpc/run_model.sh
echo; echo "########## 3. Local-inversion checks"
bash hpc/check_local.sh | head -n 25
echo; echo "########## 4. Comparison with the team run"
PYTHONPATH=src python -m ibp.stats.compare_team --loci results/loci/$COHORT --model results/model/$COHORT --out results/compare/$COHORT
echo; echo "########## 5. Figures for this run"
PYTHONPATH=src python -m ibp.stats.plot_run --loci results/loci/$COHORT --model results/model/$COHORT --out results/figures/$COHORT
