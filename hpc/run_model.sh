#!/bin/bash
# Phase 3 S0 + Phase 7 model. Run from the repo folder with ibp-pipeline active:
#   bash hpc/run_model.sh
# About 1-3 minutes; login node is fine.
set -euo pipefail
COHORT=selected_50_v1
LOCI=results/loci/$COHORT
OUT=results/model/$COHORT
export PYTHONPATH=src

python -m ibp.stats.model --loci "$LOCI" --tree data/processed/species.tree --out "$OUT" --write-null-input
echo "== S0: geometry null with the real pooled sizes, chromosome lengths and family counts"
python -m ibp.nullmodel.geometry --out "$OUT/null_slopes.tsv" --species-table "$OUT/null_species.tsv" --empirical-sizes "$OUT/null_sizes.tsv" --replicates 200 --floor-bp 6000 --min-events 3 --seed 42 | grep -E "Lengths|empirical"
echo "== Models"
python -m ibp.stats.model --loci "$LOCI" --tree data/processed/species.tree --null "$OUT/null_slopes.tsv" --out "$OUT"
