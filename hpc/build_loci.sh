#!/bin/bash
# Calls -> loci -> families for all species. Light (a few minutes); login node is fine.
# Run from the repo folder with ibp-pipeline active:  bash hpc/build_loci.sh
PYTHONPATH=src python -m ibp.loci.build --module5-root "$VSC_SCRATCH/IBP-main/module5" --prepared-root "$VSC_SCRATCH/IBP-main/prepared" --species-table data/processed/species.tsv --cohort selected_50_v1 --floor-bp 6000 --tolerance-bp 5000 --max-failed 0.20 --out results/loci/selected_50_v1
