#!/bin/bash
# Phase 4 pilot report. Run from the repo folder with ibp-pipeline active:
#   bash hpc/pilot_report.sh
PYTHONPATH=src python -m ibp.pilot.report --module5-root "$VSC_SCRATCH/IBP-main/module5" --species escherichia_coli_3697165 vibrio_cholerae_621645 staphylococcus_aureus_1027105 bordetella_pertussis_306945 --out results/pilot
