# Where the code came from

Both sources were uploaded as branch zips on 2026-10-21 (no git history in the zips, so no commit hash).

| In this repo | Source branch | Source path | Changes |
|---|---|---|---|
| `src/ibp/module1`–`module8`, `src/ibp/production` | `feature/multispecies-integration` | `src/ibp/...` | Lint fix only (variable rename in `module6/catalog.py`) |
| `tests/` | `feature/multispecies-integration` | `tests/` | None |
| `workflow/rules/*.smk`, `Snakefile`, `scripts/run_study.py` | `feature/multispecies-integration` | same | `Snakefile` now loads `config/config.yaml` by default |
| `envs/pipeline.yaml` | `feature/multispecies-integration` | `environment.yaml` | Added Slurm executor plugin |
| `envs/syri.yaml` | `leen-final-run` | `envs/syri.yaml` | None |
| `src/ibp/analysis_legacy/final_inversion_model.py` | `leen-final-run` | `scripts/steps/analysis/` | None yet; to be adapted in Phases 6–7 |
| `src/ibp/analysis_legacy/breakpoint_repeats.py` | `leen-final-run` | `scripts/steps/analysis/` | None yet; Phase 6 |
| `src/ibp/analysis_legacy/plot_position_classes.py`, `plot_species_boxplots.py` | `leen-final-run` | `scripts/steps/analysis/` | None yet; Phase 7 |
| `config/team_run/species_ncbi.tsv` | `feature/multispecies-integration` | `config/species.tsv` | Kept only to compare species sets with the team run |
| `docs/team_run/` | `feature/multispecies-integration` | `docs/` | None (reference) |

The Leen caller in `src/ibp/module5/leen_core.py` was checked function by function against
`leen-final-run/scripts/ibp/inversions/leen.py`: all 14 algorithm functions are identical
(including the refined clean / in_place / shared_breakpoint classification).

Not included: the NUCmer+SyRI caller (it lives on `abi/integrated-pipeline`, not in either zip),
the Parsnp/Gubbins phylogeny code (plan step: skipped), EDA figures and one-off VSC fix scripts.
