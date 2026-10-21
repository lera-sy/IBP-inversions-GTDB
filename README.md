# IBP-inversions-GTDB

Does inversion size scale with genome size in prokaryotes?
A clean rerun of the IBP 2026 analysis, with every definition fixed before running, GTDB species from the start, a null model, and locus-level confidence tiers.
It runs in parallel to the team run (`Integrated-Bioinformatics-Project-2026`) so the two results can be compared.

## Status

| Phase | What | Status |
|---|---|---|
| 0 | Setup: repo, envs, config, dry run | done |
| 1 | Lock definitions (`config/config.yaml`) | answered 21 Oct, see `docs/decisions.md` |
| 2 | GTDB species selection | done: 80 species, 3,987 chromosomes |
| 3 | Null model and detectability | 3A geometry done; S0 needs Phase 5 sizes |
| 4 | Pilot: 3 contrasting species + *B. pertussis* | running |
| 5 | Full run | to do |
| 6 | Validation and locus tiers | to do |
| 7 | Statistics and comparison with the team run | to do |

## Layout

```
config/config.yaml        all definitions and settings (one source of truth)
config/runs/              run-specific overrides (--configfile)
docs/decisions.md         every change to a definition, with date and reason
docs/provenance.md        which code came from which team branch
envs/                     conda envs: pipeline, syri, stats, rrna
profiles/slurm/           Snakemake profile for VSC wICE
src/ibp/module1-8         pipeline modules (from the team's multispecies branch)
src/ibp/analysis_legacy/  final model, breakpoint repeats, plots (from leen-final-run)
workflow/rules/           Snakemake rules, one file per module
tests/                    unit tests
```

## Setup on the VSC (once)

In a VSC terminal (login node is fine for this part):

```bash
cd $VSC_DATA
git clone https://github.com/lera-sy/ibp-inversions-gtdb.git
cd ibp-inversions-gtdb

source $VSC_DATA/miniforge3/bin/activate
mamba env create -f envs/pipeline.yaml   # ibp-pipeline (if it already exists: mamba env update -n ibp-pipeline -f envs/pipeline.yaml)
mamba env create -f envs/syri.yaml       # ibp-syri (skip if you already have it)
mamba env create -f envs/stats.yaml      # ibp-stats (skip if you already have it)
mamba env create -f envs/rrna.yaml       # ibp-rrna (skip if you already have it)
```

Check that everything works:

```bash
conda activate ibp-pipeline
PYTHONPATH=src python -m unittest discover -s tests
ruff check .
snakemake --snakefile workflow/rules/module5.smk \
  --configfile config/runs/module5_mock.yaml --cores 1 -n module5_complete
```

The last command is a dry run on a 3-genome mock package; it should end with `This was a dry-run`.

## Running

On Slurm, set your credit account once per session and use the profile:

```bash
export SBATCH_ACCOUNT=<your lp_ account>
snakemake --profile profiles/slurm -n      # always dry run first
```

### Phase 2: GTDB selection

Needs internet (NCBI metadata) and about 16 GB memory; 2 cores are enough.

```bash
conda activate ibp-pipeline
cd $VSC_DATA/IBP-inversions-GTDB
snakemake -n --cores 2 data/processed/selected_assemblies.tsv data/processed/species.tree
snakemake --cores 2 data/processed/selected_assemblies.tsv data/processed/species.tree
```

Outputs in `data/processed/`: `selected_assemblies.tsv`, `species.tsv`, `species.tree`,
`selection_steps.tsv` (genomes left after each filter), `selection_thresholds.tsv`
(species kept at 5/10/20/30/50 genomes), `selection_summary.tsv` (per species),
`species_not_in_tree.tsv`.

### Phase 4: pilot (4 species)

```bash
sbatch --account=lp_edu_large_omics hpc/pilot_callers.slurm
# when it has finished (seconds, login node is fine):
PYTHONPATH=src python -m ibp.pilot.report --module5-root $VSC_SCRATCH/IBP-main/module5 \
    --species escherichia_coli_3697165 vibrio_cholerae_621645 \
              staphylococcus_aureus_1027105 bordetella_pertussis_306945 \
    --out results/pilot
```

Exact commands for later phases are added here as each phase is built.
