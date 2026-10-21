# Inversion-only multi-species production runner (Modules 1–8)

This wrapper executes **unchanged** Modules 1–6 once per species using isolated
cohort IDs, then calls the existing Modules 7 and 8 on validated catalogues.
It does **not** perform phylogeny, rooting or independent event inference.
Caller-specific candidate inversion families are not biological independent events.

## Quick start (WSL, `ibp-pipeline` activated)

Module 1 must have produced `data/processed/selected_assemblies.tsv` and
`config/species.tsv` must map each selected TaxID to a unique species ID.
The production YAML deliberately omits `species`; omission means all selected
species. `run_assembly_limit: 0` retains the existing 50-assembly selection.
Use a smaller `run_assembly_limit` in a **new cohort** for a bounded test.

```bash
cd ~/IBP-module5
export VSC_SCRATCH="$HOME/ibp-local-scratch"
python scripts/run_study.py --configfile config/runs/module8_all_species_production.yaml \
  --cores 4 --species-workers 1 --plan
# For a real production run, remove --plan:
python scripts/run_study.py --configfile config/runs/module8_all_species_production.yaml \
  --cores 4 --species-workers 1
```

Optional `--species-workers 2 --cores 8` schedules two independent species with
four CPU slots each. With only four available CPU slots, use one species worker
until resource use and memory are measured. `--dry-run` evaluates each isolated
species' Snakemake DAG without running scientific tools. Re-executing the same
command is **resumable** via Snakemake `--rerun-incomplete` and distinct cohort
paths. Never launch two copies of the *same* study ID simultaneously.

## Fault tolerance

Each species has its own Snakemake run and workspace. A failed download, QC,
normalization, reference preparation, caller rule, or corrupted catalogue is
recorded in `species_status.tsv`; remaining species continue. Completed species
are aggregated; failed species **are not interpreted as zero inversions**.
Module 5's successful/failed/unresolved/unsupported comparison statuses remain
available; a successful catalogue with unresolved comparisons is marked
`partial`, not `failed`, and unknown genome states remain unknown.

If every species fails, the runner exits nonzero and produces no aggregate.
If only some fail, it exits zero with explicit failures in the status TSV;
inspect and retry before treating the aggregate as final. The summary is a
**subset** until all required species are complete. A fatal aggregate error
(e.g. a cohort identity mismatch) stops publication; it is not hidden.

## Files

`$VSC_SCRATCH/IBP-studies/<study-cohort>/` contains:

- `species_status.tsv`: status (`completed`, `partial`, `failed`), original
  source-cohort ID, comparison counts, log path and failure reason.
- `logs/<species>.log`: complete Snakemake stdout/stderr per species.
- `configs/<species>.yaml`: reproducible per-species config.
- `module7/<study-cohort>/*.tsv`: verified subset aggregate, including source
  identities normalized to the study cohort. The source IDs remain in the
  status report and original Module 6 catalogues.
- `module8/<study-cohort>/*.tsv` and `figures/*.png`: descriptive results.

Each species' raw outputs remain under
`$VSC_SCRATCH/IBP-main/runs/<study-cohort>__<species-id>/`.
Module 1's selection stays shared/read-only. The original pilot cohort
and all existing Module 1–8 Snakemake rules are unchanged.

## Why normalization is slow

The existing Module 3 normalization **already runs chromosomes in parallel**.
It submits one Dnaapler process per chromosome to a Python thread pool, with
`--workers` equal to Snakemake's `threads` (normally `--cores`). Each Dnaapler
process uses `--threads 1` to avoid oversubscribing CPU cores. The subprocess
may perform substantial MMseqs2 work, so 16 chromosomes can still take many
minutes. The change in this patch only adds per-genome wall-time logs and
out-of-order progress reporting, preserving normalized FASTA content and table
order. On a system with more physical cores, increase `--cores`, or use more
concurrent species within the same total CPU budget. Confirm that disk and
memory use remain reasonable; allocating more CPU slots is not guaranteed to
speed up I/O-limited normalization.

Normalization may need further profiling/caching for a full 50-genome × many-
species study. This runner does not claim to remove that work.

## Caller-specific agreement analysis

Module 8 also emits `caller_agreement.tsv`, `caller_agreement_calls.tsv`,
`caller_agreement_families.tsv`, `figures/caller_agreement.png`, and
`figures/agreement_family_lengths.png`. These use Module 7's existing
caller-separated original calls and provisional matches; they need no rerun
of Leen, SyRI or normalization when only analysis code changes.
Consult `docs/modules/module_7_8.md` for definitions and missingness rules.
