# Opt-in cohort-isolated Modules 2–8 execution

Legacy behaviour is unchanged unless `isolate_run: true` is supplied.

## Why

The original Modules 2–3 write shared `data/processed/*.tsv` and `IBP-main/*` files.
Changing the list of selected species in a new run previously reused those paths,
causing Snakemake to reschedule upstream downloads and normalization and possibly
invalidate a protected Module 4 package. This is unsafe for simultaneous runs.

## Isolated layout

For a config with `cohort: local_two_species_8_v1` and `isolate_run: true`:

- Module 1 metadata and the canonical full selection stay at `data/raw/` and
  `data/processed/selected_assemblies.tsv`; the cohort shares these read-only inputs.
- Selected accessions, local assembly subsampling, QC, topology and normalization
  TSVs go into `data/processed/runs/local_two_species_8_v1/`.
- Downloaded genomes, sequence reports, chromosomes, normalized FASTAs and
  Modules 4–8 go under `$VSC_SCRATCH/IBP-main/runs/local_two_species_8_v1/`.
- The original validated E. coli results under `$VSC_SCRATCH/IBP-main/` and
  the existing `data/processed/*.tsv` stay intact.

The mode is **cohort-isolation**, not per-species failure isolation. A species
whose upstream QC or package preparation errors out can still stop the cohort,
and Module 7 currently requires all expected Module 6 summaries. This needs
separate handling before executing the intended ~95-species study.

## Local two-species test

In `config/runs/module8_local_two_species.yaml` add `isolate_run: true`;
keep `run_assembly_limit: 8`, the two explicit species IDs and the unique
`cohort` value. Preview with:

```sh
VSC_SCRATCH="$HOME/ibp-local-scratch" snakemake --snakefile Snakefile \
  --configfile config/runs/module8_local_two_species.yaml \
  --cores 4 --dry-run module8_complete
```

Expect fresh downstream work for 16 genomes under the isolated workspace,
while Module 1's already-generated full selection is reused. Avoid concurrent
runs that modify Module 1's canonical metadata selection.
