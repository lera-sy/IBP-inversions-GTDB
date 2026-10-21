# Running the inversion-only IBP pipeline on VSC wICE

This runbook documents the **99-species × up to 50 chromosomes/species** production workflow used in October 2026 on VSC wICE, with a separate `ibp-pipeline` environment and `ibp-syri` environment. It is a **tree-free, reference-relative inversion analysis**: Leen and SyRI report candidate inversions, Module 6 builds caller-specific candidate families and provisional cross-caller associations, and Modules 7–8 produce descriptive cross-species analyses. **These are not phylogenetically proven independent evolutionary inversion events.**

The instructions assume a VSC wICE account with authorized credits, `batch` partition access, the repository branch `feature/multispecies-integration`, and the version-controlled `hpc/` helper scripts listed below. **Do not run compute-heavy steps on a login node.** Submit them via `sbatch`. Paths and accounting rules are cluster/account-specific.

## 1. Reproducible environment and checkout

On a VSC **login node** (`tier2-p-login-*`), from a fresh shell:

```bash
cd "$VSC_DATA/projects/IBP-module5"
git switch feature/multispecies-integration
export PATH="$VSC_DATA/bin:$PATH"
export MAMBA_ROOT_PREFIX="$VSC_DATA/micromamba"
export VSC_SCRATCH="$VSC_SCRATCH_LUSTRE1"
export PYTHONPATH="$PWD/src"
micromamba env list
```

This particular deployment uses the existing named environments:

- Main environment: `$VSC_DATA/conda_envs/ibp-pipeline` (Python 3.13, Snakemake and Modules 1–8 dependencies; see `environment.yaml`).
- SyRI environment: `$VSC_DATA/conda_envs/ibp-syri` (separate compatible SyRI dependencies).

The environments must already have been created and validated. `environment.yaml` describes the main environment; do **not** assume it installs the separate SyRI environment. Record versions for reproducibility:

```bash
micromamba run -n ibp-pipeline python --version
"$VSC_DATA/conda_envs/ibp-syri/bin/syri" --version
snakemake --version   # only if the active shell exposes the main environment
```

If an executable does not support `--version`, consult its help instead. The `syri` executable is ordinarily `$VSC_DATA/conda_envs/ibp-syri/bin/syri` on this deployment, not under `$MAMBA_ROOT_PREFIX/envs/`.

**SyRI direct-launch fix:** `src/ibp/module5/syri.py` uses the existing executable directly when `VSC_DATA` is defined and the requested named environment has an executable at `$VSC_DATA/conda_envs/<env>/bin/syri`, adding that environment's `bin` to the subprocess `PATH`. It otherwise retains the conda/micromamba fallback. Previously, `micromamba run` was executed **for each comparison**, causing a shared `~/.cache/mamba/proc` lock bottleneck under high concurrency. One production alignment was independently reanalysed with direct SyRI: **2 parent INV records in both runs, byte-identical raw SyRI output**, verifier job `62366443` (2026-10-19). This is a one-comparison regression check, not a guarantee across every possible genome.

## 2. Inputs and taxonomic mapping (Module 1)

Module 1 discovers and selects qualifying NCBI RefSeq assemblies. Before a production run, ensure these files exist and reflect the intended **same frozen selection**:

- `data/processed/selected_assemblies.tsv`
- `data/processed/selection_summary.tsv`
- `config/species.tsv` (TaxID, reviewed canonical species name, stable species ID)
- `config/runs/module8_all_species_production.yaml` (study settings)

This October 2026 cohort had **99 species and 4,950 assemblies** (50/species). This is a property of that selection, **not** an unconditional pipeline guarantee. The taxonomic mapping includes **TaxID 28447 = *Clavibacter michiganensis*** at species rank; do not substitute the first assembly's subspecies string, which may refer to a different TaxID.

Validate the frozen selection and mapping before preparation:

```bash
micromamba run -n ibp-pipeline python hpc/vsc_check_species_mapping.py
micromamba run -n ibp-pipeline env PYTHONPATH="$PWD/src" python scripts/run_study.py --configfile config/runs/module8_all_species_production.yaml --cores 8 --species-workers 1 --plan
```

The mapping checker should show no unmapped TaxIDs. Confirm the study plan selects the intended species. Metadata acquisition/Module 1 itself is documented in the project README and workflow rules; this runbook **does not** automatically refresh Module 1 metadata. Re-running that selection can change the study cohort and invalidate frozen manifests.

## 3. Prepare a new, isolated 99-species SLURM study

**Only do this for a new study.** Do not rerun preparation against the currently active `all_species_inversions_v1` jobroot. The helper refuses to overwrite an existing `configs/` directory to protect active work.

```bash
micromamba run -n ibp-pipeline env PYTHONPATH="$PWD/src" python hpc/vsc_prepare_full_study.py --configfile config/runs/module8_all_species_production.yaml --expected-species 99 --cores-per-species 8 --memory-gb 24 --hours 4 --concurrent 6
JOBROOT="$VSC_SCRATCH/IBP-jobs/all_species_inversions_v1"
head -n 6 "$JOBROOT/species_manifest.tsv"
sam-balance
sam-quote sbatch "$JOBROOT/run_species.slurm"
```

The helper creates independent Snakemake workspaces under `$JOBROOT/workspaces/`, a frozen selection snapshot under `production_snapshots/all_species_inversions_v1/`, per-species YAMLs, the array script, and a separate `finalize.slurm`. It **does not submit jobs**. It records the checksum of `selected_assemblies.tsv` and rejects changes to that input while the array runs. Workspaces link to the repository's `src`, `config`, `data`, etc.; **edits to shared source during an active array can affect future processes**, so coordinate any code hotfix and record its exact timing and validation.

Adjust `--account`, `--partition`, `--cluster`, `--hours`, `--cores-per-species`, `--memory-gb`, and `--concurrent` for the actual authorized allocation. The current script defaults to account `lp_edu_evol_quant_genetics`, cluster `wice`, partition `batch`. Be considerate of this **shared course allocation**.

## 4. Submit and monitor the species array (Modules 2–6)

For an entirely new study after the preflight checks:

```bash
sbatch "$JOBROOT/run_species.slurm"
```

An array index corresponds to exactly one species in `species_manifest.tsv`. `--array=0-98%6` means up to six simultaneous species jobs, **each** requesting 8 CPUs, 24 GB memory and a maximum of four hours. It is not one 99-species job using only eight CPUs. For a pilot, override on the command line: `sbatch --array=0-3%2 "$JOBROOT/run_species.slurm"`; once complete, submit the rest with `sbatch --array=4-98%6 "$JOBROOT/run_species.slurm"` (do not run the full range again). The October 2026 array IDs were pilot `62364152` and remaining tasks `62364379`.

For an existing live study, **do not submit overlapping task indices**. Check first:

```bash
squeue -M wice -j 62364379
sacct -M wice -j 62364379 -X --format=JobID,State,Elapsed,ExitCode
sacct -M wice -j 62364379 -X -n --format=State | sort | uniq -c
sam-balance
sam-list-usagerecords --account=lp_edu_evol_quant_genetics --start=2026-10-18 --end=2026-10-20
du -sh "$VSC_SCRATCH/IBP-main/runs"
df -h "$VSC_SCRATCH"
```

Use the appropriate live job ID and date range, not the historical ones. An empty `squeue` does **not** prove success. `sacct` records SLURM completion, which is also distinct from scientific completeness. Review each species' `$JOBROOT/logs/<species_id>.log`, the `module5/.../comparisons.tsv` and Module 6 catalogues. A species can finish all Snakemake steps but have `unsupported`, `unresolved`, or failed caller comparisons; those are **not zeros** and will be flagged by the finalizer.

Concurrency can be adjusted on an existing array when site policy permits it; for example:

```bash
scontrol -M wice update JobId=62364379 ArrayTaskThrottle=8
squeue -M wice -j 62364379
```

Some completed array tasks can trigger an `already finished` error even when the throttle change took effect. Inspect the pending-task line for `%8`. More concurrency increases instantaneous credit burn and I/O pressure; it does **not** accelerate an individual eight-core species task. Avoid changing concurrency blindly if jobs are already waiting on shared resources.

**SyRI bottleneck symptoms/diagnostics:** comparisons previously logged `syri_seconds` near 1,700 seconds while `TotalCPU / (Elapsed × AllocCPUS)` was only a few percent. `micromamba run -n ibp-syri ...` visibly blocked on `~/.cache/mamba/proc`; direct SyRI startup took about 0.12 seconds. Examine `$.../syri/<query>/run.log` for `SyRI launcher: direct ...` after the patch. Don't remove cache locks or kill other users' processes. The array's outer `micromamba run -n ibp-pipeline` is still executed once per species; the patch specifically removes the **per-comparison** SyRI environment launch.

## 5. Recovery and resumption

Wait for running tasks to finish before assessing failures. Determine the failed **indices** with `sacct` and map them using `$JOBROOT/species_manifest.tsv` (first column). Diagnose each species log and whether the cause is a time limit (`TIMEOUT`), memory exhaustion (`OUT_OF_MEMORY`), missing inputs or a caller failure. Do not mistake `COMPLETED` for complete caller coverage.

After fixing the underlying issue and ensuring no duplicate task is running, resubmit **only specific failed indices**, e.g.:

```bash
sbatch --array=17,42%2 "$JOBROOT/run_species.slurm"
```

`--rerun-incomplete` lets Snakemake resume relevant missing outputs in that species' isolated workspace. Do not delete completed results or overwrite `selected_assemblies.tsv` to make a failure disappear. Keep `production_snapshots/`, run logs and Git commit identifiers as analysis provenance. Account for partial caller coverage in the final report even if retries cannot resolve it.

## 6. Final aggregation and validation (Modules 7–8)

**The species array runs only through `module6_complete`.** Modules 7 and 8 are generated by a separate finalization job **after the array is finished** and any necessary retries are done:

```bash
JOBROOT="$VSC_SCRATCH/IBP-jobs/all_species_inversions_v1"
sbatch "$JOBROOT/finalize.slurm"
```

The finalizer checks each Module 6 catalogue and per-genome caller status, writes `$VSC_SCRATCH/IBP-studies/all_species_inversions_v1/species_status.tsv`, and aggregates available validated catalogues into:

```text
$VSC_SCRATCH/IBP-studies/all_species_inversions_v1/
  species_status.tsv
  module7/all_species_inversions_v1/
  module8/all_species_inversions_v1/
    figures/
```

Inspect the finalizer job's SLURM status, `$JOBROOT/finalize-<job_id>.out`, the species-status table and all expected Module 8 figures. **Important:** the finalizer exits with status 2 if *any* species is failed or partial, including any with `unsupported` comparisons. It can still produce **partial descriptive results** from validated catalogues. Never describe those results as a complete, fully callable 99-species study until the statuses and denominators have been reviewed. The first four-species pilot already contained two `unsupported=1` comparisons, so a nonzero finalizer exit is possible even if all 99 array tasks show `COMPLETED`.

## 7. Saving the workflow in GitHub

Commit **source code, reviewed species IDs, stable helper scripts and this runbook**, but not data downloads, generated snapshots, HPC outputs, run logs or local source backups. Before staging, inspect `git status --short` and `git diff`. Then explicitly `git add` only desired files (do **not** run `git add .` indiscriminately). Exclude `hpc/backups/`, `production_snapshots/`, and any generated `syri-verify-*.out` logs. The source file `src/ibp/module5/syri.py` was updated **on VSC**, not automatically in GitHub.

A validated production run should record the Git commit (`git rev-parse HEAD`), immutable Module 1 selection SHA-256, per-species manifests, program/environment versions, the submission script, Slurm job IDs and the final `species_status.tsv`. Avoid accidentally committing private credentials, account secrets, or large datasets.

## 8. Historical validation notes (2026-10-18/19)

- E. coli 50-assembly smoke test: `62363995`, 98/98 successful Module 5 comparisons, 101 caller-specific candidate families and 6 provisional cross-caller associations.
- Initial four-species pilot: `62364152`, 4/4 SLURM tasks completed, with some `unsupported` comparisons recorded rather than hidden.
- Main 95-species continuation: `62364379`, initially `%6`, increased successfully to `%8` concurrency.
- Direct SyRI verifier: `62366443`, `COMPLETED` in 9 seconds, identical raw SyRI output on one existing stored alignment.
- On 2026-10-19, `src/ibp/module5/syri.py` was patched to avoid per-comparison micromamba cache locking. Results produced before and after this change should retain precise version/provenance information and be audited at finalization.
