# Module 5 — Inversion Detection

## Purpose

Detect candidate chromosomal inversions by comparing
prepared query chromosomes against one selected reference
chromosome per species.

Callers:
- Leen
- minimap2 + SyRI

## Input

One Module 4 species/cohort package:

- manifest.json
- genomes.tsv
- genomes/<assembly_id>.fna

The input follows the Module 4 → 5 Contract v1.1.

Module 5 must not modify these files.

## Processing

1. Validate the input package.
2. Identify the reference chromosome.
3. Identify all accepted query chromosomes.
4. Run Leen for each reference/query pair.
5. Run minimap2 + SyRI for each pair.
6. Validate and standardize caller outputs.

Each caller and query pair is an independent job.

## Output

For each species, cohort, caller and query:

- calls.tsv: detected inversion candidates
- status.json: execution and callability status
- run.log: execution diagnostics

Successful comparisons with no inversions must be
distinguished from failed or unsupported comparisons.

Output coordinates refer to the normalized reference
chromosome supplied by Module 4.

## Responsibilities

Module 5 DOES:
- Validate prepared genome inputs
- Execute Leen and SyRI
- Record inversion candidates
- Standardize output formats
- Record execution status and diagnostics

Module 5 DOES NOT:
- Download or select genomes
- Perform genome QC or normalization
- Select references
- Group inversions into evolutionary families
- Infer independent evolutionary events
- Perform phylogenetic reconstruction
- Perform cross-species statistical analysis

## Execution

Snakemake manages job dependencies and parallel execution.

Scientific caller implementations must also be executable
independently for testing and debugging.

The same module must work locally and on VSC.

## Initial tests

1. Read a mock three-genome species package.
2. Validate reference and query identification.
3. Generate two Leen and two SyRI comparison jobs.
4. Run existing callers on synthetic inversion examples.
5. Verify standardized outputs and failure handling.
6. Verify Snakemake can execute and resume the workflow.

## Standardized inversion results

Each successful caller comparison produces:

- calls.tsv: original caller-specific inversion records
- inversions.tsv: standardized inversion records
- status.json: execution and resolution status
- run.log: execution diagnostics

Additional native output files may be retained for debugging.

The canonical reference coordinate convention is zero-based,
half-open, relative to the Module 4 normalized chromosome.

Important standardized fields:

- caller: inversion detection method
- reference_assembly: selected reference accession
- query_assembly: compared query accession
- call_id: identifier of the detected inversion
- ref_start_0 / ref_end_0: reference interval coordinates
- ref_intervals_0: reference intervals, including circular cases
- qry_start_0 / qry_end_0: query coordinates when available
- reported_length_bp: inversion length reported by the caller
- interval_span_bp: length implied by the reference intervals
- coordinate_method: origin of the breakpoint estimates

Leen's window-based breakpoint estimates may differ from its
reported inversion length. Both measurements are preserved.

SyRI native coordinates are converted from one-based inclusive
to zero-based half-open coordinates.

An empty inversions.tsv is valid only for a successfully
resolved comparison with zero inversion calls.

Unresolved, unsupported and failed comparisons must not be
interpreted as genomes without inversions.

This schema currently covers candidate inversions from both
callers. It does not establish independently evolved events.

## Integrated execution (Modules 1–5)

The repository root `Snakefile` includes `workflow/rules/module5.smk`.
The original Module 3 normalization checkpoint determines which species have
accepted chromosomes, and Module 4 `prepare_species` is a checkpoint so that
Snakemake can discover the actual accepted comparison queries **after** the
package has been published. Module 5 reads but never modifies the package.

```bash
# On VSC, VSC_SCRATCH is normally already defined. For a local test,
# set it to a location with enough free disk space first.
conda activate ibp-pipeline
snakemake --snakefile Snakefile \
  --configfile config/runs/module5_integrated.yaml \
  --cores 8 --dry-run module5_complete

# Only when the planned workload is acceptable:
snakemake --snakefile Snakefile \
  --configfile config/runs/module5_integrated.yaml \
  --cores 8 module5_complete
```

Optionally add `--config species='["escherichia_coli_562"]'` to restrict the
species **analyzed after normalization**. This does **not** restrict the global
Module 1/2 metadata selection or downloads. The first uncached integrated run
can therefore involve thousands of NCBI assemblies and should be run with
appropriate VSC scratch storage/compute, not casually on a laptop.

Each cohort report is published at:

`$VSC_SCRATCH/IBP-main/module5/<species_id>/<cohort_id>/comparisons.tsv`

The comparison-specific `inversions.tsv`, `calls.tsv`, `status.json`, and
`run.log` are stored under the same directory in
`<caller>/<query_assembly>/`. The report distinguishes `success`,
`unresolved`, `unsupported`, and `failed`; downstream code must inspect status
rather than treating an empty inversion file as a biological zero.

The combined `ibp-pipeline` environment needs minimap2. SyRI is deliberately
invoked in the separate `ibp-syri` Conda environment configured through
`syri_env`. Ensure both environments are installed on the machine running jobs.

### Standalone Module 5 mock regression run

```bash
snakemake --snakefile workflow/rules/module5.smk \
  --configfile config/runs/module5_mock.yaml \
  --cores 4 module5_complete
```

This mode keeps the original three-genome fixture and does not trigger
Modules 1–4. It is useful for verifying the Module 5 code before a costly
integrated run.

Note: a successful DAG/dry run only confirms scheduling. A complete scientific
validation still requires actual genome data and examination of caller outputs.
