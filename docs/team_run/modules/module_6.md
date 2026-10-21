# Module 6 — first integrated research slice (tree-free)

Based on legacy `prepare_leen_states.py`, `prepare_syri_states.py`,
`build_inversion_family_dataset.py` and exploratory caller comparisons.
The legacy `infer_family_parsimony.py` is retained for later phylogenetic inference,
not silently executed by this module.

## Installing

1. Extract this ZIP at the root of `refactor/modules1-5-layout`.
2. Run `python ibp_enable_module6.py` once; it adds a single `include` line
   after the Module 5 include in the root `Snakefile`. The installer is temporary.
3. Run `PYTHONPATH=src python -m unittest discover -s tests/module6 -v`
   and `ruff check src/ibp/module6 tests/module6`.
4. Dry-run on the cached local cohort:

   ```bash
   VSC_SCRATCH="$HOME/ibp-local-scratch" snakemake \
     --snakefile Snakefile \
     --configfile config/runs/module5_local_ecoli.yaml \
     --cores 4 --dry-run module6_complete
   ```

5. If the dry run is correct, remove `--dry-run` to execute Module 6.
   Module 1–5 outputs are reused; only Module 6 should run unless provenance
   changes invalidate upstream jobs.

Default output root: `<VSC_SCRATCH>/IBP-main/module6/<species>/<cohort>/`.
The `module6_complete` target uses the same checkpoint-aware species list as
Module 5. This first iteration leaves the root `all` target unchanged.

## Outputs

- `candidate_families.tsv`: complete-link grouped calls within **each caller**;
  these are reference-relative candidates, **not evolutionary events**.
- `call_family_map.tsv`: every call has exactly one caller-specific family.
- `tip_states.tsv`: `1` = detected in a successful query, `0` = **not detected**
  relative to the selected reference in a successful query (or reference itself), `?` =
  unsuccessful/unsupported/unresolved comparison. These are neither ancestral states nor proof of true biological absence.
- `cross_caller_associations.tsv`: same-query breakpoint-proximity candidates,
  marked `unique_candidate` only for reciprocal degree-one associations;
  multiple possible matches remain `ambiguous_candidate`. These associations
  are **not** validated consensus calls.
- `genome_sizes.tsv`: reference and query chromosome/genome lengths from the
  accepted Module 4 package.
- `summary.tsv`: per-caller counts, failure denominators and association counts.

## Scientific / implementation limits

- The legacy default breakpoint tolerance is 5,000 bp, configurable as
  `module6.breakpoint_tolerance_bp` (0 or larger) in YAML. Grouping uses a
  deterministic greedy **complete-link** procedure; different thresholds may
  change the family catalogue. This is exploratory, not a validated classifier.
- Geometry represented by a single reference interval, including a circular
  wrap, can be grouped. Two pieces meeting exactly at the reference origin
  are safely combined. Other disjoint segments are retained as **unmerged
  singleton complex-geometry candidates**; they are not treated as absence.
- A match requires both left and right breakpoints within tolerance, circular
  modulo chromosome length only for circular references. The caller-specific
  definitions of inversion length are preserved, not conflated.
- No tree inference, independent origin counting, ancestor reconstruction,
  rearrangement-direction inference, or cross-species statistics are asserted.
- Input `comparisons.tsv` must be the exact Module 5 schema with one status
  row per query and caller. Module 4 FASTAs are checksum-validated by the
  existing package loader.
