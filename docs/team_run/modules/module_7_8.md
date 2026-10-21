# Modules 7 and 8 — Tree-free aggregation and descriptive analysis

Module 7 aggregates the **caller-specific** families, genome lengths, unknown tip states, and provisional cross-caller associations produced by Module 6. Missing input files, species IDs, incomplete tip states, wrong schemas and duplicate species are errors, **not zero inversion counts**. Both caller rows retain their own denominators and the aggregate association count is not double-counted as two independent events.

Module 8 produces five descriptive tables and five PNG figures, adapting the earlier inversion-size/genome-size and caller-comparison analyses to the current standardized Module 6 catalogue. The **candidate family length** is the median *caller-reported length* for a family. Genome size plots use the median `genome_length_bp` of accepted assemblies; per-species summaries also retain the selected-chromosome length. Cross-species comparisons are descriptive and do not adjust for phylogenetic relatedness, sampling or caller sensitivity. At least three species, with variation and adequate QC, are needed for meaningful cross-species relationships. No inferential relationship or p value is reported here.

No output is named an independently evolved mutation: breakpoint families can merge independent changes or split one complex rearrangement, and matching Leen/SyRI breakpoints are only **provisional associations**.

## Run

After applying this patch and running its one-time installer:

```bash
VSC_SCRATCH="$HOME/ibp-local-scratch" snakemake --snakefile Snakefile \
  --configfile config/runs/module5_local_ecoli.yaml \
  --config cohort=local_ecoli_8_refactor_v1 --cores 4 --dry-run module8_complete

VSC_SCRATCH="$HOME/ibp-local-scratch" snakemake --snakefile Snakefile \
  --configfile config/runs/module5_local_ecoli.yaml \
  --config cohort=local_ecoli_8_refactor_v1 --cores 4 module8_complete
```

Outputs under `$VSC_SCRATCH/IBP-main/module7/<cohort>/` and `$VSC_SCRATCH/IBP-main/module8/<cohort>/`.

The `module8_complete` target is the one-command entry point for all **successfully prepared and selected species in the config**. As currently wired, upstream Snakemake species failures still halt the overall DAG: unattended per-species retry/isolation is *not* implemented by this patch. A multi-species run and failure recovery still require separate validation. This patch deliberately doesn't modify Modules 1–6 or require a phylogenetic tree.

## Caller agreement (Module 8 extension)

Module 8 additionally reads Module 7's `call_family_map.tsv` and
`cross_caller_associations.tsv`, including the standardized per-query call IDs,
plus genome-level `leen_callable` / `syri_callable` flags. Neither Leen nor SyRI
families are merged or destroyed. The following are **additional** outputs:

- `caller_agreement.tsv`: one row per species with joint-callable **query** genomes,
  total Leen/SyRI calls, uniquely matched **pairs** (counted once), caller-only
  calls, separate ambiguous Leen and SyRI *call* counts, unassessable calls,
  and counts of caller-specific families containing at least one unique match.
- `caller_agreement_calls.tsv`: one row per original call, retaining species,
  caller, query accession, original call ID, family ID, length, assessment
  category and possible matched call IDs. Categories are exclusive:
  `unique_candidate`, `ambiguous_candidate`, `caller_only`, `unassessable`.
- `caller_agreement_families.tsv`: one row per caller-specific candidate family,
  with counts of member calls in each category. The `support_class` is a
  convenience label with precedence *unique* > *ambiguous* > *caller-only* >
  *not evaluable*. A `has_unique_candidate` family may ALSO have unmatched or
  ambiguous member calls, which remain explicit in its count columns.
- `figures/caller_agreement.png`: stacked cross-species comparison for the 25
  species with the most assessable calls. The TSV contains **all species**.
  Unique associations count as pairs; other categories count as individual
  Leen/SyRI calls. Unassessable calls are excluded from the plot and included
  in the table. This is **not** a single additive biological event count.
- `figures/agreement_family_lengths.png`: per-caller length distributions of
  deduplicated families stratified by the family support classes. The figure
  uses the existing family median caller-reported length (not independent events).

"Caller-only" is allowed only when both callers finished a comparison for
that particular query genome; if the other caller failed or was unsupported,
the call is **unassessable**, never a negative finding. An ambiguous one-to-many
set is reported as its constituent ambiguous calls, not as one event.
`unique_candidate` is a reciprocal one-to-one breakpoint-proximity association
**within a query genome**, NOT a proven biological consensus inversion or an
inferred independent evolutionary inversion. These observations must not be
used as a sensitivity/specificity estimate in the absence of truth data.

The new TSVs and figures are additional to the original Module 8 outputs;
Modules 1–7 and the species-level runner's behavior are unchanged.
