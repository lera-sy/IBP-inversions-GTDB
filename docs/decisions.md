# Decisions log

Every change to a definition in `config/config.yaml` after Phase 1 is logged here: date, what changed, why.
Newest first.

| Date | Decision | Reason |
|---|---|---|
| 2026-10-21 | Pipeline is fully GTDB-driven: GTDB r232 decides species, genome list, quality and tree; NCBI is used only to download the FASTA files and sequence reports by accession. Team's NCBI TaxID selection (module 1) is replaced in Phase 2. | Independent pipeline; avoids split and pooled species. GTDB does not host sequences for all genomes. |
| 2026-10-21 | Main genome size = whole genome length (GTDB `genome_size`); analysed chromosome length is the sensitivity check. | The research question is about genome size. The size cap still uses chromosome length inside the null model. |
| 2026-10-21 | Minimum loci per species: 1 for the mixed model, 3 for PGLS on species means; dropped species are listed. | The mixed model weighs species by how many loci they have; a mean from 1-2 loci is too noisy for PGLS. |
| 2026-10-21 | 50 genomes per species (min and max); may lower the minimum to 20 after the Phase 2 counts. | Same depth as the team run; 20 only if too few species reach 50. |
| 2026-10-21 | No 90/10 train/test split. Tuning set = pilot species + simulations; confirmation set = all other species, run only after the config is frozen. Robustness: leave one phylum out. | A random split does not fit a regression test; the real risk is tuning definitions on the data being tested. |
| 2026-10-21 | No fixed deadline: long-read validation and the wider genome size range stay in scope. | A few days available. |
| 2026-10-21 | Repo started from the team's `feature/multispecies-integration` (Modules 1–8) plus analysis scripts from `leen-final-run`. | One clean modular pipeline to build on; caller algorithm is identical in both branches. |
| 2026-10-21 | Phase 1 definitions entered as **proposed** in `config/config.yaml`. | To be confirmed together before any run. |
