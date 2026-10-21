# Decisions log

Every change to a definition in `config/config.yaml` after Phase 1 is logged here: date, what changed, why.
Newest first.

| Date | Decision | Reason |
|---|---|---|
| 2026-10-21 | Preparation run: 3,987 chromosomes normalized (3,934 circular, 53 linear), 80 species packages. Linear chromosomes are kept in the packages for now and excluded at the analysis step (shorter-arc rule does not apply to them). | config `exclude_linear_chromosomes: true`; the exclusion is applied when loci are built. |
| 2026-10-21 | After the whitelist fix: 80 species at >= 50 (4,000 genomes, 7 phyla, 1.05-7.18 Mb); 201 at >= 20 (11 phyla, 0.81-9.14 Mb). First full run uses min 50; a separate min-20 run follows in its own folders. | The 50 run was already past download and QC; 20 widens the size range and doubles the species, so it is run next and compared. |
| 2026-10-21 | First selection run: 76 GTDB species at >= 50 genomes (3,800 genomes, 7 phyla, 1.05-7.18 Mb median size); 178 species at >= 20. Whitelist extended with Autocycler, Trycycler, Dragonflye, Hybracter, SMRT Analysis/SMRT Link Microbial Assembly. Left out: SOAPdenovo, SPAdes, 'various', unknown pipelines. | 12,076 genomes were 'unknown method'; the top ones are long-read de novo tools. Short-read assemblers and unclear labels stay out. |
| 2026-10-21 | Module 1 rebuilt: GTDB r232 filters (RefSeq, complete, not MAG, CheckM2) joined with the NCBI long-read / de novo filter; one genome per BioSample per GTDB species; >= 50 genomes, 50 sampled (seed 42). | GTDB decides species; the long-read / de novo filter is kept because reference-guided or short-read assemblies can hide or fake inversions. |
| 2026-10-21 | Species key = number in the GTDB representative accession (RS_GCF_000005845.2 -> 5845), used where downstream code says "taxid". Species folder = GTDB name + key, e.g. `escherichia_coli_5845`. | Downstream modules need an integer key; this one is unique and ties each species to its tree tip. NCBI TaxID kept only as a reference column. |
| 2026-10-21 | GTDB and NCBI joined on accession without version; the current NCBI version is downloaded. | GTDB r232 can list an older version of the same assembly. |
| 2026-10-21 | Bacteria only (include_archaea: false). | No ar53 files on the VSC. |
| 2026-10-21 | Pipeline is fully GTDB-driven: GTDB r232 decides species, genome list, quality and tree; NCBI is used only to download the FASTA files and sequence reports by accession. Team's NCBI TaxID selection (module 1) is replaced in Phase 2. | Independent pipeline; avoids split and pooled species. GTDB does not host sequences for all genomes. |
| 2026-10-21 | Main genome size = whole genome length (GTDB `genome_size`); analysed chromosome length is the sensitivity check. | The research question is about genome size. The size cap still uses chromosome length inside the null model. |
| 2026-10-21 | Minimum loci per species: 1 for the mixed model, 3 for PGLS on species means; dropped species are listed. | The mixed model weighs species by how many loci they have; a mean from 1-2 loci is too noisy for PGLS. |
| 2026-10-21 | 50 genomes per species (min and max); may lower the minimum to 20 after the Phase 2 counts. | Same depth as the team run; 20 only if too few species reach 50. |
| 2026-10-21 | No 90/10 train/test split. Tuning set = pilot species + simulations; confirmation set = all other species, run only after the config is frozen. Robustness: leave one phylum out. | A random split does not fit a regression test; the real risk is tuning definitions on the data being tested. |
| 2026-10-21 | No fixed deadline: long-read validation and the wider genome size range stay in scope. | A few days available. |
| 2026-10-21 | Repo started from the team's `feature/multispecies-integration` (Modules 1–8) plus analysis scripts from `leen-final-run`. | One clean modular pipeline to build on; caller algorithm is identical in both branches. |
| 2026-10-21 | Phase 1 definitions entered as **proposed** in `config/config.yaml`. | To be confirmed together before any run. |
