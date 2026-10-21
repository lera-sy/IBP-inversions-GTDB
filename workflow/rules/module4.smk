# Module 4 — reference selection and validated packages.
rule select_reference:
    input:
        code = "src/ibp/module4/select_reference.py",
        fasta_code = "src/ibp/module3/normalize_chromosomes.py",
        metadata = f"{RUN_PROCESSED_DIR}/normalized_assemblies.tsv",
        report = f"{RUN_PROCESSED_DIR}/normalization.tsv",
        fastas = reference_fastas
    output:
        reference = REFERENCE_DIR + "/reference.txt",
        ranking = REFERENCE_DIR + "/reference_ranking.tsv",
        distances = REFERENCE_DIR + "/mash_distances.tsv",
        sketch = REFERENCE_DIR + "/chromosomes.msh"
    params:
        outdir = REFERENCE_DIR
    threads: 1
    shell:
        """
        PYTHONPATH=src python -u -m ibp.module4.select_reference --metadata {input.metadata:q} \
            --normalization {input.report:q} --taxid {wildcards.species:q} \
            --fastas {input.fastas:q} --output-dir {params.outdir:q}
        """


rule select_references:
    input:
        reference_targets


rule prepare_packages:
    input:
        prepared_targets


checkpoint prepare_species:
    input:
        code = "src/ibp/module4/prepare_species.py",
        fasta_code = "src/ibp/module3/normalize_chromosomes.py",
        selected = RUN_ASSEMBLIES,
        accepted = f"{RUN_PROCESSED_DIR}/normalized_assemblies.tsv",
        qc = f"{RUN_PROCESSED_DIR}/genome_qc.tsv",
        normalization = f"{RUN_PROCESSED_DIR}/normalization.tsv",
        names = SPECIES_TABLE,
        reference = lambda wildcards: REFERENCE_DIR.format(species=wildcards.species.rsplit("_", 1)[-1]) + "/reference.txt",
        fastas = reference_fastas
    output:
        package = protected(directory(PREPARED_DIR))
    threads: 1
    shell:
        """
        PYTHONPATH=src python -u -m ibp.module4.prepare_species --selected {input.selected:q} \
            --accepted {input.accepted:q} --qc {input.qc:q} \
            --normalization {input.normalization:q} --reference {input.reference:q} \
            --species-table {input.names:q} --species {wildcards.species:q} \
            --cohort {COHORT:q} --fastas {input.fastas:q} --output-dir {output.package:q}
        """

