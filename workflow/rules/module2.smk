# Module 2 — downloads and genome quality control.
rule write_accessions:
    input:
        RUN_ASSEMBLIES
    output:
        f"{RUN_PROCESSED_DIR}/selected_accessions.txt"
    run:
        selected = pd.read_csv(input[0], sep="\t", dtype=str)
        accessions = selected["Assembly Accession"].dropna().drop_duplicates()
        if accessions.empty:
            raise ValueError("No selected assembly accessions to download.")
        accessions.to_csv(output[0], index=False, header=False)


rule download_genome_package:
    input:
        f"{RUN_PROCESSED_DIR}/selected_accessions.txt"
    output:
        DOWNLOAD_DIR + "/genomes.zip"
    shell:
        """
        # One bulk request for the small package of metadata and download links.
        for delay in 0 30 60; do
            sleep "$delay"
            if datasets download genome accession --inputfile {input:q} \
                --include genome --dehydrated --filename {output:q}; then
                break
            fi
            if [ "$delay" -eq 60 ]; then exit 1; fi
        done
        """


rule download_genomes:
    input:
        package = DOWNLOAD_DIR + "/genomes.zip",
        accessions = f"{RUN_PROCESSED_DIR}/selected_accessions.txt"
    output:
        directory(DOWNLOAD_DIR + "/genomes")
    threads: 2
    shell:
        """
        python -m zipfile -e {input.package:q} {output:q}

        # Limit concurrent downloads, retrying after 30 and then 60 seconds.
        for delay in 0 30 60; do
            sleep "$delay"
            if datasets rehydrate --directory {output:q} --max-workers {threads}; then
                break
            fi
            if [ "$delay" -eq 60 ]; then exit 1; fi
        done

        # Require a nonempty FASTA for every selected accession.
        while read -r accession; do
            test -s {output:q}/ncbi_dataset/data/"$accession"/*_genomic.fna
        done < {input.accessions:q}
        """

# Preparation rule for the genome_qc

rule fetch_sequence_reports:
    input:
        f"{RUN_PROCESSED_DIR}/selected_accessions.txt"
    output:
        package = DOWNLOAD_DIR + "/sequence_reports.zip",
        reports = directory(DOWNLOAD_DIR + "/sequence_reports")
    shell:
        """
        # Download only sequence reports, using the same selected accessions.
        for delay in 0 30 60; do
            sleep "$delay"
            if datasets download genome accession --inputfile {input:q} \
                --include seq-report --filename {output.package:q}; then
                break
            fi
            if [ "$delay" -eq 60 ]; then exit 1; fi
        done

        python -m zipfile -e {output.package:q} {output.reports:q}

        # Require a nonempty report for every selected assembly.
        while read -r accession; do
            test -s {output.reports:q}/ncbi_dataset/data/"$accession"/sequence_report.jsonl
        done < {input:q}
        """


rule genome_qc:
    input:
        code = "src/ibp/module2/genome_qc.py",
        metadata = RUN_ASSEMBLIES,
        genomes = DOWNLOAD_DIR + "/genomes",
        reports = DOWNLOAD_DIR + "/sequence_reports"
    output:
        qc = f"{RUN_PROCESSED_DIR}/genome_qc.tsv",
        passed = f"{RUN_PROCESSED_DIR}/qc_passed_assemblies.tsv",
        chromosomes = directory(DOWNLOAD_DIR + "/chromosomes")
    shell:
        """
        python -u {input.code:q} --metadata {input.metadata:q} \
            --genomes {input.genomes:q} --reports {input.reports:q} \
            --chromosomes {output.chromosomes:q} --qc {output.qc:q} --passed {output.passed:q}
        """

