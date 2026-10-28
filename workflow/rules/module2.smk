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
    # Resumable: progress is kept in genomes_work (not an output), only missing files
    # are retried, and the few genomes NCBI no longer serves are listed and skipped
    # (genome QC then drops them). Fails only if > 5% are missing.
    input:
        code = "src/ibp/module2/robust_download.py",
        package = DOWNLOAD_DIR + "/genomes.zip",
        accessions = f"{RUN_PROCESSED_DIR}/selected_accessions.txt"
    output:
        genomes = directory(DOWNLOAD_DIR + "/genomes"),
        missing = f"{RUN_PROCESSED_DIR}/download_missing_genomes.tsv"
    params:
        work = DOWNLOAD_DIR + "/genomes_work"
    threads: 2
    shell:
        """
        PYTHONPATH=src python -u -m ibp.module2.robust_download genomes \
            --package {input.package:q} --accessions {input.accessions:q} \
            --work {params.work:q} --out {output.genomes:q} \
            --missing {output.missing:q} --workers {threads} --rounds 4 --max-missing 0.05
        """

# Preparation rule for the genome_qc

rule fetch_sequence_reports:
    # Same idea: batches of 500 into sequence_reports_work, per accession for a batch
    # that fails; missing reports are listed and those assemblies fail QC.
    input:
        code = "src/ibp/module2/robust_download.py",
        accessions = f"{RUN_PROCESSED_DIR}/selected_accessions.txt"
    output:
        reports = directory(DOWNLOAD_DIR + "/sequence_reports"),
        missing = f"{RUN_PROCESSED_DIR}/download_missing_reports.tsv"
    params:
        work = DOWNLOAD_DIR + "/sequence_reports_work"
    shell:
        """
        PYTHONPATH=src python -u -m ibp.module2.robust_download reports \
            --accessions {input.accessions:q} --work {params.work:q} \
            --out {output.reports:q} --missing {output.missing:q} --batch 500 --max-missing 0.05
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

