import shlex

# Module 3 — chromosome topology and normalization.
rule fetch_chromosome_topology:
    input:
        code = "src/ibp/module3/fetch_chromosome_topology.py",
        qc = f"{RUN_PROCESSED_DIR}/genome_qc.tsv"
    output:
        f"{RUN_PROCESSED_DIR}/chromosome_topology.tsv"
    shell:
        """
        python -u {input.code:q} --qc {input.qc:q} --output {output:q}
        """


checkpoint normalize_chromosomes:
    input:
        code = "src/ibp/module3/normalize_chromosomes.py",
        metadata = f"{RUN_PROCESSED_DIR}/qc_passed_assemblies.tsv",
        topology = f"{RUN_PROCESSED_DIR}/chromosome_topology.tsv",
        chromosomes = DOWNLOAD_DIR + "/chromosomes"
    output:
        report = f"{RUN_PROCESSED_DIR}/normalization.tsv",
        passed = f"{RUN_PROCESSED_DIR}/normalized_assemblies.tsv",
        normalized = directory(DOWNLOAD_DIR + "/normalized_chromosomes")
    params:
        logs = DOWNLOAD_DIR + "/normalization_logs",
        # Resumable: every result is also written to normalization_cache (outside
        # Snakemake's output, so it survives a crash or time limit) and reused on a
        # rerun after the same sequence checks. `normalization_reuse` in the config
        # adds earlier runs' FASTA:LOG folders (e.g. the min-50 run).
        cache_out = DOWNLOAD_DIR + "/normalization_cache",
        cache = lambda wildcards: " ".join(
            "--cache " + shlex.quote(os.path.expandvars(entry))
            for entry in [f"{DOWNLOAD_DIR}/normalization_cache:{DOWNLOAD_DIR}/normalization_logs",
                          *config.get("normalization_reuse", [])])
    threads: workflow.cores
    shell:
        """
        python -u {input.code:q} --metadata {input.metadata:q} --topology {input.topology:q} \
            --chromosomes {input.chromosomes:q} --normalized {output.normalized:q} \
            --logs {params.logs:q} --report {output.report:q} --passed {output.passed:q} \
            --workers {threads} --cache-out {params.cache_out:q} {params.cache}
        """

