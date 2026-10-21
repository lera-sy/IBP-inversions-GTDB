# Module 1 — NCBI metadata retrieval, filtering, and species/assembly selection.
# Scientific scripts and selection thresholds are unchanged.

MIN_ASSEMBLIES = 50
ASSEMBLIES_PER_TAXON = 50
RANDOM_SEED = 42

FIELDS="accession,organism-name,organism-tax-id,assminfo-level,assminfo-sequencing-tech,assminfo-assembly-method,assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50,assmstats-scaffold-n50,assminfo-biosample-submission-date,assminfo-biosample-accession,assminfo-biosample-bioproject-accession"
FLAGS="--assembly-level complete --assembly-source refseq --mag exclude --as-json-lines"


rule fetch_metadata:
    output:
        "data/raw/all_prokaryotes_metadata.tsv"
    params:
        raw_dir = "data/raw",
        fields = FIELDS,
        flags = FLAGS
    shell:
        """
        mkdir -p {params.raw_dir}

        datasets summary genome taxon "Bacteria" {params.flags} > {params.raw_dir}/bacteria_metadata.jsonl
        dataformat tsv genome --inputfile {params.raw_dir}/bacteria_metadata.jsonl --fields "{params.fields}" > {output}

        datasets summary genome taxon "Archaea" {params.flags} > {params.raw_dir}/archaea_metadata.jsonl
        dataformat tsv genome --inputfile {params.raw_dir}/archaea_metadata.jsonl --fields "{params.fields}" | tail -n +2 >> {output}

        #Test if the output file is not empty
        n=$(( $(wc -l < {output}) - 1 ))
        if [ "$n" -le 0 ]; then
            echo "Error: No assemblies were retrieved. Please check the datasets command and flags."
            exit 1
        fi

        echo "Done retrieving metadata: $n assemblies --> {output}"
        """


rule filter_assemblies:
    input:
        metadata = "data/raw/all_prokaryotes_metadata.tsv"
    output:
        candidates = "data/processed/candidates.tsv",
        unknown = "data/processed/unknown.tsv",
        reference_guided = "data/processed/reference_guided.tsv"
    shell:
        """
        python -u src/ibp/module1/filter_assemblies.py -i {input.metadata} -c {output.candidates} -u {output.unknown} -r {output.reference_guided}
        
        # Test if the output files are not empty
        n=$(( $(wc -l < {output.candidates}) - 1 ))
        if [ "$n" -le 0 ]; then
            echo "Error: {output.candidates} is empty. No candidate assemblies were found."
            exit 1
        fi
        """

rule select_multiassembly_species:
    input:
        candidates = "data/processed/candidates.tsv"
    output:
        selected = "data/processed/selected_assemblies.tsv",
        summary = "data/processed/selection_summary.tsv"
    params:
        min_assemblies = MIN_ASSEMBLIES,
        assemblies_per_taxon = ASSEMBLIES_PER_TAXON,
        seed = RANDOM_SEED
    shell:
        """
        python -u src/ibp/module1/select_multiassembly_species.py -i {input.candidates} -o {output.selected} -s {output.summary} -m {params.min_assemblies} -n {params.assemblies_per_taxon} --seed {params.seed}
        
        # Test if the output files are not empty
        n=$(( $(wc -l < {output.selected}) - 1 ))
        if [ "$n" -le 0 ]; then
            echo "Error: {output.selected} is empty. No selected assemblies were found."
            exit 1
        fi
        """
