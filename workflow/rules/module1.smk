# Module 1 — selection. GTDB is the only source of species identity.
# NCBI metadata is fetched only for sequencing technology and assembly method
# (long-read / de novo filter); GTDB decides species, quality and tree tip.

SELECTION = config["selection"]
GTDB_METADATA = os.path.expandvars(str(SELECTION.get("gtdb_metadata") or ""))
GTDB_TREE = os.path.expandvars(str(SELECTION.get("gtdb_tree") or ""))
if not GTDB_METADATA or not GTDB_TREE:
    raise ValueError("Set selection.gtdb_metadata and selection.gtdb_tree in config/config.yaml")
NCBI_TAXA = ["Bacteria", "Archaea"] if SELECTION.get("include_archaea") else ["Bacteria"]

FIELDS="accession,organism-name,organism-tax-id,assminfo-level,assminfo-sequencing-tech,assminfo-assembly-method,assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50,assmstats-scaffold-n50,assminfo-biosample-submission-date,assminfo-biosample-accession,assminfo-biosample-bioproject-accession"
FLAGS="--assembly-level complete --assembly-source refseq --mag exclude --as-json-lines"


rule fetch_metadata:
    output:
        "data/raw/all_prokaryotes_metadata.tsv"
    params:
        raw_dir = "data/raw",
        fields = FIELDS,
        flags = FLAGS,
        taxa = " ".join(NCBI_TAXA)
    shell:
        """
        mkdir -p {params.raw_dir}

        echo -n > {output}
        for taxon in {params.taxa}; do
            datasets summary genome taxon "$taxon" {params.flags} > {params.raw_dir}/${{taxon}}_metadata.jsonl
            if [ -s {output} ]; then
                dataformat tsv genome --inputfile {params.raw_dir}/${{taxon}}_metadata.jsonl --fields "{params.fields}" | tail -n +2 >> {output}
            else
                dataformat tsv genome --inputfile {params.raw_dir}/${{taxon}}_metadata.jsonl --fields "{params.fields}" > {output}
            fi
        done

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
        metadata = "data/raw/all_prokaryotes_metadata.tsv",
        code = "src/ibp/module1/filter_assemblies.py"
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

rule gtdb_select:
    # GTDB decides species, quality and tree tip; NCBI candidates only add the
    # long-read / de novo filter. See src/ibp/module1/gtdb_select.py.
    input:
        candidates = "data/processed/candidates.tsv",
        gtdb = GTDB_METADATA,
        code = "src/ibp/module1/gtdb_select.py"
    output:
        selected = f"{SELECTION_DIR}/selected_assemblies.tsv",
        species = f"{SELECTION_DIR}/species.tsv",
        summary = f"{SELECTION_DIR}/selection_summary.tsv",
        steps = f"{SELECTION_DIR}/selection_steps.tsv",
        thresholds = f"{SELECTION_DIR}/selection_thresholds.tsv"
    params:
        completeness = SELECTION["checkm2_completeness_min"],
        contamination = SELECTION["checkm2_contamination_max"],
        min_genomes = SELECTION["min_genomes_per_species"],
        max_genomes = SELECTION["max_genomes_per_species"],
        seed = SELECTION["seed"]
    resources:
        mem_mb = 16000
    shell:
        """
        PYTHONPATH=src python -u -m ibp.module1.gtdb_select \
            --gtdb-metadata {input.gtdb:q} --candidates {input.candidates:q} \
            --selected {output.selected:q} --species {output.species:q} \
            --summary {output.summary:q} --steps {output.steps:q} --thresholds {output.thresholds:q} \
            --min-completeness {params.completeness} --max-contamination {params.contamination} \
            --min-genomes {params.min_genomes} --max-genomes {params.max_genomes} --seed {params.seed}
        """


rule gtdb_tree:
    input:
        species = f"{SELECTION_DIR}/species.tsv",
        tree = GTDB_TREE,
        code = "src/ibp/module1/gtdb_tree.py"
    output:
        tree = f"{SELECTION_DIR}/species.tree",
        missing = f"{SELECTION_DIR}/species_not_in_tree.tsv"
    resources:
        mem_mb = 16000
    shell:
        """
        PYTHONPATH=src python -u -m ibp.module1.gtdb_tree --tree {input.tree:q} \
            --species {input.species:q} --out {output.tree:q} --missing {output.missing:q}
        """
