# Phase 3: null models.

NULL = config.get("null_model", {})


rule null_geometry:
    # Phase 3A: slope produced by the shorter-arc cap and detection floor alone.
    input:
        code = "src/ibp/nullmodel/geometry.py"
    output:
        "results/null_model/geometry_slopes.tsv"
    params:
        floor = config["definitions"]["detection_floor_bp"]["main"],
        min_events = config["definitions"]["min_loci_per_species"]["pgls_means"],
        replicates = NULL.get("replicates", 200),
        seed = NULL.get("seed", 42)
    shell:
        """
        PYTHONPATH=src python -u -m ibp.nullmodel.geometry --out {output:q} \
            --lengths-mb 1 10 --n-species 80 --n-events 30 \
            --replicates {params.replicates} --floor-bp {params.floor} \
            --min-events {params.min_events} --seed {params.seed}
        """
