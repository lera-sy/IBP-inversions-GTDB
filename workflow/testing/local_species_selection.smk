# Integration-test-only selection; loaded only with run_assembly_limit.

rule select_run_assemblies:
    """Select the assemblies entering Modules 2-5, retaining Module 1 outputs."""
    input:
        selected = f"{SELECTION_DIR}/selected_assemblies.tsv",
        species_table = SPECIES_TABLE
    output:
        RUN_ASSEMBLIES
    params:
        species = config.get("species", []),
        max_per_species = int(config.get("run_assembly_limit", 0)),
        seed = RANDOM_SEED
    run:
        requested = params.species
        if not isinstance(requested, list) or not requested:
            raise ValueError("Specify a nonempty species list for run selection.")
        if params.max_per_species < 0:
            raise ValueError("run_assembly_limit must be >= 0 (0 means all).")

        names = pd.read_csv(input.species_table, sep="\t", dtype=str)
        if names["species_id"].duplicated().any():
            raise ValueError("Duplicate species IDs in species table.")
        lookup = names.set_index("species_id")["taxid"].to_dict()
        taxids = [lookup.get(str(species), str(species)) for species in requested]
        if len(set(taxids)) != len(taxids):
            raise ValueError("Requested species contains duplicates.")

        selected = pd.read_csv(input.selected, sep="\t", dtype=str)
        column = "Organism Taxonomic ID"
        if column not in selected or "Assembly Accession" not in selected:
            raise ValueError("Module 1 selected_assemblies.tsv lacks required columns.")

        subsets = []
        for taxid in taxids:
            subset = selected.loc[selected[column] == taxid]
            if subset.empty:
                raise ValueError(f"Species TaxID {taxid} is absent from Module 1 selection.")
            if 0 < params.max_per_species < len(subset):
                subset = subset.sample(n=params.max_per_species, random_state=params.seed)
            subsets.append(subset)

        result = pd.concat(subsets).sort_values([column, "Assembly Accession"])
        result.to_csv(output[0], sep="\t", index=False)
        print(f"Run selection: {len(result)} assemblies from {len(taxids)} species.")
