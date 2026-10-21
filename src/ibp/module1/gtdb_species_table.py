"""
Group complete RefSeq genomes into GTDB species clusters and prepare the species tree.


Steps:
1. Load the GTDB bacterial metadata and keep RefSeq, complete, non-MAG genomes.
2. Apply the CheckM2 completeness / contamination filter.
3. Keep one genome per BioSample within each GTDB species.
4. Keep species with at least --min-genomes genomes.
5. Prune the GTDB tree to those species (one tip per species, named after it).

Works from the command line or inside a notebook: `from gtdb_species_table import run`.

Note: genome_size in GTDB is the total assembly length (plasmids included). The
chromosome length used in the final model has to come from the sequences later.
"""

import argparse

import pandas as pd

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
THRESHOLDS = [5, 10, 20, 50, 100]

# Column names differ slightly between GTDB releases: first match wins.
COLUMN_OPTIONS = {
    "accession": ["accession"],
    "taxonomy": ["gtdb_taxonomy"],
    "representative": ["gtdb_genome_representative"],
    "completeness": ["checkm2_completeness", "checkm_completeness"],
    "contamination": ["checkm2_contamination", "checkm_contamination"],
    "genome_size": ["genome_size"],
    "assembly_level": ["ncbi_assembly_level"],
    "genome_category": ["ncbi_genome_category"],
    "biosample": ["ncbi_biosample"],
}
EXCLUDED_CATEGORIES = r"metagenome|environmental|single cell"


def resolve_columns(path):
    """Map our internal names to the column names present in this GTDB release."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    columns, missing = {}, []
    for key, options in COLUMN_OPTIONS.items():
        found = [c for c in options if c in header]
        if found:
            columns[key] = found[0]
        else:
            missing.append(f"{key} (looked for: {', '.join(options)})")
    if missing:
        raise ValueError("Columns not found in GTDB metadata:\n  " + "\n  ".join(missing))
    if columns["completeness"] != "checkm2_completeness":
        print(f"WARNING: no CheckM2 columns, using {columns['completeness']} instead.")
    return columns


def load_and_filter(path, min_completeness, max_contamination):
    """Return complete, non-MAG RefSeq genomes passing the quality filter."""
    cols = resolve_columns(path)
    df = pd.read_csv(path, sep="\t", usecols=list(cols.values()), dtype=str)
    df = df.rename(columns={v: k for k, v in cols.items()})
    print(f"Loaded {len(df):,} genomes from GTDB metadata.")

    df = df[df["accession"].str.startswith("RS_")]
    print(f"RefSeq only: {len(df):,}")

    df = df[df["assembly_level"].str.strip().str.lower() == "complete genome"]
    print(f"Complete genomes: {len(df):,}")

    is_excluded = df["genome_category"].fillna("").str.contains(EXCLUDED_CATEGORIES, case=False)
    df = df[~is_excluded]
    print(f"Not MAG / SAG / environmental: {len(df):,}")

    for col in ["completeness", "contamination", "genome_size"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    n_missing = int(df[["completeness", "contamination"]].isna().any(axis=1).sum())
    if n_missing:
        print(f"WARNING: {n_missing:,} genomes have no quality values and are dropped.")
    df = df[(df["completeness"] >= min_completeness) & (df["contamination"] <= max_contamination)]
    print(f"Completeness >= {min_completeness}, contamination <= {max_contamination}: {len(df):,}")

    ranks = df["taxonomy"].str.split(";", expand=True)
    if ranks.shape[1] != len(RANKS):
        raise ValueError("gtdb_taxonomy does not have the expected 7 ranks.")
    for i, rank in enumerate(RANKS):
        df[rank] = ranks[i].str.strip().str.replace(r"^[a-z]__", "", regex=True)
    df = df[df["species"] != ""].drop(columns="taxonomy")

    # GTDB prefixes accessions with RS_ / GB_, NCBI does not
    df["assembly_accession"] = df["accession"].str.replace(r"^(RS|GB)_", "", regex=True)
    return df


def deduplicate_biosamples(df):
    """Keep one genome per BioSample within each species (best quality first)."""
    df = df.copy()
    key = df["biosample"].fillna("").str.strip()
    df["biosample_key"] = key.where(~key.isin(["", "none"]), df["accession"])
    df["quality"] = df["completeness"] - 5 * df["contamination"]
    df = df.sort_values(["species", "quality", "accession"], ascending=[True, False, True])
    out = df.drop_duplicates(["species", "biosample_key"]).drop(columns="biosample_key")
    print(f"One genome per BioSample per species: {len(out):,}")
    return out


def summarise_species(df):
    """One row per species: counts, genome size, tree tip and alignment reference."""
    rows = []
    for species, group in df.groupby("species", sort=True):
        representative = group["representative"].iloc[0]
        rep_passes = representative in set(group["accession"])
        # group is sorted best quality first, so iloc[0] is the fallback reference
        reference = representative if rep_passes else group["accession"].iloc[0]
        rows.append({
            "species": species,
            "tree_label": species.replace(" ", "_"),
            "n_genomes": len(group),
            "median_genome_size": int(group["genome_size"].median()),
            "min_genome_size": int(group["genome_size"].min()),
            "max_genome_size": int(group["genome_size"].max()),
            "gtdb_representative": representative,
            "representative_passes_filters": rep_passes,
            "alignment_reference": reference,
            **{rank: group[rank].iloc[0] for rank in RANKS[:-1]},
        })
    return pd.DataFrame(rows)


def threshold_table(species):
    """How many species (and what genome-size range) survive each threshold."""
    rows = []
    for t in THRESHOLDS:
        kept = species[species["n_genomes"] >= t]
        sizes = kept["median_genome_size"] / 1e6
        rows.append({
            "min_genomes": t,
            "species": len(kept),
            "phyla": kept["phylum"].nunique(),
            "smallest_Mb": round(sizes.min(), 2) if len(kept) else None,
            "largest_Mb": round(sizes.max(), 2) if len(kept) else None,
        })
    return pd.DataFrame(rows)


def prune_tree(tree_path, species, out_path):
    """Prune the GTDB tree to one tip per species and rename tips to species names."""
    import dendropy  # imported here so the selection step does not need it

    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    print(f"GTDB tree loaded: {len(tree.taxon_namespace):,} tips.")
    labels = dict(zip(species["gtdb_representative"], species["tree_label"]))
    in_tree = {t.label for t in tree.taxon_namespace}
    missing = sorted(set(labels) - in_tree)
    if missing:
        print(f"WARNING: {len(missing)} representatives are not in the tree, e.g. {missing[:3]}")
    keep = set(labels) & in_tree
    if len(keep) < 3:
        raise RuntimeError("Fewer than 3 species found in the tree: check the GTDB release.")
    tree.retain_taxa_with_labels(keep)
    tree.purge_taxon_namespace()
    for taxon in tree.taxon_namespace:
        taxon.label = labels[taxon.label]
    tree.write(path=out_path, schema="newick", suppress_internal_node_labels=True,
               suppress_rooting=True, unquoted_underscores=True)
    print(f"Pruned tree: {len(keep)} species --> {out_path}")
    return missing


def run(metadata, tree, out_genomes, out_species, out_tree,
        min_completeness=95.0, max_contamination=5.0, min_genomes=50):
    genomes = load_and_filter(metadata, min_completeness, max_contamination)
    genomes = deduplicate_biosamples(genomes)
    species = summarise_species(genomes)

    table = threshold_table(species)
    print("\nSpecies kept at each minimum number of genomes:")
    print(table.to_string(index=False))

    species = species[species["n_genomes"] >= min_genomes].reset_index(drop=True)
    if species.empty:
        raise RuntimeError(f"No species with >= {min_genomes} genomes.")
    genomes = genomes[genomes["species"].isin(species["species"])].drop(columns="quality")
    print(f"\nKept {len(species)} species / {len(genomes):,} genomes with >= {min_genomes} genomes.")

    missing = prune_tree(tree, species, out_tree)
    species["in_tree"] = ~species["gtdb_representative"].isin(missing)

    genomes.to_csv(out_genomes, sep="\t", index=False)
    species.to_csv(out_species, sep="\t", index=False)
    print(f"Written: {out_genomes}, {out_species}")
    return genomes, species, table


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metadata", required=True, help="GTDB bac120 metadata (.tsv or .tsv.gz)")
    parser.add_argument("-t", "--tree", required=True, help="GTDB bac120 tree (newick)")
    parser.add_argument("-g", "--out-genomes", required=True)
    parser.add_argument("-s", "--out-species", required=True)
    parser.add_argument("-o", "--out-tree", required=True)
    parser.add_argument("--min-completeness", type=float, default=95.0)
    parser.add_argument("--max-contamination", type=float, default=5.0)
    parser.add_argument("--min-genomes", type=int, default=50)
    args = parser.parse_args()
    run(args.metadata, args.tree, args.out_genomes, args.out_species, args.out_tree,
        args.min_completeness, args.max_contamination, args.min_genomes)
