"""Prune the GTDB tree to one tip per selected species, tips named by species_id."""

import argparse

import pandas as pd


def prune(tree_path, species_table, out_path, missing_path):
    import dendropy  # only needed here; loading 190k tips takes a few minutes

    species = pd.read_csv(species_table, sep="\t", dtype=str)
    labels = dict(zip(species["gtdb_representative"], species["species_id"]))
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    in_tree = {t.label for t in tree.taxon_namespace}
    print(f"GTDB tree loaded: {len(in_tree):,} tips.")
    missing = sorted(set(labels) - in_tree)
    pd.DataFrame({"gtdb_representative": missing,
                  "species_id": [labels[m] for m in missing]}).to_csv(missing_path, sep="\t", index=False)
    keep = set(labels) & in_tree
    if len(keep) < 3:
        raise RuntimeError("Fewer than 3 selected species are in the tree: check the GTDB release.")
    tree.retain_taxa_with_labels(keep)
    tree.purge_taxon_namespace()
    for taxon in tree.taxon_namespace:
        taxon.label = labels[taxon.label]
    tree.write(path=out_path, schema="newick", suppress_internal_node_labels=True,
               suppress_rooting=True, unquoted_underscores=True)
    print(f"Pruned tree: {len(keep)} species ({len(missing)} not in tree) -> {out_path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--missing", required=True)
    a = p.parse_args()
    prune(a.tree, a.species, a.out, a.missing)
