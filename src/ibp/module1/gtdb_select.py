"""
Select species and genomes with GTDB as the only source of species identity.

GTDB decides: species, quality (CheckM2), which genomes count, the tree tip.
NCBI is only used for what GTDB does not have: sequencing technology and
assembly method (long-read / de novo filter, `filter_assemblies.py`) and,
later, the FASTA files themselves.

Steps
1. GTDB metadata: RefSeq, complete genome, not MAG/SAG, CheckM2 filter
   (`gtdb_species_table.load_and_filter`).
2. Keep only genomes that also passed the NCBI long-read / de novo filter
   (join on assembly accession).
3. One genome per BioSample within each GTDB species (best quality first).
4. Keep species with at least --min-genomes genomes; sample --max-genomes
   per species with a fixed seed.
5. Give each GTDB species a stable integer key: the number in its GTDB
   representative accession (RS_GCF_000005845.2 -> 5845). The downstream
   modules call this column "Organism Taxonomic ID" / "taxid"; here it is a
   GTDB species key, NOT an NCBI TaxID. The NCBI TaxID is kept as
   "NCBI Organism Taxonomic ID" for reference only.

Outputs
- selected_assemblies.tsv  input for modules 2-4 (one row per selected genome)
- species.tsv              taxid (GTDB key), scientific_name (GTDB species), species_id, taxonomy
- selection_summary.tsv    per species: genomes at each step
- selection_steps.tsv      genomes left after each filter (whole dataset)
- thresholds.tsv           species kept at 5/10/20/50 genomes, to decide 50 vs 20
"""

import argparse
import re

import pandas as pd

from ibp.module1.gtdb_species_table import RANKS, load_and_filter

THRESHOLDS = [5, 10, 20, 30, 50]


def species_key(representative):
    """RS_GCF_000005845.2 -> 5845 (unique per representative genome)."""
    match = re.fullmatch(r"(?:RS_|GB_)?GC[AF]_(\d{9})\.\d+", str(representative))
    if not match:
        raise ValueError(f"Unexpected GTDB representative accession: {representative}")
    return int(match.group(1))


def slug(name):
    text = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    if not text:
        raise ValueError(f"Cannot make a folder name from species {name!r}")
    return text


def join_long_read(gtdb, candidates):
    """Keep GTDB genomes that passed the NCBI long-read / de novo filter."""
    ncbi = candidates.rename(columns={"Organism Taxonomic ID": "NCBI Organism Taxonomic ID",
                                      "Organism Name": "NCBI Organism Name"})
    # Join without the version suffix: GTDB r232 may list GCF_x.1 while NCBI now has GCF_x.2.
    # The current NCBI version is the one downloaded later.
    ncbi["_base"] = ncbi["Assembly Accession"].str.replace(r"\.\d+$", "", regex=True)
    ncbi = ncbi.sort_values("Assembly Accession").drop_duplicates("_base", keep="last")
    gtdb = gtdb.assign(_base=gtdb["assembly_accession"].str.replace(r"\.\d+$", "", regex=True))
    merged = gtdb.merge(ncbi, on="_base", how="inner", validate="one_to_one")
    merged["version_changed"] = merged["assembly_accession"] != merged["Assembly Accession"]
    merged["assembly_accession"] = merged["Assembly Accession"]
    return merged.drop(columns="_base")


def deduplicate(df):
    df = df.copy()
    key = df["biosample"].fillna("").str.strip()
    ncbi_key = df.get("Assembly BioSample Accession", pd.Series("", index=df.index)).fillna("").str.strip()
    key = key.where(~key.isin(["", "none"]), ncbi_key)
    df["biosample_key"] = key.where(key != "", df["assembly_accession"])
    df["quality"] = df["completeness"] - 5 * df["contamination"]
    df = df.sort_values(["species", "quality", "assembly_accession"], ascending=[True, False, True])
    return df.drop_duplicates(["species", "biosample_key"]).drop(columns="biosample_key")


def threshold_table(counts, sizes, phyla):
    rows = []
    for t in THRESHOLDS:
        kept = counts[counts >= t].index
        mb = sizes.loc[kept] / 1e6
        rows.append({"min_genomes": t, "species": len(kept), "phyla": phyla.loc[kept].nunique(),
                     "smallest_median_Mb": round(mb.min(), 2) if len(kept) else None,
                     "largest_median_Mb": round(mb.max(), 2) if len(kept) else None})
    return pd.DataFrame(rows)


def select(gtdb, candidates, min_genomes, max_genomes, seed):
    steps = [("gtdb_refseq_complete_checkm2", len(gtdb))]
    joined = join_long_read(gtdb, candidates)
    steps.append(("and_ncbi_long_read_de_novo", len(joined)))
    unique = deduplicate(joined)
    steps.append(("one_per_biosample", len(unique)))

    counts = unique.groupby("species").size()
    sizes = unique.groupby("species")["genome_size"].median()
    phyla = unique.groupby("species")["phylum"].first()
    thresholds = threshold_table(counts, sizes, phyla)

    eligible = counts[counts >= min_genomes].index
    if len(eligible) == 0:
        raise RuntimeError(f"No GTDB species with >= {min_genomes} genomes after filtering.")
    picked = [g if len(g) <= max_genomes else g.sample(n=max_genomes, random_state=seed)
              for _, g in unique[unique["species"].isin(eligible)].groupby("species", sort=True)]
    selected = pd.concat(picked, ignore_index=True)
    steps.append((f"species_with_at_least_{min_genomes}", int(counts.loc[eligible].sum())))
    steps.append((f"selected_max_{max_genomes}_per_species", len(selected)))

    reps = selected.groupby("species")["representative"].first()
    keys = reps.map(species_key)
    if keys.duplicated().any():
        raise RuntimeError("Two GTDB species got the same key; check representatives.")
    selected["Organism Taxonomic ID"] = selected["species"].map(keys).astype(str)
    selected["Organism Name"] = selected["species"]

    species = (selected.groupby("species", sort=True)
               .agg(n_selected=("assembly_accession", "size"),
                    median_genome_size=("genome_size", "median"),
                    gtdb_representative=("representative", "first"),
                    **{rank: (rank, "first") for rank in RANKS[:-1]})
               .reset_index())
    species.insert(0, "taxid", species["species"].map(keys).astype(str))
    species.insert(1, "scientific_name", species["species"])
    species.insert(2, "species_id", [f"{slug(s)}_{k}" for s, k in zip(species["species"], species["taxid"])])
    species["median_genome_size"] = species["median_genome_size"].astype(int)
    species = species.drop(columns="species")

    summary = pd.DataFrame({
        "gtdb_species": counts.index,
        "genomes_after_filters": counts.values,
        "eligible": counts.index.isin(eligible),
    }).merge(species[["scientific_name", "species_id", "n_selected"]],
             left_on="gtdb_species", right_on="scientific_name", how="left").drop(columns="scientific_name")
    summary["n_selected"] = summary["n_selected"].fillna(0).astype(int)

    order = ["Assembly Accession", "Organism Taxonomic ID", "Organism Name",
             "NCBI Organism Taxonomic ID", "species", "representative",
             "completeness", "contamination", "genome_size", *RANKS[:-1]]
    rest = [c for c in selected.columns if c not in order and c not in {"assembly_accession", "accession", "quality"}]
    selected = (selected[order + rest]
                .rename(columns={"species": "GTDB Species", "representative": "GTDB Representative",
                                 "completeness": "CheckM2 Completeness",
                                 "contamination": "CheckM2 Contamination",
                                 "genome_size": "GTDB Genome Size"})
                .sort_values(["Organism Taxonomic ID", "Assembly Accession"]))
    return selected, species, summary.sort_values("gtdb_species"), pd.DataFrame(steps, columns=["step", "genomes"]), thresholds


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--gtdb-metadata", required=True)
    p.add_argument("--candidates", required=True, help="NCBI long-read / de novo candidates.tsv")
    p.add_argument("--selected", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--steps", required=True)
    p.add_argument("--thresholds", required=True)
    p.add_argument("--min-completeness", type=float, default=95.0)
    p.add_argument("--max-contamination", type=float, default=5.0)
    p.add_argument("--min-genomes", type=int, default=50)
    p.add_argument("--max-genomes", type=int, default=50)
    p.add_argument("--seed", type=int, default=42)
    a = p.parse_args()

    gtdb = load_and_filter(a.gtdb_metadata, a.min_completeness, a.max_contamination)
    candidates = pd.read_csv(a.candidates, sep="\t", dtype=str)
    selected, species, summary, steps, thresholds = select(
        gtdb, candidates, a.min_genomes, a.max_genomes, a.seed)

    selected.to_csv(a.selected, sep="\t", index=False)
    species.to_csv(a.species, sep="\t", index=False)
    summary.to_csv(a.summary, sep="\t", index=False)
    steps.to_csv(a.steps, sep="\t", index=False)
    thresholds.to_csv(a.thresholds, sep="\t", index=False)
    print("\nGenomes left after each step:\n" + steps.to_string(index=False))
    print("\nSpecies kept at each minimum number of genomes:\n" + thresholds.to_string(index=False))
    print(f"\nSelected {len(selected):,} genomes from {len(species)} GTDB species "
          f"(>= {a.min_genomes}, max {a.max_genomes} each).")


if __name__ == "__main__":
    main()
