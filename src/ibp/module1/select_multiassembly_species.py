"""
Second filter for IBP data, select (reproducible manner) species that have multiple assemblies

Only the species with >= 50 assemblies are kept.
Organism Taxonomic IDs are also resolved to specie-level via the NCBI taxonomy API.
"""

import argparse
import time

import pandas as pd
import requests

species_names = {}


def resolve_species_taxids(tax_ids, batch_size = 100):
    ids = sorted({str(int(float(t))) for t in tax_ids if pd.notna(t)})
    mapping = {}
    for i in range(0, len(ids), batch_size):
        batch = ids[i : i + batch_size]
        url = ("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/" + ",".join(batch) + "/dataset_report?page_size=1000")
        for attempt in range(5):
            r = requests.get(url)
            if r.status_code == 429 or r.status_code >= 500:
                print(f"API returned {r.status_code}, retrying in {2 ** attempt} seconds...")
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            break
        else:
            raise RuntimeError(f"NCBI API kept returning 429 for batch starting at index {i}")
        for report in r.json().get("reports", []):
            tax = report.get("taxonomy", {})
            sp = (tax.get("classification") or {}).get("species") or {}
            if sp.get("id"):
                mapping[str(tax["tax_id"])] = str(sp["id"])
                species_names[str(sp["id"])] = sp.get("name", "")
        time.sleep(0.35) # 3 requests per second
    print(f"Resolved {len(mapping)} / {len(ids)} taxonomic IDs to species-level.")
    return mapping


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("-s", "--summary", required=True)
parser.add_argument("-m", "--min-assemblies", type=int, default=50)
parser.add_argument("-n", "--assemblies-per-taxon", type=int, default=50)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

if args.assemblies_per_taxon > args.min_assemblies:
    raise ValueError("--assemblies-per-taxon cannot exceed --min-assemblies")

df = pd.read_csv(args.input, sep="\t", dtype=str)
print(f"Input: {len(df)} assemblies from {args.input}.")

accession_column = "Assembly Accession"
tax_id_column = "Organism Taxonomic ID"
biosample_column = "Assembly BioSample Accession"
required_columns = [accession_column, tax_id_column, biosample_column]
missing_columns = [c for c in required_columns if c not in df.columns]
if missing_columns:
    raise ValueError(f"Missing required columns: {missing_columns}")

df["Original Organism Taxonomic ID"] = df[tax_id_column]
species_map = resolve_species_taxids(df[tax_id_column])
df[tax_id_column] = df[tax_id_column].map(species_map)
# drop unresolved tax IDs
n_old = len(df)
df = df.dropna(subset=[tax_id_column])
df = df[~df[tax_id_column].map(species_names).fillna("").str.contains(r"\bsp\.|uncultured|unclassified|\bbacterium\b|\barchaeon\b", case=False)]
print(f"Dropped {n_old - len(df)} assemblies that do not have a proper species-level taxon.")

# Only keep one record per accession, then one record per biosample within each taxon (to avoid duplicates)
df = df.sort_values(accession_column).drop_duplicates(accession_column)
df["Biosample Selection Key"] = df[biosample_column].fillna("").str.strip()
df.loc[df["Biosample Selection Key"] == "", "Biosample Selection Key"] = df[accession_column]

candidate_counts = df.groupby(tax_id_column).size().rename("Candidate Assemblies")
unique_samples = df.drop_duplicates([tax_id_column, "Biosample Selection Key"])
sample_counts = unique_samples.groupby(tax_id_column).size().rename("Unique Biosamples")
eligible_taxa = sample_counts[sample_counts >= args.min_assemblies].index
eligible = unique_samples[unique_samples[tax_id_column].isin(eligible_taxa)]
if eligible.empty:
    raise RuntimeError("No taxa meet the minimum number of unique BioSamples.")

selected_groups = [
    group.sample(n=args.assemblies_per_taxon, random_state=args.seed)
    for _, group in eligible.groupby(tax_id_column, sort=True)
]
selected = (
    pd.concat(selected_groups, ignore_index=True)
    .sort_values([tax_id_column, accession_column])
    .drop(columns="Biosample Selection Key", errors="ignore")
)

selected_counts = selected.groupby(tax_id_column).size().rename("Selected Assemblies")
summary = pd.concat([candidate_counts, sample_counts, selected_counts], axis=1).fillna(0).reset_index()
for col in ["Candidate Assemblies", "Unique Biosamples", "Selected Assemblies"]:
    summary[col] = summary[col].astype(int)

selected.to_csv(args.output, sep="\t", index=False)
summary.to_csv(args.summary, sep="\t", index=False)

print(f"Selected {len(selected)} assemblies: {args.assemblies_per_taxon} each from {selected[tax_id_column].nunique()} taxa with >= {args.min_assemblies} unique BioSamples.")

print(f"Written to {args.output} and summary to {args.summary}.")

