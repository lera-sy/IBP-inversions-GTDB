"""Select the accepted chromosome with the lowest mean Mash distance per species."""

import argparse
import math
import subprocess
from pathlib import Path

import pandas as pd

from ibp.module3.normalize_chromosomes import read_chromosome


def rank_distances(path, assemblies, eligible=None):
    distances = {}
    for line in path.read_text().splitlines():
        first, second, value, _, _ = line.split("\t")
        pair = (Path(first).stem, Path(second).stem)
        distance = float(value)
        if pair in distances or not math.isfinite(distance) or not 0 <= distance <= 1:
            raise ValueError("Duplicate comparison or invalid Mash distance.")
        distances[pair] = distance
    expected = {(first, second) for first in assemblies for second in assemblies}
    if set(distances) != expected:
        raise ValueError("Mash distances do not match the complete accepted cohort.")
    for (first, second), distance in distances.items():
        if distance != distances[second, first] or (first == second and distance != 0):
            raise ValueError("Mash distances must be symmetric with zero self-distances.")

    rows = []
    for accession in sorted(assemblies):
        neighbors = [distances[accession, other] for other in sorted(assemblies) if other != accession]
        mean = sum(neighbors) / len(neighbors) if neighbors else None
        rows.append({"assembly_id": accession, "mean_mash_distance": mean})
    ranking = pd.DataFrame(rows).sort_values(["mean_mash_distance", "assembly_id"])
    # The reference must be circular (the callers need it); fall back to all if none is.
    allowed = set(eligible) if eligible else set(assemblies)
    ranking["eligible_reference"] = ranking["assembly_id"].isin(allowed)
    if not ranking["eligible_reference"].any():
        ranking["eligible_reference"] = True
    best = ranking[ranking["eligible_reference"]].iloc[0]["assembly_id"]
    ranking["is_reference"] = (ranking["assembly_id"] == best).map({True: "true", False: "false"})
    return ranking


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--normalization", required=True)
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--fastas", nargs="+", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    key = "Assembly Accession"
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str, keep_default_na=False)
    report = pd.read_csv(args.normalization, sep="\t", dtype=str, keep_default_na=False)
    if metadata[key].duplicated().any() or report[key].duplicated().any():
        raise ValueError("Assembly accessions must be unique in metadata and normalization report.")
    selected = metadata.loc[metadata["Organism Taxonomic ID"] == args.taxid, key].tolist()
    if not selected or not set(selected).issubset(set(report[key])):
        raise ValueError("Species must have accepted assemblies present in the normalization report.")
    details = report.set_index(key).loc[selected]
    if not (details["Normalization Status"].eq("PASS").all()
            and details["Organism Taxonomic ID"].eq(args.taxid).all()
            and details["Topology"].isin(["circular", "linear"]).all()):
        raise ValueError("Reference selection requires passing chromosomes from one species.")
    if len(args.fastas) != len(selected) or {path.stem for path in args.fastas} != set(selected):
        raise ValueError("FASTA inputs must match the accepted assemblies exactly.")
    for path in args.fastas:
        row = details.loc[path.stem]
        if path.resolve() != Path(row["Normalized FASTA"]).resolve():
            raise ValueError(f"FASTA path differs from the normalization report: {path}")
        sequence_id, sequence = read_chromosome(path)
        if sequence_id != path.stem or len(sequence) != int(row["Chromosome Length"]) or len(sequence) != int(row["NCBI Length"]):
            raise ValueError(f"FASTA ID or length differs from the normalization report: {path}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_dir / "chromosomes"
    sketch = args.output_dir / "chromosomes.msh"
    distances = args.output_dir / "mash_distances.tsv"
    subprocess.run(["mash", "sketch", "-k", "21", "-s", "10000", "-S", "42", "-p", "1",
                    "-o", str(prefix), *[str(path.resolve()) for path in sorted(args.fastas)]], check=True)
    with distances.open("w") as handle:
        subprocess.run(["mash", "dist", "-p", "1", str(sketch), str(sketch)], stdout=handle, check=True)
    circular = details.index[details["Topology"] == "circular"].tolist()
    ranking = rank_distances(distances, selected, circular)
    ranking.to_csv(args.output_dir / "reference_ranking.tsv", sep="\t", index=False)
    reference = ranking.loc[ranking["is_reference"] == "true", "assembly_id"].iloc[0]
    (args.output_dir / "reference.txt").write_text(reference + "\n")
    print(f"Species {args.taxid}: selected {reference} from {len(selected)} accepted chromosomes.")


if __name__ == "__main__":
    main()
