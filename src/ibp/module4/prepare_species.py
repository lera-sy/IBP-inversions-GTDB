"""Publish one validated, portable genome package per species and cohort."""

import argparse
import errno
import hashlib
import json
import os
import re
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from ibp.module3.normalize_chromosomes import read_chromosome

COLUMNS = ["assembly_id", "fasta_path", "sequence_id", "chromosome_accession",
           "chromosome_length_bp", "genome_length_bp", "topology", "is_reference", "sha256"]


def validate_package(package):
    manifest = json.loads((package / "manifest.json").read_text())
    if (manifest["schema_version"] != "1.0" or manifest["analysis_scope"] != "primary_chromosome"
            or manifest["genomes_file"] != "genomes.tsv" or not manifest["scientific_name"]
            or type(manifest["taxid"]) is not int or manifest["taxid"] <= 0
            or type(manifest["selected_count"]) is not int or type(manifest["accepted_count"]) is not int
            or not manifest["cohort_id"] or "/" in manifest["cohort_id"]
            or not re.fullmatch(r"[a-z0-9_]+_" + str(manifest["taxid"]), manifest["species_id"])):
        raise ValueError("Invalid species package manifest.")
    genomes = pd.read_csv(package / "genomes.tsv", sep="\t", dtype=str, keep_default_na=False)
    if (list(genomes.columns) != COLUMNS or genomes.empty or genomes["assembly_id"].duplicated().any()
            or len(genomes) != manifest["accepted_count"] or manifest["selected_count"] < len(genomes)
            or not genomes["is_reference"].isin(["true", "false"]).all()):
        raise ValueError("Invalid genome table or package counts.")
    references = genomes.loc[genomes["is_reference"] == "true", "assembly_id"].tolist()
    if references != [manifest["reference_assembly_id"]]:
        raise ValueError("Package must contain exactly one matching reference.")
    if {path.name for path in (package / "genomes").iterdir()} != set(genomes["assembly_id"] + ".fna"):
        raise ValueError("Package FASTAs do not match the genome table.")
    for row in genomes.to_dict("records"):
        accession = row["assembly_id"]
        if (not re.fullmatch(r"GC[AF]_\d+\.\d+", accession)
                or row["fasta_path"] != f"genomes/{accession}.fna"
                or row["topology"] not in {"circular", "linear"} or not row["chromosome_accession"]):
            raise ValueError(f"Invalid assembly, path or topology: {accession}")
        fasta = package / row["fasta_path"]
        if fasta.is_symlink():
            raise ValueError(f"Package FASTAs must not be symlinks: {accession}")
        sequence_id, sequence = read_chromosome(fasta)
        length = int(row["chromosome_length_bp"])
        if (sequence_id != row["sequence_id"] or sequence_id != accession
                or len(sequence) != length or not 0 < length <= int(row["genome_length_bp"])):
            raise ValueError(f"FASTA ID or length mismatch: {accession}")
        if hashlib.sha256(fasta.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(f"FASTA checksum mismatch: {accession}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["selected", "accepted", "qc", "normalization", "reference", "species-table"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--fastas", nargs="+", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError("Published packages are not overwritten; choose a new cohort version.")
    if args.output_dir.name != args.cohort or args.output_dir.parent.name != args.species:
        raise ValueError("Package path must end with <species>/<cohort>.")

    names = pd.read_csv(args.species_table, sep="\t", dtype=str, keep_default_na=False)
    if names["species_id"].duplicated().any() or names["taxid"].duplicated().any():
        raise ValueError("Species IDs and TaxIDs must be unique in the species table.")
    species = names.set_index("species_id").loc[args.species]
    taxid = species["taxid"]
    key = "Assembly Accession"
    tables = {}
    for name in ["selected", "accepted", "qc", "normalization"]:
        table = pd.read_csv(getattr(args, name), sep="\t", dtype=str, keep_default_na=False)
        if table[key].duplicated().any():
            raise ValueError(f"Duplicate assembly accession in {name}.")
        tables[name] = table[table["Organism Taxonomic ID"] == taxid].set_index(key)
    selected, accepted, qc, normalization = [tables[name] for name in ["selected", "accepted", "qc", "normalization"]]
    if (selected.empty or accepted.empty or set(qc.index) != set(selected.index)
            or set(normalization.index) != set(qc.index[qc["QC Status"] == "PASS"])
            or set(accepted.index) != set(normalization.index[normalization["Normalization Status"] == "PASS"])):
        raise ValueError("Selected, QC, normalization and accepted memberships disagree.")
    fastas = {path.stem: path for path in args.fastas}
    if len(args.fastas) != len(accepted) or set(fastas) != set(accepted.index):
        raise ValueError("FASTA inputs must match the accepted assemblies exactly.")
    reference = Path(args.reference).read_text().strip()
    if reference not in accepted.index:
        raise ValueError("Reference is not an accepted assembly of this species.")

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".prepare-", dir=args.output_dir.parent) as temporary:
        package = Path(temporary) / "package"
        (package / "genomes").mkdir(parents=True)
        rows = []
        for accession in sorted(accepted.index):
            source = fastas[accession]
            record = normalization.loc[accession]
            checked = qc.loc[accession]
            if (source.resolve() != Path(record["Normalized FASTA"]).resolve()
                    or record["Chromosome Accession"] != checked["Chromosome Accession"]
                    or int(record["Chromosome Length"]) != int(checked["Chromosome Length"])
                    or int(record["NCBI Length"]) != int(checked["Chromosome Length"])):
                raise ValueError(f"QC and normalization disagree: {accession}")
            target = package / "genomes" / f"{accession}.fna"
            checksum = hashlib.sha256(source.read_bytes()).hexdigest()
            try:
                os.link(source, target)
            except OSError as error:
                if error.errno != errno.EXDEV:
                    raise
                shutil.copyfile(source, target)
            rows.append({"assembly_id": accession, "fasta_path": f"genomes/{accession}.fna",
                         "sequence_id": accession, "chromosome_accession": checked["Chromosome Accession"],
                         "chromosome_length_bp": int(checked["Chromosome Length"]),
                         "genome_length_bp": int(checked["Assembly Length"]), "topology": record["Topology"],
                         "is_reference": "true" if accession == reference else "false", "sha256": checksum})
        pd.DataFrame(rows, columns=COLUMNS).to_csv(package / "genomes.tsv", sep="\t", index=False)
        manifest = {"schema_version": "1.0", "species_id": args.species, "taxid": int(taxid),
                    "scientific_name": species["scientific_name"], "cohort_id": args.cohort,
                    "genomes_file": "genomes.tsv", "reference_assembly_id": reference,
                    "analysis_scope": "primary_chromosome", "selected_count": len(selected), "accepted_count": len(accepted)}
        (package / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        preparation = selected.copy()
        for table in [qc, normalization]:
            preparation = preparation.join(table.drop(columns=table.columns.intersection(preparation.columns)), how="left")
        preparation["Preparation Status"] = ["ACCEPTED" if accession in accepted.index else "EXCLUDED" for accession in preparation.index]
        preparation["Normalization Method"] = preparation["Topology"].map({
            "circular": "dnaapler_all_dnaa_autocomplete_none", "linear": "unchanged_linear"}).fillna("")
        preparation = preparation.drop(columns=["Chromosome FASTA", "Normalized FASTA", "Dnaapler Log"], errors="ignore")
        preparation.to_csv(package / "preparation.tsv", sep="\t")
        validate_package(package)
        package.rename(args.output_dir)
    print(f"Published {args.species}/{args.cohort}: {len(accepted)} accepted assemblies; reference {reference}.")


if __name__ == "__main__":
    main()
