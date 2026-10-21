"""Verify downloaded assemblies and extract their largest annotated chromosome."""

import argparse
import json
from pathlib import Path

import pandas as pd


def read_fasta(path):
    records = []
    header = None
    parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts).upper()))
                header = line[1:].strip()
                if not header:
                    raise ValueError("empty_fasta_header")
                parts = []
            else:
                if header is None:
                    raise ValueError("sequence_before_fasta_header")
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts).upper()))
    if not records or any(not sequence for _, sequence in records):
        raise ValueError("empty_fasta_record")
    return records


def qc_assembly(accession, genomes, reports, chromosomes):
    row = {
        "Assembly Accession": accession, "QC Status": "REVIEW", "QC Notes": "",
        "FASTA Sequences": "", "Assembly Length": "", "Chromosomes": "", "Plasmids": "",
        "Chromosome Accession": "", "Chromosome Length": "", "GC Percent": "", "N Percent": "",
        "Excluded Chromosomes": "", "Excluded Plasmids": "", "Chromosome FASTA": "",
    }
    try:
        files = list((genomes / "ncbi_dataset/data" / accession).glob("*_genomic.fna"))
        if len(files) != 1:
            raise ValueError("missing_or_multiple_fasta_files")
        fasta = read_fasta(files[0])
        row["FASTA Sequences"] = len(fasta)
        row["Assembly Length"] = sum(len(sequence) for _, sequence in fasta)

        report_path = reports / "ncbi_dataset/data" / accession / "sequence_report.jsonl"
        if not report_path.is_file():
            raise ValueError("missing_sequence_report")
        report = [json.loads(line) for line in report_path.read_text().splitlines() if line.strip()]
        row["Chromosomes"] = sum(r.get("assignedMoleculeLocationType") == "Chromosome" for r in report)
        row["Plasmids"] = sum(r.get("assignedMoleculeLocationType") == "Plasmid" for r in report)
        if len(fasta) != len(report):
            raise ValueError("fasta_report_sequence_count_mismatch")

        # Both RefSeq and GenBank sequence IDs can identify the same record.
        sequence_ids = {}
        for index, record in enumerate(report):
            if record.get("assemblyAccession") != accession:
                raise ValueError("report_assembly_accession_mismatch")
            if record.get("assignedMoleculeLocationType") not in {"Chromosome", "Plasmid"}:
                raise ValueError("unclassified_sequence")
            ids = {record[key] for key in ("refseqAccession", "genbankAccession") if record.get(key)}
            if not ids or sequence_ids.keys() & ids:
                raise ValueError("missing_or_duplicate_report_sequence_id")
            for sequence_id in ids:
                sequence_ids[sequence_id] = index

        matched = set()
        chromosome_records = []
        plasmid_ids = []
        for header, sequence in fasta:
            sequence_id = header.split()[0]
            if sequence_id not in sequence_ids:
                raise ValueError("fasta_sequence_id_missing_from_report")
            index = sequence_ids[sequence_id]
            if index in matched:
                raise ValueError("duplicate_fasta_sequence")
            matched.add(index)
            record = report[index]
            if not record.get("length") or int(record["length"]) <= 0:
                raise ValueError("missing_or_invalid_report_length")
            if len(sequence) != int(record["length"]):
                raise ValueError("sequence_length_mismatch")
            if record["assignedMoleculeLocationType"] == "Chromosome":
                chromosome_records.append((header, sequence))
            else:
                plasmid_ids.append(sequence_id)

        if not chromosome_records:
            raise ValueError("no_annotated_chromosome")
        header, sequence = max(chromosome_records, key=lambda record: len(record[1]))
        if sum(len(seq) == len(sequence) for _, seq in chromosome_records) != 1:
            raise ValueError("largest_chromosome_length_tie")
        if set(sequence) - set("ACGTN"):
            raise ValueError("non_ACGTN_chromosome_bases")

        # Only the selected chromosome is written; the original assembly stays intact.
        target = chromosomes / f"{accession}.fna"
        with target.open("w") as handle:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
        row.update({
            "QC Status": "PASS", "Chromosome Accession": header.split()[0],
            "Chromosome Length": len(sequence),
            "GC Percent": round(100 * (sequence.count("G") + sequence.count("C")) / len(sequence), 4),
            "N Percent": round(100 * sequence.count("N") / len(sequence), 4),
            "Excluded Chromosomes": ",".join(h.split()[0] for h, _ in chromosome_records if h != header),
            "Excluded Plasmids": ",".join(plasmid_ids), "Chromosome FASTA": str(target),
        })
    except ValueError as error:
        row["QC Notes"] = str(error)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--genomes", type=Path, required=True)
    parser.add_argument("--reports", type=Path, required=True)
    parser.add_argument("--chromosomes", type=Path, required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--passed", required=True)
    args = parser.parse_args()

    selected = pd.read_csv(args.metadata, sep="\t", dtype=str)
    accessions = selected["Assembly Accession"]
    if selected.empty or accessions.isna().any() or accessions.duplicated().any():
        raise ValueError("Expected nonempty metadata with unique assembly accessions.")
    args.chromosomes.mkdir(parents=True, exist_ok=True)
    # Remove only previously extracted FASTAs, so a rerun cannot keep stale results.
    for old_file in args.chromosomes.glob("*.fna"):
        old_file.unlink()
    rows = []
    for number, accession in enumerate(accessions, 1):
        rows.append(qc_assembly(accession, args.genomes, args.reports, args.chromosomes))
        if number % 500 == 0:
            print(f"Checked {number}/{len(selected)} assemblies.", flush=True)
    qc = pd.DataFrame(rows)
    qc.insert(1, "Organism Taxonomic ID", selected["Organism Taxonomic ID"])
    qc.to_csv(args.qc, sep="\t", index=False)
    selected.loc[qc["QC Status"] == "PASS"].to_csv(args.passed, sep="\t", index=False)
    print(qc["QC Status"].value_counts().to_string())


if __name__ == "__main__":
    main()
