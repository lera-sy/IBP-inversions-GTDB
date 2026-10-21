"""Normalize circular chromosomes to dnaA and retain linear chromosomes unchanged."""

import argparse
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import pandas as pd

DNAAPLER_VERSION = "1.4.0"


def read_chromosome(path):
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith(">") or sum(line.startswith(">") for line in lines) != 1:
        raise ValueError("expected_one_fasta_sequence")
    header = lines[0][1:].split()
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    if not header or not sequence or set(sequence) - set("ACGTN"):
        raise ValueError("invalid_chromosome_fasta")
    return header[0], sequence


def check_transformation(original, normalized):
    if len(original) != len(normalized):
        raise ValueError("chromosome_length_changed")
    if original == normalized:
        return "unchanged", 0
    offset = (original + original).find(normalized)
    if 0 <= offset < len(original):
        return "rotation", offset
    reverse = original.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    offset = (reverse + reverse).find(normalized)
    if 0 <= offset < len(original):
        return "reverse_complement_and_rotation", offset
    raise ValueError("unexpected_sequence_change")


def normalize_one(row, args):
    row = row.copy()
    row.update({"Normalization Status": "REVIEW", "Normalization Notes": "",
                "Anchor Result": "", "Transformation": "", "Rotation Offset": "",
                "Normalized FASTA": "", "Dnaapler Log": "", "Dnaapler Version": ""})
    if row["Topology Status"] != "PASS":
        row["Normalization Notes"] = row["Topology Notes"]
        return row

    accession = row["Assembly Accession"]
    source = args.chromosomes / f"{accession}.fna"
    log_path = args.logs / f"{accession}.log"
    try:
        sequence_id, original = read_chromosome(source)
        if sequence_id != row["Chromosome Accession"]:
            raise ValueError("input_chromosome_accession_mismatch")
        if row["Topology"] not in {"circular", "linear"} or len(original) != int(row["NCBI Length"]) or len(original) != int(row["Chromosome Length"]):
            raise ValueError("input_topology_or_length_mismatch")

        if row["Topology"] == "linear":
            normalized = original
            transformation, offset = "unchanged_linear", 0
            row["Anchor Result"] = "not_applicable"
        else:
            row["Dnaapler Log"] = str(log_path)
            row["Dnaapler Version"] = DNAAPLER_VERSION
            # Temporary MMseqs2 files are removed; the log and summary are retained.
            summary_copy = args.logs / f"{accession}.tsv"
            summary_copy.unlink(missing_ok=True)
            with TemporaryDirectory(dir=args.logs) as temporary:
                work = Path(temporary) / "dnaapler"
                command = ["dnaapler", "all", "--db", "dnaa", "--autocomplete", "none",
                           "--input", str(source), "--output", str(work), "--prefix", accession,
                           "--threads", "1"]
                with log_path.open("w") as log:
                    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
                if result.returncode != 0:
                    raise RuntimeError(f"Dnaapler exited with code {result.returncode}; see {log_path}")
                summary_path = work / f"{accession}_all_reorientation_summary.tsv"
                output = work / f"{accession}_reoriented.fasta"
                if not summary_path.is_file() or not output.is_file():
                    raise ValueError("missing_dnaapler_output_or_summary")
                shutil.copyfile(summary_path, summary_copy)
                summary = pd.read_csv(summary_path, sep="\t", dtype=str, keep_default_na=False)
                if len(summary) != 1 or summary.iloc[0]["Contig"].split()[0] != sequence_id:
                    raise ValueError("dnaapler_summary_sequence_mismatch")
                anchor = summary.iloc[0]["Gene_Reoriented"]
                row["Anchor Result"] = anchor
                if anchor not in {"dnaA", "Contig_already_reoriented"}:
                    raise ValueError("dnaA_anchor_not_confirmed")
                if anchor == "dnaA" and "DNAA" not in summary.iloc[0]["Top_Hit"].upper():
                    raise ValueError("unexpected_anchor_hit")
                output_id, normalized = read_chromosome(output)
                if output_id != sequence_id:
                    raise ValueError("normalized_chromosome_accession_mismatch")
                transformation, offset = check_transformation(original, normalized)

        # Write the accepted chromosome with its assembly accession as the ID.
        target = args.normalized / f"{accession}.fna"
        with target.open("w") as handle:
            handle.write(f">{accession}\n")
            for start in range(0, len(normalized), 80):
                handle.write(normalized[start:start + 80] + "\n")
        row.update({"Normalization Status": "PASS", "Transformation": transformation,
                    "Rotation Offset": offset, "Normalized FASTA": str(target)})
    except ValueError as error:
        row["Normalization Notes"] = str(error)
    return row


def timed_normalize(row, args):
    """Measure actual per-genome time, without changing scientific outputs."""
    started = perf_counter()
    result = normalize_one(row, args)
    return result, perf_counter() - started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--chromosomes", type=Path, required=True)
    parser.add_argument("--normalized", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--passed", required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    if args.workers < 1:
        raise ValueError("--workers must be at least 1.")
    if not shutil.which("dnaapler") or not shutil.which("mmseqs"):
        raise RuntimeError("Activate the pipeline environment containing Dnaapler and MMseqs2.")
    version = subprocess.run(["dnaapler", "--version"], check=True, capture_output=True, text=True)
    if version.stdout.strip().split()[-1] != DNAAPLER_VERSION:
        raise RuntimeError(f"This script expects Dnaapler {DNAAPLER_VERSION} from environment.yaml.")
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str)
    topology = pd.read_csv(args.topology, sep="\t", dtype=str, keep_default_na=False)
    key = "Assembly Accession"
    if metadata.empty or metadata[key].duplicated().any() or topology[key].duplicated().any() or set(metadata[key]) != set(topology[key]):
        raise ValueError("Metadata and topology must contain the same unique QC-passing assemblies.")

    args.normalized.mkdir(parents=True, exist_ok=True)
    args.logs.mkdir(parents=True, exist_ok=True)
    # A fresh run must not retain normalized FASTAs from an older cohort.
    for old_file in args.normalized.glob("*.fna"):
        old_file.unlink()
    # Preserve deterministic report order while reporting jobs as they finish.
    # Dnaapler already runs one CPU thread per subprocess; --workers controls
    # how many independent chromosomes normalize concurrently.
    rows = [None] * len(topology)
    total_started = perf_counter()
    print(f"Normalizing {len(topology)} chromosomes with {args.workers} "
          "concurrent single-thread Dnaapler processes.", flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(timed_normalize, row, args): idx
            for idx, row in enumerate(topology.to_dict("records"))
        }
        try:
            for number, future in enumerate(as_completed(futures), 1):
                idx = futures[future]
                rows[idx], seconds = future.result()
                print(f"Normalized {rows[idx]['Assembly Accession']}: "
                      f"{rows[idx]['Normalization Status']} "
                      f"({seconds:.1f}s); {number}/{len(topology)} finished", flush=True)
        except Exception:
            for future in futures:
                future.cancel()
            raise
    print(f"Normalization wall time: {perf_counter() - total_started:.1f}s", flush=True)
    report = pd.DataFrame(rows)
    report.to_csv(args.report, sep="\t", index=False)
    passed = report.loc[report["Normalization Status"] == "PASS", key]
    metadata[metadata[key].isin(passed)].to_csv(args.passed, sep="\t", index=False)
    print(report["Normalization Status"].value_counts().to_string())


if __name__ == "__main__":
    main()
