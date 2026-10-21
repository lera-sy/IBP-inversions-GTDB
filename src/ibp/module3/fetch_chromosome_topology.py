"""Fetch chromosome topology for assemblies that passed genome QC."""

import argparse
import time

import pandas as pd
import requests


def fetch_summaries(accessions):
    params = {
        "db": "nuccore", "id": ",".join(accessions),
        "retmode": "json", "version": "2.0", "tool": "ibp_pipeline",
    }
    for delay in (1, 30, 60):
        time.sleep(delay)
        try:
            response = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params=params, timeout=60,
            )
            response.raise_for_status()
            data = response.json()
            if data.get("error") or "uids" not in data.get("result", {}):
                raise ValueError("NCBI returned an invalid sequence-summary response.")
            return data["result"]
        except (requests.RequestException, ValueError):
            if delay == 60:
                raise
            print("NCBI request failed; retrying after a pause.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    qc = pd.read_csv(args.qc, sep="\t", dtype=str)
    selected = qc.loc[qc["QC Status"] == "PASS", [
        "Assembly Accession", "Organism Taxonomic ID", "Chromosome Accession", "Chromosome Length",
    ]].copy()
    if selected.empty or selected.isna().any().any() or selected["Assembly Accession"].duplicated().any():
        raise ValueError("Expected QC-passing assemblies with unique accessions and chromosome details.")

    accessions = selected["Chromosome Accession"].drop_duplicates().tolist()
    summaries = {}
    for start in range(0, len(accessions), 200):
        batch = accessions[start:start + 200]
        result = fetch_summaries(batch)
        for uid in result["uids"]:
            record = result[uid]
            if record.get("accessionversion"):
                # Match the full accession version, never a newer version instead.
                summaries[record["accessionversion"]] = record
        print(f"Fetched summaries for {min(start + 200, len(accessions))}/{len(accessions)} chromosomes.", flush=True)

    rows = []
    for row in selected.to_dict("records"):
        record = summaries.get(row["Chromosome Accession"], {})
        topology = (record.get("topology") or "").lower()
        row["Topology"] = topology if topology in {"circular", "linear"} else "unknown"
        row["NCBI Length"] = record.get("slen", "")
        if not record:
            reason = "missing_exact_sequence_summary"
        elif not record.get("slen"):
            reason = "missing_ncbi_length"
        elif int(record["slen"]) != int(row["Chromosome Length"]):
            reason = "chromosome_length_mismatch"
        elif row["Topology"] == "unknown":
            reason = "unknown_topology"
        else:
            reason = ""
        row["Topology Status"] = "REVIEW" if reason else "PASS"
        row["Topology Notes"] = reason
        rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(args.output, sep="\t", index=False)
    print(result["Topology Status"].value_counts().to_string())


if __name__ == "__main__":
    main()
