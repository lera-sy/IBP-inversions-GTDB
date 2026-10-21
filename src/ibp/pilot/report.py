"""
Phase 4 pilot report: does the pipeline behave on contrasting species?

Reads Module 5 outputs (status.json + inversions.tsv per caller and query)
and writes, per species:

- callers.tsv   per caller: comparisons by status, failure fraction, calls,
                calls >= floor, calls per genome, call types, SyRI run time
- loci.tsv      per species: Leen calls merged into loci (overlapping or nested
                calls in one genome = one locus), calls per locus, locus sizes,
                and agreement with SyRI (a locus agrees if a SyRI call in the
                same genome has both breakpoints within the tolerance of the
                locus outer span)

Light: runs in seconds on a login node.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_intervals(text):
    out = []
    for part in str(text).split(";"):
        if part:
            a, b = map(int, part.split("-"))
            out.append((a, b))
    return out


def linear_segments(intervals, length):
    """Split circular intervals that wrap the origin into linear pieces."""
    segs = []
    for a, b in intervals:
        if b > a:
            segs.append((a, b))
        else:
            segs.append((a, length))
            if b > 0:
                segs.append((0, b))
    return segs


def merge_loci(calls, length):
    """Group calls whose reference intervals overlap into loci.

    Calls are joined when any of their pieces overlap; pieces meeting exactly
    at the origin are joined too. A locus's outer span on the circle is the
    chromosome minus the largest uncovered gap, so a locus across the origin
    is reported with start > end. Returns (start, end, n_calls) per locus.
    """
    pieces = [(a, b, i) for i, text in enumerate(calls)
              for a, b in linear_segments(parse_intervals(text), length)]
    if not pieces:
        return []
    parent = list(range(len(calls)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        parent[find(i)] = find(j)

    pieces.sort()
    reach, owner = -1, None
    for a, b, i in pieces:
        if owner is not None and a <= reach:
            union(i, owner)
        if b > reach:
            reach, owner = b, i
    at_zero = [i for a, b, i in pieces if a == 0]
    at_end = [i for a, b, i in pieces if b == length]
    for i in at_zero:
        for j in at_end:
            union(i, j)

    groups = {}
    for a, b, i in pieces:
        groups.setdefault(find(i), []).append((a, b))
    loci = []
    for root, segs in groups.items():
        segs.sort()
        covered = []
        for a, b in segs:
            if covered and a <= covered[-1][1]:
                covered[-1][1] = max(covered[-1][1], b)
            else:
                covered.append([a, b])
        # gaps between covered arcs, including the one across the origin
        gaps = [(covered[k][1], covered[k + 1][0]) for k in range(len(covered) - 1)]
        gaps.append((covered[-1][1], covered[0][0] + length))
        g0, g1 = max(gaps, key=lambda g: g[1] - g[0])
        start, end = g1 % length, g0 % length
        if end == 0:
            end = length
        n_calls = len({i for a, b, i in pieces if find(i) == root})
        loci.append((start, end, n_calls))
    return sorted(loci)


def span(start, end, length):
    return end - start if end > start else length - start + end


def circ(a, b, length):
    d = abs(a - b) % length
    return min(d, length - d)


def read_comparisons(root, species, cohort, caller):
    rows, calls = [], []
    for status_path in sorted((root / species / cohort / caller).glob("*/status.json")):
        status = json.loads(status_path.read_text())
        row = {"query": status.get("query_assembly"), "status": status.get("status"),
               "reason": status.get("reason", status.get("error", "")),
               "ref_length": status.get("reference_length_bp"),
               "seconds": (status.get("minimap_seconds") or 0) + (status.get("syri_seconds") or 0)}
        rows.append(row)
        inv = status_path.parent / "inversions.tsv"
        if status.get("status") == "success" and inv.is_file():
            df = pd.read_csv(inv, sep="\t", dtype=str)
            if len(df):
                calls.append(df)
    calls = pd.concat(calls, ignore_index=True) if calls else pd.DataFrame()
    return pd.DataFrame(rows), calls


def caller_summary(species, caller, comps, calls, floor):
    n = len(comps)
    counts = comps["status"].value_counts() if n else pd.Series(dtype=int)
    ok = int(counts.get("success", 0))
    row = {"species": species, "caller": caller, "comparisons": n,
           **{s: int(counts.get(s, 0)) for s in ["success", "unresolved", "unsupported", "failed"]},
           "fail_fraction": round(1 - ok / n, 3) if n else np.nan,
           "calls": len(calls)}
    if len(calls):
        length = calls["reported_length_bp"].astype(float)
        ref_len = calls["query_assembly"].map(dict(zip(comps["query"], comps["ref_length"]))).astype(float)
        row["calls_longer_arc"] = int((length > ref_len / 2).sum())
        length = np.minimum(length, ref_len - length)  # shorter-arc rule for every caller
        row["calls_ge_floor"] = int((length >= floor).sum())
        row["calls_per_success"] = round(len(calls) / ok, 2) if ok else np.nan
        row["median_call_kb"] = round(float(length.median()) / 1000, 1)
        if caller == "leen":
            for t, c in calls["call_type"].value_counts().items():
                row[f"type_{t}"] = int(c)
    if caller == "syri" and n:
        row["median_seconds"] = round(float(comps["seconds"].median()), 1)
        row["max_seconds"] = round(float(comps["seconds"].max()), 1)
    return row


def locus_summary(species, leen, syri, comps, floor, tolerance):
    lengths = dict(zip(comps["query"], comps["ref_length"]))
    rows = []
    if not len(leen):
        return {"species": species, "leen_loci": 0}
    for query, group in leen.groupby("query_assembly"):
        length = int(lengths.get(query) or 0)
        if not length:
            continue
        s = syri[syri["query_assembly"] == query] if len(syri) else pd.DataFrame()
        for a, b, n in merge_loci(group["ref_intervals_0"].tolist(), length):
            size = span(a, b, length)
            match = False
            for _, c in s.iterrows():
                # Unordered breakpoints: SyRI may report the other arc of the
                # circle, which has the same two breakpoints in swapped order.
                x, y = int(c["ref_start_0"]), int(c["ref_end_0"])
                same = circ(a, x, length) <= tolerance and circ(b, y, length) <= tolerance
                swapped = circ(a, y, length) <= tolerance and circ(b, x, length) <= tolerance
                if same or swapped:
                    match = True
                    break
            # Shorter-arc rule applied to the merged locus as well.
            rows.append({"query": query, "n_calls": n, "size": min(size, length - size),
                         "outer_span": size, "syri_match": match})
    loci = pd.DataFrame(rows)
    big = loci[loci["size"] >= floor]
    return {"species": species, "leen_calls": len(leen), "leen_loci": len(loci),
            "calls_per_locus": round(len(leen) / len(loci), 2) if len(loci) else np.nan,
            "loci_multi_call": int((loci["n_calls"] > 1).sum()),
            "loci_ge_floor": len(big),
            "loci_outer_span_gt_half": int((loci["outer_span"] > loci["size"]).sum()),
            "median_locus_kb": round(float(big["size"].median()) / 1000, 1) if len(big) else np.nan,
            "max_locus_kb": round(float(big["size"].max()) / 1000, 1) if len(big) else np.nan,
            "loci_matched_by_syri": int(big["syri_match"].sum()) if len(big) else 0,
            "fraction_matched": round(float(big["syri_match"].mean()), 2) if len(big) else np.nan}


def main():
    p = argparse.ArgumentParser(description="Phase 4 pilot report")
    p.add_argument("--module5-root", type=Path, required=True)
    p.add_argument("--species", nargs="+", required=True)
    p.add_argument("--cohort", default="selected_50_v1")
    p.add_argument("--floor-bp", type=int, default=6000)
    p.add_argument("--tolerance-bp", type=int, default=5000)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()

    a.out.mkdir(parents=True, exist_ok=True)
    caller_rows, locus_rows = [], []
    for sp in a.species:
        data = {c: read_comparisons(a.module5_root, sp, a.cohort, c) for c in ["leen", "syri"]}
        for caller, (comps, calls) in data.items():
            caller_rows.append(caller_summary(sp, caller, comps, calls, a.floor_bp))
        locus_rows.append(locus_summary(sp, data["leen"][1], data["syri"][1], data["leen"][0],
                                        a.floor_bp, a.tolerance_bp))
    callers = pd.DataFrame(caller_rows)
    loci = pd.DataFrame(locus_rows)
    callers.to_csv(a.out / "callers.tsv", sep="\t", index=False)
    loci.to_csv(a.out / "loci.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print("Per caller:\n" + callers.to_string(index=False))
        print("\nLeen loci (overlapping/nested calls merged):\n" + loci.to_string(index=False))


if __name__ == "__main__":
    main()
