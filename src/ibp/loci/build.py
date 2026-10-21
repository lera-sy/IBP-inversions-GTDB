"""
Raw caller output -> calls -> loci (per genome) -> families (per species).

Applies the Phase 1 / Phase 4 definitions (config/config.yaml, docs/decisions.md):

calls     Successful comparisons only; circular query and reference only
          (linear chromosomes excluded). Every call, from either caller, is
          put on the SHORTER arc of the circle (SyRI often reports the longer
          one) and dropped if that arc is below the detection floor.
loci      Within one genome, calls from both callers whose reference intervals
          overlap or nest form one locus (one event). Locus size = shorter arc
          of the locus outer span (main) or the largest member call
          (sensitivity). Position class as the team: `local` if both outer
          breakpoints lie on one replication arm (origin at 0 after dnaapler,
          terminus at half the chromosome), else `spanning`.
tiers     Per locus, provisional until the raw-alignment check exists:
            high    both callers, with matching breakpoints (either order,
                    within the tolerance, circular distance)
            medium  both callers overlap, breakpoints do not match
            low     one caller only
families  Within a species, loci from different genomes with matching outer
          breakpoints (complete link) are one inversion family; each family is
          counted once per species. Family tier = best tier of its loci.

Outputs (one folder): calls.tsv, loci.tsv, families.tsv, species.tsv,
comparisons.tsv, pooled_sizes.tsv (family sizes for the empirical null).
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TIER_RANK = {"high": 3, "medium": 2, "low": 1}


# ---------------------------------------------------------------------------
# Circle geometry
# ---------------------------------------------------------------------------

def circ(a, b, length):
    d = abs(a - b) % length
    return min(d, length - d)


def arc(start, end, length):
    """Length of the arc from start to end going forward (wrap allowed)."""
    return end - start if end > start else length - start + end


def parse_intervals(text):
    return [tuple(map(int, p.split("-"))) for p in str(text).split(";") if p]


def linear_pieces(intervals, length):
    out = []
    for a, b in intervals:
        if b > a:
            out.append((a, b))
        else:
            out.append((a, length))
            if b > 0:
                out.append((0, b))
    return out


def outer_span(pieces, length):
    """Outer arc covering all pieces: the circle minus its largest gap."""
    pieces = sorted(pieces)
    covered = []
    for a, b in pieces:
        if covered and a <= covered[-1][1]:
            covered[-1][1] = max(covered[-1][1], b)
        else:
            covered.append([a, b])
    gaps = [(covered[k][1], covered[k + 1][0]) for k in range(len(covered) - 1)]
    gaps.append((covered[-1][1], covered[0][0] + length))
    g0, g1 = max(gaps, key=lambda g: g[1] - g[0])
    start, end = g1 % length, g0 % length
    return start, (end if end != 0 else length)


def group_overlapping(interval_texts, length):
    """Indices grouped by overlapping intervals; pieces meeting at the origin join.

    Returns a list of (start, end, member_indices) with the group's outer span.
    """
    pieces = [(a, b, i) for i, t in enumerate(interval_texts)
              for a, b in linear_pieces(parse_intervals(t), length)]
    if not pieces:
        return []
    parent = list(range(len(interval_texts)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    pieces.sort()
    reach, owner = -1, None
    for a, b, i in pieces:
        if owner is not None and a < reach:
            parent[find(i)] = find(owner)
        if b > reach:
            reach, owner = b, i
    zero = [i for a, b, i in pieces if a == 0]
    end = [i for a, b, i in pieces if b == length]
    for i in zero:
        for j in end:
            parent[find(i)] = find(j)
    groups = {}
    for a, b, i in pieces:
        groups.setdefault(find(i), {"pieces": [], "members": set()})
        groups[find(i)]["pieces"].append((a, b))
        groups[find(i)]["members"].add(i)
    out = [(*outer_span(g["pieces"], length), sorted(g["members"])) for g in groups.values()]
    return sorted(out)


def position_class(start, end, length):
    return "spanning" if (start < length / 2) != ((end % length) < length / 2) else "local"


def same_breakpoints(a, b, x, y, length, tolerance):
    """Breakpoint pairs match in either order (the other arc swaps them)."""
    return ((circ(a, x, length) <= tolerance and circ(b, y, length) <= tolerance)
            or (circ(a, y, length) <= tolerance and circ(b, x, length) <= tolerance))


# ---------------------------------------------------------------------------
# Reading Module 4 / 5 outputs
# ---------------------------------------------------------------------------

def read_package(prepared_root, species, cohort):
    g = pd.read_csv(Path(prepared_root) / species / cohort / "genomes.tsv", sep="\t", dtype=str)
    g["chromosome_length_bp"] = g["chromosome_length_bp"].astype(int)
    g["genome_length_bp"] = g["genome_length_bp"].astype(int)
    g["is_reference"] = g["is_reference"].str.lower() == "true"
    return g.set_index("assembly_id")


def read_caller(module5_root, species, cohort, caller):
    statuses, calls = [], []
    for path in sorted((Path(module5_root) / species / cohort / caller).glob("*/status.json")):
        st = json.loads(path.read_text())
        statuses.append({"species_id": species, "caller": caller,
                         "query_assembly": st.get("query_assembly"), "status": st.get("status"),
                         "reason": st.get("reason", st.get("error", ""))})
        inv = path.parent / "inversions.tsv"
        if st.get("status") == "success" and inv.is_file():
            df = pd.read_csv(inv, sep="\t", dtype=str)
            if len(df):
                calls.append(df)
    return statuses, (pd.concat(calls, ignore_index=True) if calls else pd.DataFrame())


# ---------------------------------------------------------------------------
# Calls -> loci -> families
# ---------------------------------------------------------------------------

def prepare_calls(raw, length, floor_bp):
    """Shorter-arc interval and size for every call; drop calls below the floor."""
    rows = []
    for _, c in raw.iterrows():
        if c["caller"] == "leen":
            pieces = linear_pieces(parse_intervals(c["ref_intervals_0"]), length)
            a, b = outer_span(pieces, length)
            text = c["ref_intervals_0"]
        else:
            a, b = int(c["ref_start_0"]), int(c["ref_end_0"])
            text = f"{a}-{b}"
        span = arc(a, b, length)
        if span > length / 2:   # report the shorter arc: same breakpoints, other way round
            a, b, span = b, a, length - span
            text = f"{a}-{b}"
        rows.append({"caller": c["caller"], "query_assembly": c["query_assembly"],
                     "call_id": c["call_id"], "call_type": c.get("call_type", ""),
                     "start": a, "end": b, "intervals": text, "size_bp": span,
                     "reported_length_bp": int(c["reported_length_bp"]),
                     "kept": span >= floor_bp})
    out = pd.DataFrame(rows)
    return out


def build_loci(calls, length, tolerance):
    rows = []
    for query, g in calls.groupby("query_assembly", sort=True):
        g = g.reset_index(drop=True)
        for k, (a, b, members) in enumerate(group_overlapping(g["intervals"].tolist(), length)):
            m = g.loc[members]
            leen, syri = m[m["caller"] == "leen"], m[m["caller"] == "syri"]
            agree = any(same_breakpoints(x.start, x.end, y.start, y.end, length, tolerance)
                        for x in leen.itertuples() for y in syri.itertuples())
            tier = "high" if agree else ("medium" if len(leen) and len(syri) else "low")
            span = arc(a, b, length)
            rows.append({"query_assembly": query, "locus": f"{query}:{k + 1}",
                         "start": a, "end": b,
                         "size_bp": min(span, length - span),
                         "largest_call_bp": int(m["size_bp"].max()),
                         "n_calls": len(m), "n_leen": len(leen), "n_syri": len(syri),
                         "breakpoints_agree": agree, "tier": tier,
                         "position_class": position_class(a, b, length),
                         "leen_call_types": ",".join(sorted(set(leen["call_type"].dropna()))),
                         "call_ids": ",".join(m["call_id"])})
    return pd.DataFrame(rows)


def build_families(loci, length, tolerance):
    """Complete-link grouping of loci across genomes by outer breakpoints."""
    if loci.empty:
        return pd.DataFrame(), loci.assign(family=pd.Series(dtype=str))
    order = loci.sort_values(["size_bp", "query_assembly", "start"]).index
    families = []
    for idx in order:
        a, b = loci.at[idx, "start"], loci.at[idx, "end"]
        for fam in families:
            if all(same_breakpoints(a, b, loci.at[j, "start"], loci.at[j, "end"], length, tolerance)
                   for j in fam):
                fam.append(idx)
                break
        else:
            families.append([idx])
    loci = loci.copy()
    rows = []
    for n, fam in enumerate(families, 1):
        fid = f"F{n:04d}"
        loci.loc[fam, "family"] = fid
        f = loci.loc[fam]
        tiers = f["tier"].value_counts()
        rows.append({"family": fid, "n_genomes": f["query_assembly"].nunique(),
                     "size_bp": float(np.median(f["size_bp"])),
                     "largest_call_bp": float(np.median(f["largest_call_bp"])),
                     "tier": max(f["tier"], key=TIER_RANK.get),
                     **{f"n_{t}": int(tiers.get(t, 0)) for t in TIER_RANK},
                     "position_class": f["position_class"].mode().iloc[0],
                     "start": int(f["start"].iloc[0]), "end": int(f["end"].iloc[0])})
    return pd.DataFrame(rows), loci


def process_species(module5_root, prepared_root, species, cohort, floor_bp, tolerance, max_failed):
    package = read_package(prepared_root, species, cohort)
    ref = package[package["is_reference"]]
    if len(ref) != 1:
        raise ValueError(f"{species}: expected one reference, found {len(ref)}")
    ref_id = ref.index[0]
    length = int(ref["chromosome_length_bp"].iloc[0])
    circular = package["topology"].eq("circular")

    statuses, raws = [], []
    for caller in ["leen", "syri"]:
        st, raw = read_caller(module5_root, species, cohort, caller)
        statuses += st
        if len(raw):
            raws.append(raw)
    comps = pd.DataFrame(statuses)
    comps["query_circular"] = comps["query_assembly"].map(circular).fillna(False)

    summary = {"species_id": species, "reference": ref_id, "reference_length_bp": length,
               "reference_circular": bool(circular.get(ref_id, False)),
               "genomes_accepted": len(package), "queries_circular": int(circular.sum() - circular.get(ref_id, False)),
               "median_genome_length_bp": int(package["genome_length_bp"].median()),
               "median_chromosome_length_bp": int(package["chromosome_length_bp"].median())}
    for caller in ["leen", "syri"]:
        c = comps[(comps["caller"] == caller) & comps["query_circular"]]
        summary[f"{caller}_success"] = int((c["status"] == "success").sum())
        summary[f"{caller}_fail_fraction"] = round(1 - (c["status"] == "success").mean(), 3) if len(c) else np.nan
    reasons = []
    if not summary["reference_circular"]:
        reasons.append("linear reference")
    for caller in ["leen", "syri"]:
        if not summary[f"{caller}_fail_fraction"] <= max_failed:
            reasons.append(f"{caller} failed on > {max_failed:.0%} of genomes")
    summary["excluded"] = "; ".join(reasons)

    # Only genomes where BOTH callers succeeded, so tiers mean the same everywhere.
    ok = comps[comps["query_circular"]].pivot_table(index="query_assembly", columns="caller",
                                                   values="status", aggfunc="first")
    both = set(ok.index[(ok.get("leen") == "success") & (ok.get("syri") == "success")])
    summary["genomes_both_callers"] = len(both)

    raw = pd.concat(raws, ignore_index=True) if raws else pd.DataFrame()
    if summary["excluded"] or raw.empty or not both:
        empty = pd.DataFrame()
        return summary, comps, empty, empty, empty
    raw = raw[raw["query_assembly"].isin(both)]
    calls = prepare_calls(raw, length, floor_bp)
    loci = build_loci(calls[calls["kept"]], length, tolerance) if calls["kept"].any() else pd.DataFrame()
    families, loci = build_families(loci, length, tolerance)
    summary.update({"calls": len(calls), "calls_kept": int(calls["kept"].sum()), "loci": len(loci),
                    "families": len(families)})
    if len(families):
        for t in TIER_RANK:
            summary[f"families_{t}"] = int((families["tier"] == t).sum())
    for df in (calls, loci, families):
        if len(df):
            df.insert(0, "species_id", species)
    return summary, comps, calls, loci, families


def main():
    p = argparse.ArgumentParser(description="Build calls, loci and families for all species")
    p.add_argument("--module5-root", type=Path, required=True)
    p.add_argument("--prepared-root", type=Path, required=True)
    p.add_argument("--species-table", type=Path, required=True)
    p.add_argument("--cohort", default="selected_50_v1")
    p.add_argument("--floor-bp", type=int, default=6000)
    p.add_argument("--tolerance-bp", type=int, default=5000)
    p.add_argument("--max-failed", type=float, default=0.20)
    p.add_argument("--only", nargs="*", help="species_id subset")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()

    species = pd.read_csv(a.species_table, sep="\t", dtype=str)
    ids = a.only or species["species_id"].tolist()
    out = {k: [] for k in ["species", "comparisons", "calls", "loci", "families"]}
    for sp in ids:
        if not (a.prepared_root / sp / a.cohort / "genomes.tsv").is_file():
            out["species"].append({"species_id": sp, "excluded": "no package"})
            continue
        res = process_species(a.module5_root, a.prepared_root, sp, a.cohort,
                              a.floor_bp, a.tolerance_bp, a.max_failed)
        out["species"].append(res[0])
        for key, df in zip(["comparisons", "calls", "loci", "families"], res[1:]):
            if len(df):
                out[key].append(df)

    a.out.mkdir(parents=True, exist_ok=True)
    tables = {"species": pd.DataFrame(out["species"])}
    for key in ["comparisons", "calls", "loci", "families"]:
        tables[key] = pd.concat(out[key], ignore_index=True) if out[key] else pd.DataFrame()
    info = species.set_index("species_id")
    sp = tables["species"].join(info[["scientific_name", "phylum"]], on="species_id")
    tables["species"] = sp
    for key, df in tables.items():
        df.to_csv(a.out / f"{key}.tsv", sep="\t", index=False)
    fam = tables["families"]
    if len(fam):
        keep = fam[fam["tier"].isin(["high", "medium"])]
        keep[["species_id", "family", "size_bp"]].to_csv(a.out / "pooled_sizes.tsv", sep="\t", index=False)

    included = sp[sp["excluded"].fillna("") == ""]
    print(f"Species: {len(sp)}; included {len(included)}; excluded {len(sp) - len(included)}")
    for _, r in sp[sp["excluded"].fillna("") != ""].iterrows():
        print(f"  excluded {r['species_id']}: {r['excluded']}")
    if len(fam):
        print(f"Calls kept (>= floor): {int(tables['calls']['kept'].sum()):,} of {len(tables['calls']):,}")
        print(f"Loci: {len(tables['loci']):,}; families: {len(fam):,}")
        print("Families by tier:\n" + fam["tier"].value_counts().reindex(list(TIER_RANK)).fillna(0).astype(int).to_string())
        print("Families by position class (high + medium):\n"
              + fam[fam["tier"] != "low"]["position_class"].value_counts().to_string())
        print(f"Species with >= 1 / >= 3 high+medium families: "
              f"{fam[fam['tier'] != 'low'].groupby('species_id').size().ge(1).sum()} / "
              f"{fam[fam['tier'] != 'low'].groupby('species_id').size().ge(3).sum()}")
    print(f"Written to {a.out}")


if __name__ == "__main__":
    main()
