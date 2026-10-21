#!/usr/bin/env python3
"""Final step: inversion size vs genome size across species, with phylogeny.

Reads the aggregated multi-species output of the inversion-only pipeline and:

1. Builds one event table per species from the two callers' candidate families.
   A Leen family and a SyRI family whose two breakpoints lie within --tolerance bp
   of each other (circular distance) are the same event. Each event gets a
   confidence tier:
     high    both callers found it
     medium  one caller only, and not "low"
     low     Leen only, and most of its calls are uncertain: with the refined
             caller (call_type column) that means labelled "shared_breakpoint"
             (indistinguishable from a moved-and-flipped segment); with the old
             caller it means overlapping other calls in the same genome
   Each family is counted once per species, so an inversion inherited by many
   genomes is one observation, not many.
2. Adds genome size per species: median total assembly length (whole genome) and
   median length of the analysed chromosome.
3. Maps every species to a GTDB species through its genome accessions, records
   CheckM2 quality, and takes the GTDB representative genome as the tree tip.
4. Prunes the GTDB bac120 tree to those tips and fits a PGLS of
   mean log10(event size) on log10(genome size), with Pagel's lambda estimated
   by maximum likelihood, and tests the slope against each null slope (S0, S1).
5. Fits a phylogenetic mixed model on the individual events (species effect
   split into a tree-shared and a species-unique part, plus within-species
   variation), the main model: it uses every event and the tree together.
6. Exploratory: the same models separately for local inversions (both
   breakpoints on one replication arm) and inversions spanning the origin or
   terminus (column position_class in events.tsv).

Usage:
  python final_inversion_model.py \
      --study results/multispecies/live95/final_reports \
      --cohorts results/multispecies/live95/generated_cohorts \
      --gtdb-metadata /path/bac120_metadata_r232.tsv.gz \
      --gtdb-tree /path/bac120_r232.tree \
      --pilot-root results/pilot \
      --out results/final_model
"""

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

OVERLAP_FLAG = "HIGHLY_OVERLAPPING_SAME_GENOME"
# name -> (tiers, position class or None for all)
ANALYSES = {
    "main_high_medium": (["high", "medium"], None),
    "sensitivity_all_tiers": (["high", "medium", "low"], None),
    "sensitivity_high_only": (["high"], None),
    "exploratory_local": (["high", "medium"], "local"),
    "exploratory_spanning": (["high", "medium"], "spanning"),
}


# ----------------------------------------------------------------------------
# 1. Events and confidence tiers
# ----------------------------------------------------------------------------

def circular_distance(a, b, length):
    d = abs(a - b) % length
    return min(d, length - d)


def breakpoint_cost(f1, f2, length):
    """Largest breakpoint distance under the best pairing of the two breakpoints."""
    a, b = f1
    c, d = f2
    direct = max(circular_distance(a, c, length), circular_distance(b, d, length))
    crossed = max(circular_distance(a, d, length), circular_distance(b, c, length))
    return min(direct, crossed)


def leen_uncertain_fraction(flags, leen_calls):
    """Per (species, event_id): fraction of member calls that are uncertain.

    Refined caller (calls have call_type): uncertain = "shared_breakpoint".
    Old caller (no call_type): uncertain = flagged as overlapping other calls.
    """
    if leen_calls is not None and "call_type" in leen_calls.columns:
        calls = leen_calls.copy()
        calls["uncertain"] = calls["call_type"] == "shared_breakpoint"
    else:
        calls = flags.copy()
        calls["uncertain"] = calls["qc_flags"].fillna("").str.contains(OVERLAP_FLAG)
    return (calls.groupby(["species", "event_id"])["uncertain"]
            .mean().rename("uncertain_fraction").reset_index())


def read_leen_calls(pilot_root, species_list):
    """Leen inversion_calls.tsv of every species, if the per-species results are available."""
    if pilot_root is None:
        return None
    parts = []
    for species in species_list:
        path = pilot_root / species / "callers" / "leen" / "inversion_calls.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"Missing Leen calls for {species}: {path}")
        calls = pd.read_csv(path, sep="\t")
        calls["species"] = species
        parts.append(calls)
    return pd.concat(parts, ignore_index=True)


def build_events(families, flags, ref_lengths, tolerance, overlap_majority, leen_calls=None):
    overlap = leen_uncertain_fraction(flags, leen_calls)
    events = []
    for species, fam in families.groupby("species"):
        length = ref_lengths[species]
        leen = fam[fam["caller"] == "leen"].copy()
        syri = fam[fam["caller"] == "minimap2_syri"].copy()
        leen["event_id"] = leen["family_id"].str.lstrip("F").astype(int)
        leen = leen.merge(overlap[overlap["species"] == species].drop(columns="species"),
                          on="event_id", how="left")
        leen["uncertain_fraction"] = leen["uncertain_fraction"].fillna(0.0)

        # Greedy one-to-one matching on breakpoint cost, cheapest pairs first.
        pairs = []
        for i, lrow in leen.iterrows():
            for j, srow in syri.iterrows():
                cost = breakpoint_cost((lrow.ref_start, lrow.ref_end),
                                       (srow.ref_start, srow.ref_end), length)
                if cost <= tolerance:
                    pairs.append((cost, i, j))
        used_l, used_s, matched = set(), set(), []
        for cost, i, j in sorted(pairs):
            if i not in used_l and j not in used_s:
                used_l.add(i)
                used_s.add(j)
                matched.append((i, j, cost))

        def record(tier, size, leen_row=None, syri_row=None, cost=None,
                   species=species, length=length):
            base = syri_row if syri_row is not None else leen_row
            events.append({
                "species": species,
                "tier": tier,
                "size_bp": int(size),
                "leen_family": None if leen_row is None else leen_row.family_id,
                "syri_family": None if syri_row is None else syri_row.family_id,
                "leen_size_bp": None if leen_row is None else int(leen_row.reference_relative_size_bp),
                "syri_size_bp": None if syri_row is None else int(syri_row.reference_relative_size_bp),
                "ref_start": int(base.ref_start),
                "ref_end": int(base.ref_end),
                "breakpoint_difference_bp": cost,
                "leen_uncertain_fraction": None if leen_row is None else round(leen_row.uncertain_fraction, 3),
                "carriers": int(max(leen_row.observed_carriers if leen_row is not None else 0,
                                    syri_row.observed_carriers if syri_row is not None else 0)),
                "reference_length_bp": length,
            })

        for i, j, cost in matched:
            # SyRI gives base-pair breakpoints; Leen's length is rounded to windows.
            record("high", syri.loc[j, "reference_relative_size_bp"],
                   leen.loc[i], syri.loc[j], int(cost))
        for j, srow in syri.iterrows():
            if j not in used_s:
                record("medium", srow.reference_relative_size_bp, syri_row=srow)
        for i, lrow in leen.iterrows():
            if i not in used_l:
                tier = "low" if lrow.uncertain_fraction > overlap_majority else "medium"
                record(tier, lrow.reference_relative_size_bp, leen_row=lrow)
    events = pd.DataFrame(events)
    events["log10_size"] = np.log10(events["size_bp"])
    return add_position_class(events)


def add_position_class(events, symmetric_tolerance=0.10):
    """Where an inversion sits relative to the replication origin and terminus.

    Genomes were rotated by dnaapler to start at dnaA, so the origin is at
    position 0; the terminus is taken as halfway round (an approximation).
    local     - both breakpoints on the same replication arm
    spanning  - breakpoints on opposite arms (the inversion spans the origin or
                the terminus); position_detail says whether it is symmetric
                (distances to the origin differ by < symmetric_tolerance of the
                chromosome) and around which point.
    """
    G = events["reference_length_bp"].astype(float)
    p1, p2 = events["ref_start"] % G, events["ref_end"] % G
    arm1, arm2 = p1 < G / 2, p2 < G / 2
    d1, d2 = np.where(arm1, p1, G - p1), np.where(arm2, p2, G - p2)
    opposite = arm1 != arm2
    symmetric = opposite & (np.abs(d1 - d2) < symmetric_tolerance * G)
    near_origin = (d1 + d2) / 2 < G / 4
    events = events.copy()
    events["position_class"] = np.where(opposite, "spanning", "local")
    events["position_detail"] = np.where(
        ~opposite, "local",
        np.where(symmetric, np.where(near_origin, "symmetric_origin", "symmetric_terminus"),
                 "asymmetric_opposite_arms"))
    return events


# ----------------------------------------------------------------------------
# 2. Genome size
# ----------------------------------------------------------------------------

LENGTH_COL = "Assembly Stats Total Sequence Length"
ACC_COL = "Assembly Accession"


def genome_sizes(cohort_dir, chromosome_sizes):
    rows = []
    for species, chrom in chromosome_sizes.groupby("species"):
        path = cohort_dir / f"{species}.selected.tsv"
        if not path.is_file():  # generated cohorts drop the _m<snapshot> suffix
            path = cohort_dir / f"{re.sub(r'_m[0-9a-f]{10}$', '', species)}.selected.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"Cohort table missing for {species}: {path}")
        cohort = pd.read_csv(path, sep="\t")
        cohort = cohort[cohort[ACC_COL].isin(chrom["accession"])]
        if cohort.empty:
            raise ValueError(f"No analysed genome of {species} found in {path}")
        rows.append({
            "species": species,
            "analysed_genomes": len(chrom),
            "genome_size_bp": float(cohort[LENGTH_COL].median()),
            "chromosome_size_bp": float(chrom["chromosome_length_bp"].median()),
            "accessions": list(chrom["accession"]),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# 3. GTDB mapping
# ----------------------------------------------------------------------------

def accession_key(accession):
    """GCF_000005845.2 / RS_GCF_000005845.2 / GCA_000005845.1 -> 000005845."""
    match = re.search(r"GC[AF]_(\d{9})", str(accession))
    return match.group(1) if match else None


def load_gtdb(metadata_path, wanted_keys):
    columns = ["accession", "gtdb_taxonomy", "gtdb_genome_representative",
               "checkm2_completeness", "checkm2_contamination"]
    parts = []
    for chunk in pd.read_csv(metadata_path, sep="\t", usecols=columns,
                             chunksize=200_000, low_memory=False):
        chunk["key"] = chunk["accession"].map(accession_key)
        parts.append(chunk[chunk["key"].isin(wanted_keys)])
    gtdb = pd.concat(parts, ignore_index=True)
    # Prefer the RefSeq (RS_) row if a genome is listed twice.
    gtdb["is_refseq"] = gtdb["accession"].str.startswith("RS_")
    gtdb = gtdb.sort_values("is_refseq", ascending=False).drop_duplicates("key")
    gtdb["gtdb_species"] = gtdb["gtdb_taxonomy"].str.extract(r"s__([^;]*)$")[0]
    return gtdb.set_index("key")


def map_to_gtdb(sizes, gtdb, min_completeness, max_contamination):
    rows = []
    for _, sp in sizes.iterrows():
        keys = [accession_key(a) for a in sp["accessions"]]
        found = gtdb.reindex([k for k in keys if k in gtdb.index])
        row = {"species": sp["species"], "genomes": len(keys), "genomes_in_gtdb": len(found)}
        if found.empty:
            rows.append({**row, "gtdb_species": None, "gtdb_representative": None})
            continue
        counts = found["gtdb_species"].value_counts()
        top = counts.index[0]
        reps = found.loc[found["gtdb_species"] == top, "gtdb_genome_representative"]
        passing = ((found["checkm2_completeness"] >= min_completeness)
                   & (found["checkm2_contamination"] <= max_contamination))
        rows.append({
            **row,
            "gtdb_species": top,
            "gtdb_species_agreement": round(counts.iloc[0] / len(found), 3),
            "other_gtdb_species": "; ".join(f"{s} ({n})" for s, n in counts.iloc[1:].items()),
            "gtdb_representative": reps.value_counts().index[0],
            "checkm2_pass_fraction": round(passing.mean(), 3),
            "median_checkm2_completeness": float(found["checkm2_completeness"].median()),
            "median_checkm2_contamination": float(found["checkm2_contamination"].median()),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# 4. Tree: parse, prune to tips, shared-path covariance
# ----------------------------------------------------------------------------

def parse_newick_tree(text):
    """Iterative Newick parser. Returns parent[], branch_length[], {tip name: node}.

    Handles quoted internal labels such as GTDB's '100.0:p__Bacillota'.
    """
    parent, length, tips = [-1], [0.0], {}
    stack = [0]          # open clades; node 0 is the root clade
    last = 0             # node that a following ':length' or label applies to
    i, n = 0, len(text)
    first = True
    while i < n:
        ch = text[i]
        if ch == "(":
            if first:
                first = False          # the root clade is node 0
            else:
                node = len(parent)
                parent.append(stack[-1])
                length.append(0.0)
                stack.append(node)
            i += 1
        elif ch == ",":
            i += 1
        elif ch == ")":
            last = stack.pop()
            i += 1
            # optional internal label (may be quoted)
            if i < n and text[i] == "'":
                j = text.index("'", i + 1)
                i = j + 1
            else:
                while i < n and text[i] not in ",():;":
                    i += 1
        elif ch == ":":
            j = i + 1
            while j < n and text[j] not in ",();":
                j += 1
            length[last] = float(text[i + 1:j])
            i = j
        elif ch == ";":
            break
        elif ch.isspace():
            i += 1
        else:
            if ch == "'":
                j = text.index("'", i + 1)
                label, i = text[i + 1:j], j + 1
            else:
                j = i
                while j < n and text[j] not in ",():;":
                    j += 1
                label, i = text[i:j].strip(), j
            node = len(parent)
            parent.append(stack[-1])
            length.append(0.0)
            tips[label] = node
            last = node
    return np.array(parent), np.array(length), tips


def tip_covariance(parent, length, tip_nodes):
    """Brownian covariance: shared root-to-MRCA path length for each pair of tips."""
    paths = []
    for node in tip_nodes:
        anc = {}
        chain = []
        while node != -1:
            chain.append(node)
            node = parent[node]
        # depth from root along the chain
        total = 0.0
        for node in reversed(chain):
            if parent[node] != -1:
                total += length[node]
            anc[node] = total
        paths.append(anc)
    k = len(tip_nodes)
    cov = np.zeros((k, k))
    for a in range(k):
        cov[a, a] = paths[a][tip_nodes[a]]
        for b in range(a + 1, k):
            shared = set(paths[a]).intersection(paths[b])
            mrca_depth = max(paths[a][x] for x in shared)
            cov[a, b] = cov[b, a] = mrca_depth
    return cov


def pruned_newick(parent, length, tips, names):
    """Newick of the subtree spanned by the given tips (unary nodes collapsed)."""
    nodes = {tips[nm]: nm for nm in names}
    keep = set()
    for node in nodes:
        while node != -1 and node not in keep:
            keep.add(node)
            node = parent[node]
    children = {}
    for node in keep:
        p = parent[node]
        if p != -1:
            children.setdefault(p, []).append(node)
    root = next(nd for nd in keep if parent[nd] == -1)

    def render(node, extra=0.0):
        kids = children.get(node, [])
        while len(kids) == 1 and node not in nodes:      # collapse unary chain
            extra += length[node] if parent[node] != -1 else 0.0
            node = kids[0]
            kids = children.get(node, [])
        bl = (length[node] if parent[node] != -1 else 0.0) + extra
        if node in nodes:
            return f"{nodes[node]}:{bl:.6f}"
        return "(" + ",".join(render(c) for c in kids) + f"):{bl:.6f}"

    return render(root) + ";"


# ----------------------------------------------------------------------------
# 5. PGLS with Pagel's lambda
# ----------------------------------------------------------------------------

def gls(y, X, V):
    Vi = np.linalg.inv(V)
    XtVi = X.T @ Vi
    cov_beta_unscaled = np.linalg.inv(XtVi @ X)
    beta = cov_beta_unscaled @ XtVi @ y
    r = y - X @ beta
    n = len(y)
    rss = float(r @ Vi @ r)
    _sign, logdet = np.linalg.slogdet(V)
    loglik = -0.5 * (n * math.log(2 * math.pi * rss / n) + logdet + n)
    return beta, cov_beta_unscaled, rss, loglik


def lambda_matrix(C, lam):
    V = C * lam
    np.fill_diagonal(V, np.diag(C))
    return V


def fit_pgls(y, x, C, null_slopes, fixed_lambda=None):
    X = np.column_stack([np.ones_like(x), x])
    n, k = X.shape
    if fixed_lambda is None:
        res = optimize.minimize_scalar(lambda lam: -gls(y, X, lambda_matrix(C, lam))[3],
                                       bounds=(0.0, 1.0), method="bounded")
        lam = float(res.x)
        # compare against the bounds (the optimum is often exactly 0 or 1)
        for edge in (0.0, 1.0):
            if gls(y, X, lambda_matrix(C, edge))[3] > gls(y, X, lambda_matrix(C, lam))[3]:
                lam = edge
    else:
        lam = fixed_lambda
    V = C if lam is None else lambda_matrix(C, lam)
    beta, cov_u, rss, loglik = gls(y, X, V)
    s2 = rss / (n - k)
    se = np.sqrt(np.diag(cov_u) * s2)
    tcrit = stats.t.ppf(0.975, n - k)
    out = {
        "n_species": n, "lambda": round(lam, 4) if lam is not None else None,
        "intercept": beta[0], "slope": beta[1], "slope_se": se[1],
        "slope_ci_low": beta[1] - tcrit * se[1], "slope_ci_high": beta[1] + tcrit * se[1],
        "loglik": loglik,
    }
    # Pseudo R2 against an intercept-only GLS on the same covariance.
    _, _, rss0, _ = gls(y, X[:, :1], V)
    out["r2_gls"] = 1 - rss / rss0
    for s in null_slopes:
        t = (beta[1] - s) / se[1]
        out[f"p_slope_vs_{s:g}"] = 2 * stats.t.sf(abs(t), n - k)
    return out


def fit_pmm(ybar, n, ssw, x, C, null_slopes, phylogenetic=True):
    """Phylogenetic mixed model on individual events, fitted by maximum likelihood.

    y_ij = a + b*x_i + p_i + s_i + e_ij for event j of species (tip) i, with
      p ~ N(0, s2p * C)  species effect shared along the tree (C = tree covariance)
      s ~ N(0, s2s * I)  species effect unique to each species
      e ~ N(0, s2e)      variation between events of the same species.
    Because x is a species-level value, the likelihood splits exactly into
    within-species deviations (depend on s2e only) and species means
    ybar_i ~ N(a + b*x_i, s2p*C + s2s*I + s2e*diag(1/n_i)), so every event is
    used without building an events x events matrix.
    phylogenetic=False fixes s2p = 0 (plain mixed model with a species effect).
    """
    k, N = len(ybar), int(n.sum())
    X = np.column_stack([np.ones(k), x])
    Cn = C / np.mean(np.diag(C))           # tree covariance scaled to mean variance 1
    inv_n = np.diag(1.0 / n)

    def unpack(theta):
        s2p = np.exp(theta[0]) if phylogenetic else 0.0
        return s2p, np.exp(theta[-2]), np.exp(theta[-1])

    def negloglik(theta, want=False):
        s2p, s2s, s2e = unpack(theta)
        V = s2p * Cn + s2s * np.eye(k) + s2e * inv_n
        try:
            L = np.linalg.cholesky(V)
        except np.linalg.LinAlgError:
            return 1e12
        Li = np.linalg.inv(L)
        Xw, yw = Li @ X, Li @ ybar
        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        r = yw - Xw @ beta
        ll_between = -0.5 * (k * np.log(2 * np.pi) + 2 * np.log(np.diag(L)).sum() + r @ r)
        ll_within = -0.5 * ((N - k) * np.log(2 * np.pi * s2e) + ssw / s2e)
        if want:
            cov = np.linalg.inv(Xw.T @ Xw)
            return beta, cov, ll_between + ll_within, (s2p, s2s, s2e)
        return -(ll_between + ll_within)

    total_var = np.var(ybar) + ssw / max(N - k, 1)
    start = np.log([total_var / 3] * (3 if phylogenetic else 2))
    best = None
    for shift in (0.0, -2.0, 2.0):          # a few starts: variance components can sit on a boundary
        res = optimize.minimize(negloglik, start + shift, method="L-BFGS-B",
                                bounds=[(-25, 5)] * len(start))
        if best is None or res.fun < best.fun:
            best = res
    beta, cov, loglik, (s2p, s2s, s2e) = negloglik(best.x, want=True)
    se = np.sqrt(np.diag(cov))
    tcrit = stats.t.ppf(0.975, k - 2)
    out = {"n_species": k, "n_events": N, "lambda": None,
           "intercept": beta[0], "slope": beta[1], "slope_se": se[1],
           "slope_ci_low": beta[1] - tcrit * se[1], "slope_ci_high": beta[1] + tcrit * se[1],
           "loglik": loglik, "var_tree": s2p, "var_species": s2s, "var_within": s2e,
           "tree_share_of_species_variance": s2p / (s2p + s2s) if s2p + s2s > 0 else np.nan}
    for s0 in null_slopes:
        t = (beta[1] - s0) / se[1]
        out[f"p_slope_vs_{s0:g}"] = 2 * stats.t.sf(abs(t), k - 2)
    return out


def tip_event_summary(events, tiers, table, merged):
    """Per tree tip: number of events, mean log10 size, within-tip sum of squares."""
    tip_of = dict(zip(table["species"], table["gtdb_representative"]))
    ev = events[events["tier"].isin(tiers) & events["species"].isin(tip_of)].copy()
    ev["tip"] = ev["species"].map(tip_of)
    g = ev.groupby("tip")["log10_size"]
    summary = pd.DataFrame({"n": g.size(), "ybar": g.mean(),
                            "ssw": g.apply(lambda v: float(((v - v.mean()) ** 2).sum()))})
    return summary.reindex(merged["gtdb_representative"]).reset_index(drop=True)


def species_table(events, tiers, sizes_mapped, size_column):
    sel = events[events["tier"].isin(tiers)]
    agg = (sel.groupby("species")
           .agg(n_events=("size_bp", "size"), mean_log10_size=("log10_size", "mean"),
                median_size_bp=("size_bp", "median"))
           .reset_index())
    table = sizes_mapped.merge(agg, on="species", how="left")
    table["n_events"] = table["n_events"].fillna(0).astype(int)
    table["log10_genome_size"] = np.log10(table[size_column])
    return table


def merge_by_tip(table):
    """NCBI species that GTDB lumps into one species share a tip: pool them."""
    def pool(g):
        w = g["n_events"].clip(lower=1)
        return pd.Series({
            "ncbi_species": "; ".join(g["species"]),
            "gtdb_species": g["gtdb_species"].iloc[0],
            "n_events": int(g["n_events"].sum()),
            "mean_log10_size": np.average(g["mean_log10_size"], weights=w),
            "log10_genome_size": np.average(g["log10_genome_size"], weights=w),
        })
    return (table.groupby("gtdb_representative", group_keys=False)
            .apply(pool, include_groups=False).reset_index())


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--study", type=Path, required=True,
                   help="Folder with the aggregated tables (candidate_families.tsv, ...)")
    p.add_argument("--cohorts", type=Path, required=True, help="generated_cohorts folder")
    p.add_argument("--gtdb-metadata", type=Path, required=True)
    p.add_argument("--gtdb-tree", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--tolerance", type=int, default=5000,
                   help="Max breakpoint distance (bp) for a Leen and a SyRI family to be one event")
    p.add_argument("--pilot-root", type=Path,
                   help="results/pilot folder; lets the refined caller's call_type set the tiers")
    p.add_argument("--overlap-majority", type=float, default=0.5,
                   help="Leen-only family is 'low' if more than this fraction of its calls are uncertain")
    p.add_argument("--null-slopes", default="0,1", help="Slopes to test against (S0,S1)")
    p.add_argument("--min-events", type=int, default=1, help="Min events for a species to enter a model")
    p.add_argument("--min-gtdb-agreement", type=float, default=0.8,
                   help="Drop species whose genomes split across GTDB species below this agreement")
    p.add_argument("--min-checkm2-pass", type=float, default=0.8,
                   help="Drop species with a lower fraction of genomes passing CheckM2")
    p.add_argument("--checkm2-completeness", type=float, default=95.0)
    p.add_argument("--checkm2-contamination", type=float, default=5.0)
    args = p.parse_args(argv)
    null_slopes = [float(s) for s in args.null_slopes.split(",")]
    args.out.mkdir(parents=True, exist_ok=True)

    families = pd.read_csv(args.study / "candidate_families.tsv", sep="\t")
    flags = pd.read_csv(args.study / "leen_call_flags.tsv", sep="\t")
    qc = pd.read_csv(args.study / "qc_summary.tsv", sep="\t")
    chrom = pd.read_csv(args.study / "genome_sizes.tsv", sep="\t")
    ref_lengths = (qc[qc["metric"] == "reference_length_bp"]
                   .set_index("species")["value"].astype(int).to_dict())

    print("1/4 Building events and confidence tiers", flush=True)
    leen_calls = read_leen_calls(args.pilot_root, sorted(families["species"].unique()))
    print("   tiers from " + ("refined caller call types" if leen_calls is not None
                              and "call_type" in leen_calls.columns else "overlap flags"), flush=True)
    events = build_events(families, flags, ref_lengths, args.tolerance, args.overlap_majority,
                          leen_calls)
    events.to_csv(args.out / "events.tsv", sep="\t", index=False)
    print(events.groupby("tier").size().to_string(), flush=True)

    print("2/4 Genome sizes", flush=True)
    sizes = genome_sizes(args.cohorts, chrom)

    print("3/4 GTDB mapping (reading metadata, takes a minute)", flush=True)
    keys = {accession_key(a) for accs in sizes["accessions"] for a in accs}
    gtdb = load_gtdb(args.gtdb_metadata, keys)
    mapping = map_to_gtdb(sizes, gtdb, args.checkm2_completeness, args.checkm2_contamination)
    sizes_mapped = sizes.drop(columns="accessions").merge(mapping, on="species")

    parent, length, tips = parse_newick_tree(args.gtdb_tree.read_text())
    sizes_mapped["in_tree"] = sizes_mapped["gtdb_representative"].isin(tips)
    keep = (sizes_mapped["in_tree"]
            & (sizes_mapped["gtdb_species_agreement"] >= args.min_gtdb_agreement)
            & (sizes_mapped["checkm2_pass_fraction"] >= args.min_checkm2_pass))
    sizes_mapped["excluded_reason"] = np.select(
        [~sizes_mapped["in_tree"],
         sizes_mapped["gtdb_species_agreement"] < args.min_gtdb_agreement,
         sizes_mapped["checkm2_pass_fraction"] < args.min_checkm2_pass],
        ["representative_not_in_tree", "splits_across_gtdb_species", "low_checkm2_quality"],
        default="")
    sizes_mapped.to_csv(args.out / "species_gtdb_mapping.tsv", sep="\t", index=False)
    print(f"   {keep.sum()}/{len(sizes_mapped)} species usable; excluded: "
          f"{sizes_mapped.loc[~keep, 'excluded_reason'].value_counts().to_dict()}", flush=True)
    usable = sizes_mapped[keep]

    print("4/4 PGLS", flush=True)
    results = []
    tables = []
    for size_column in ("genome_size_bp", "chromosome_size_bp"):
        for name, (tiers, position) in ANALYSES.items():
            subset = events if position is None else events[events["position_class"] == position]
            if size_column == "chromosome_size_bp" and name != "main_high_medium":
                continue
            t = species_table(subset, tiers, usable, size_column)
            t = t[t["n_events"] >= args.min_events]
            merged = merge_by_tip(t)
            if len(merged) < 4:
                print(f"   {name} / {size_column}: only {len(merged)} species, skipped", flush=True)
                continue
            C = tip_covariance(parent, length, [tips[r] for r in merged["gtdb_representative"]])
            y = merged["mean_log10_size"].to_numpy()
            x = merged["log10_genome_size"].to_numpy()
            label = {"analysis": name, "genome_size": size_column}
            results.append({**label, "model": "PGLS_pagel_lambda", **fit_pgls(y, x, C, null_slopes)})
            results.append({**label, "model": "PGLS_brownian", **fit_pgls(y, x, C, null_slopes, 1.0)})
            ols = fit_pgls(y, x, np.eye(len(y)), null_slopes, 1.0)
            ols["lambda"] = None
            results.append({**label, "model": "OLS_no_phylogeny", **ols})
            per_tip = tip_event_summary(subset, tiers, t, merged)
            n, ssw = per_tip["n"].to_numpy(float), float(per_tip["ssw"].sum())
            results.append({**label, "model": "PMM_phylogenetic_mixed",
                            **fit_pmm(per_tip["ybar"].to_numpy(), n, ssw, x, C, null_slopes)})
            results.append({**label, "model": "MM_species_effect_only",
                            **fit_pmm(per_tip["ybar"].to_numpy(), n, ssw, x, C, null_slopes, False)})
            merged.insert(0, "genome_size", size_column)
            merged.insert(0, "analysis", name)
            tables.append(merged)

    if not results:
        sys.exit("Not enough species for any model.")
    res = pd.DataFrame(results)
    res.to_csv(args.out / "pgls_results.tsv", sep="\t", index=False, float_format="%.6g")
    model_data = pd.concat(tables)
    model_data.to_csv(args.out / "model_species_table.tsv", sep="\t", index=False)
    main_tips = model_data.loc[(model_data["analysis"] == "main_high_medium")
                               & (model_data["genome_size"] == "genome_size_bp"),
                               "gtdb_representative"]
    (args.out / "species_tree_pruned.nwk").write_text(
        pruned_newick(parent, length, tips, list(main_tips)) + "\n")

    show = ["analysis", "genome_size", "model", "n_species", "lambda", "slope",
            "slope_ci_low", "slope_ci_high"] + [c for c in res.columns if c.startswith("p_slope")]
    print(res[show].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    plot(model_data, res, args.out)
    print(f"Written to {args.out}")


def plot(model_data, res, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed: no plot made")
        return
    d = model_data[(model_data["analysis"] == "main_high_medium")
                   & (model_data["genome_size"] == "genome_size_bp")]
    fit = res[(res["analysis"] == "main_high_medium") & (res["genome_size"] == "genome_size_bp")
              & (res["model"] == "PGLS_pagel_lambda")].iloc[0]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.scatter(d["log10_genome_size"], d["mean_log10_size"],
               s=12 + 4 * np.sqrt(d["n_events"]), alpha=0.7, color="#2a6f97")
    xs = np.linspace(d["log10_genome_size"].min(), d["log10_genome_size"].max(), 50)
    ax.plot(xs, fit["intercept"] + fit["slope"] * xs, color="#c44536",
            label=f"PGLS slope {fit['slope']:.2f} "
                  f"[{fit['slope_ci_low']:.2f}, {fit['slope_ci_high']:.2f}], λ={fit['lambda']:.2f}")
    ax.set_xlabel("log10 genome size (bp)")
    ax.set_ylabel("mean log10 inversion size (bp)")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.set_title("Inversion size vs genome size (high + medium confidence)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "inversion_vs_genome_size.png", dpi=200)


if __name__ == "__main__":
    main()
