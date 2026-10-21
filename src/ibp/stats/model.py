"""
Phase 7: does inversion size scale with genome size?

y = log10(family size, kb), x = log10(genome size, Mb), one family counted once
per species (results of src/ibp/loci/build.py), GTDB tree for relatedness.

Models (fitting code reused from the team's final_inversion_model.py):
- PMM   phylogenetic mixed model on every family (main): species effect split
        into a tree-shared and a species-unique part, plus within-species
        variation.
- MM    same without the tree (species effect only).
- PGLS  one geometric mean per species (species with >= --min-pgls families),
        Pagel's lambda by maximum likelihood (sensitivity).

Each slope is tested against 0, 1 (proportional scaling, S1) and S0, the slope
the measurement rules alone produce on these species (geometry null with the
real pooled sizes, real chromosome lengths and real family counts; see
src/ibp/nullmodel/geometry.py). Local and spanning inversions are tested
against their own S0.

Robustness: main PMM refitted leaving out one phylum at a time.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ibp.analysis_legacy.final_inversion_model import (fit_pgls, fit_pmm, parse_newick_tree,
                                                       tip_covariance)

# name: (tiers, position class or None, size column, genome size column)
ANALYSES = {
    "main": (["high", "medium"], None, "size_bp", "median_genome_length_bp"),
    "sens_high_only": (["high"], None, "size_bp", "median_genome_length_bp"),
    "sens_all_tiers": (["high", "medium", "low"], None, "size_bp", "median_genome_length_bp"),
    "sens_largest_call": (["high", "medium"], None, "largest_call_bp", "median_genome_length_bp"),
    "sens_chromosome_length": (["high", "medium"], None, "size_bp", "median_chromosome_length_bp"),
    "expl_local": (["high", "medium"], "local", "size_bp", "median_genome_length_bp"),
    "expl_spanning": (["high", "medium"], "spanning", "size_bp", "median_genome_length_bp"),
}


def null_input(families, species, out_dir):
    """Tables for the empirical geometry null (S0)."""
    hm = families[families["tier"].isin(["high", "medium"])]
    counts = hm.groupby("species_id").size()
    tab = species.assign(n_events=species["species_id"].map(counts).fillna(0).astype(int))
    tab = tab[tab["n_events"] > 0]
    tab = tab.rename(columns={"reference_length_bp": "chromosome_length_bp",
                              "median_genome_length_bp": "genome_length_bp"})
    tab[["species_id", "chromosome_length_bp", "genome_length_bp", "n_events"]].to_csv(
        out_dir / "null_species.tsv", sep="\t", index=False)
    hm[["species_id", "family", "size_bp"]].to_csv(out_dir / "null_sizes.tsv", sep="\t", index=False)


def read_s0(null_path):
    """S0 per subset from the empirical fixed-size null (events fit)."""
    if not null_path or not Path(null_path).is_file():
        return {}
    t = pd.read_csv(null_path, sep="\t")
    t = t[(t["size_distribution"] == "empirical_pooled_sizes") & (t["fit"] == "events")]
    return {r.subset: (r["slope_mean"], r["slope_2.5"], r["slope_97.5"]) for _, r in t.iterrows()}


def species_data(families, species, tiers, position, size_col, x_col):
    ev = families[families["tier"].isin(tiers)]
    if position:
        ev = ev[ev["position_class"] == position]
    ev = ev.assign(y=np.log10(ev[size_col].astype(float) / 1000))
    g = ev.groupby("species_id")["y"]
    d = pd.DataFrame({"n": g.size(), "ybar": g.mean(),
                      "ssw": g.apply(lambda v: float(((v - v.mean()) ** 2).sum()))})
    d = d.join(species.set_index("species_id")[[x_col, "phylum"]], how="inner")
    d["x"] = np.log10(d[x_col].astype(float) / 1e6)
    return d


def fit_all(d, tree, min_pgls, s0):
    parent, length, tips = tree
    d = d[d.index.isin(tips)]
    nulls = [0.0, 1.0] + ([s0] if s0 is not None else [])
    C = tip_covariance(parent, length, [tips[s] for s in d.index])
    out = {}
    out["PMM"] = fit_pmm(d["ybar"].to_numpy(), d["n"].to_numpy(float), float(d["ssw"].sum()),
                         d["x"].to_numpy(), C, nulls, phylogenetic=True)
    out["MM"] = fit_pmm(d["ybar"].to_numpy(), d["n"].to_numpy(float), float(d["ssw"].sum()),
                        d["x"].to_numpy(), C, nulls, phylogenetic=False)
    m = d[d["n"] >= min_pgls]
    if len(m) >= 5:
        Cm = tip_covariance(parent, length, [tips[s] for s in m.index])
        out["PGLS"] = fit_pgls(m["ybar"].to_numpy(), m["x"].to_numpy(), Cm, nulls)
    return out


def tidy(name, model, res, s0):
    row = {"analysis": name, "model": model, "n_species": res["n_species"],
           "n_families": res.get("n_events", np.nan), "slope": res["slope"],
           "ci_low": res["slope_ci_low"], "ci_high": res["slope_ci_high"],
           "p_vs_0": res["p_slope_vs_0"], "p_vs_1": res["p_slope_vs_1"],
           "S0": np.nan, "S0_2.5": np.nan, "S0_97.5": np.nan, "p_vs_S0": np.nan,
           "lambda": res.get("lambda"), "tree_share": res.get("tree_share_of_species_variance")}
    if s0 is not None:
        row.update({"S0": s0[0], "S0_2.5": s0[1], "S0_97.5": s0[2],
                    "p_vs_S0": res.get(f"p_slope_vs_{s0[0]:g}")})
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--loci", type=Path, required=True, help="results/loci/<cohort>")
    p.add_argument("--tree", type=Path, required=True, help="data/processed/species.tree")
    p.add_argument("--null", type=Path, help="geometry null output with empirical sizes")
    p.add_argument("--min-pgls", type=int, default=3)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--write-null-input", action="store_true",
                   help="only write the tables the null model needs, then stop")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    families = pd.read_csv(a.loci / "families.tsv", sep="\t")
    species = pd.read_csv(a.loci / "species.tsv", sep="\t")
    species = species[species["excluded"].fillna("") == ""]
    if a.write_null_input:
        null_input(families, species, a.out)
        print(f"Null input written to {a.out}")
        return

    tree = parse_newick_tree(a.tree.read_text().strip())
    s0s = read_s0(a.null)
    subset_of = {None: "all", "local": "local", "spanning": "spanning"}
    rows, data = [], {}
    for name, (tiers, pos, size_col, x_col) in ANALYSES.items():
        d = species_data(families, species, tiers, pos, size_col, x_col)
        data[name] = d
        s0 = s0s.get(subset_of[pos])
        for model, res in fit_all(d, tree, a.min_pgls, s0[0] if s0 else None).items():
            rows.append(tidy(name, model, res, s0))
    results = pd.DataFrame(rows)

    # leave one phylum out (main PMM)
    lopo = []
    d = data["main"]
    s0 = s0s.get("all")
    for phylum in sorted(d["phylum"].dropna().unique()):
        dd = d[d["phylum"] != phylum]
        if len(dd) < 10:
            continue
        res = fit_all(dd, tree, a.min_pgls, s0[0] if s0 else None)["PMM"]
        lopo.append({"left_out": phylum, "species_left_out": int((d["phylum"] == phylum).sum()),
                     **{k: v for k, v in tidy("main", "PMM", res, s0).items() if k not in ("analysis", "model")}})
    lopo = pd.DataFrame(lopo)

    results.to_csv(a.out / "model_results.tsv", sep="\t", index=False)
    lopo.to_csv(a.out / "leave_one_phylum_out.tsv", sep="\t", index=False)
    data["main"].to_csv(a.out / "species_data_main.tsv", sep="\t")
    show = results.copy()
    for c in ["slope", "ci_low", "ci_high", "S0", "lambda", "tree_share"]:
        show[c] = show[c].astype(float).round(3)
    for c in ["p_vs_0", "p_vs_1", "p_vs_S0"]:
        show[c] = show[c].map(lambda v: f"{v:.2g}" if pd.notna(v) else "")
    cols = ["analysis", "model", "n_species", "n_families", "slope", "ci_low", "ci_high",
            "S0", "p_vs_S0", "p_vs_0", "p_vs_1", "lambda", "tree_share"]
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(show[cols].to_string(index=False))
        if len(lopo):
            print("\nLeave one phylum out (main PMM):")
            print(lopo[["left_out", "species_left_out", "n_species", "slope", "ci_low", "ci_high",
                        "p_vs_S0"]].round(3).to_string(index=False))
    print(f"\nWritten to {a.out}")


if __name__ == "__main__":
    main()
