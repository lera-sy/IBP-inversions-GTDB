"""
The research question asked separately for the two kinds of inversions:
  spanning  breakpoints on opposite replication arms (crosses origin or terminus)
  local     both breakpoints on one replication arm

For each group, the full analysis of the main question:
  - mixed model with tree (main), mixed model without tree, PGLS on species means
  - each slope tested against 0, 1 and the group's own geometry-only slope (S0)
  - sensitivity: high tier only, all tiers, chromosome length instead of genome size
  - leave one phylum out (main model)
Figures:
  split_boxes.png    two panels, one box per species at its real genome size
  split_scatter.png  species geometric means of both groups with their fitted
                     lines and their geometry-only lines
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from ibp.analysis_legacy.final_inversion_model import parse_newick_tree
from ibp.stats.model import fit_all, read_s0, species_data

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GENOME, CHROM = "median_genome_length_bp", "median_chromosome_length_bp"
GROUPS = {"spanning": ("Spanning origin or terminus", "#e8693b"),
          "local": ("Local (one replication arm)", "#2f7ad6")}
MB_TICKS = [1, 1.5, 2, 3, 4, 5, 7, 10]
KB_TICKS = [3, 10, 30, 100, 300, 1000, 3000]


def fmt_p(p):
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.2g}"


def row(group, analysis, model, r, s0):
    return {"group": group, "analysis": analysis, "model": model, "n_species": r["n_species"],
            "n_families": r.get("n_events", np.nan), "slope": r["slope"],
            "ci_low": r["slope_ci_low"], "ci_high": r["slope_ci_high"],
            "p_vs_0": r["p_slope_vs_0"], "p_vs_1": r["p_slope_vs_1"],
            "S0": s0, "p_vs_S0": r.get(f"p_slope_vs_{s0:g}", np.nan)}


def analyse(fam, sp, tree, s0s):
    rows, lopo = [], []
    for group in GROUPS:
        s0 = s0s.get(group, (None,))[0]
        for analysis, tiers, xcol in [("main", ["high", "medium"], GENOME),
                                      ("high_only", ["high"], GENOME),
                                      ("all_tiers", ["high", "medium", "low"], GENOME),
                                      ("chromosome_length", ["high", "medium"], CHROM)]:
            d = species_data(fam, sp, tiers, group, "size_bp", xcol)
            for model, r in fit_all(d, tree, 3, s0).items():
                rows.append(row(group, analysis, model, r, s0))
        d = species_data(fam, sp, ["high", "medium"], group, "size_bp", GENOME)
        for ph in sorted(d["phylum"].dropna().unique()):
            dd = d[d["phylum"] != ph]
            if len(dd) >= 10:
                r = fit_all(dd, tree, 3, s0)["PMM"]
                lopo.append({"group": group, "left_out": ph,
                             "species_left_out": int((d["phylum"] == ph).sum()),
                             **{k: v for k, v in row(group, "main", "PMM", r, s0).items()
                                if k not in ("group", "analysis", "model")}})
    return pd.DataFrame(rows), pd.DataFrame(lopo)


def ticks(ax):
    lo, hi = ax.get_xlim()
    t = [v for v in MB_TICKS if lo <= np.log10(v) <= hi]
    ax.set_xticks(np.log10(t))
    ax.set_xticklabels([f"{v:g} Mb" for v in t])
    lo, hi = ax.get_ylim()
    t = [v for v in KB_TICKS if lo <= np.log10(v) <= hi]
    ax.set_yticks(np.log10(t))
    ax.set_yticklabels([f"{v / 1000:g} Mb" if v >= 1000 else f"{v:g} kb" for v in t])
    ax.grid(alpha=0.25, lw=0.6)


def events(fam, sp, group):
    ev = fam[fam["tier"].isin(["high", "medium"]) & (fam["position_class"] == group)]
    ev = ev.join(sp.set_index("species_id")[[GENOME, "phylum"]], on="species_id", how="inner")
    return ev.assign(x=np.log10(ev[GENOME] / 1e6), y=np.log10(ev["size_bp"] / 1000))


def line(ax, g, slope, colour, ls="-", lw=2.5, label=None):
    xs = np.linspace(g["x"].min() - 0.02, g["x"].max() + 0.02, 50)
    xm, ym = np.average(g["x"], weights=g["n"]), np.average(g["y"], weights=g["n"])
    ax.plot(xs, ym + slope * (xs - xm), color=colour, ls=ls, lw=lw, label=label)


def plot_boxes(fam, sp, res, s0s, out):
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True, sharey=True)
    for ax, (group, (name, colour)) in zip(axes, GROUPS.items()):
        ev = events(fam, sp, group)
        g = ev.groupby("species_id").agg(x=("x", "first"), y=("y", "mean"), n=("y", "size")).sort_values("x")
        data = [ev.loc[ev["species_id"] == s, "y"].to_numpy() for s in g.index]
        bp = ax.boxplot(data, positions=g["x"].to_numpy(), widths=0.012, patch_artist=True, showfliers=False,
                        manage_ticks=False, medianprops={"color": "black", "lw": 1.1},
                        whiskerprops={"color": colour, "lw": 0.8}, capprops={"color": colour, "lw": 0.8})
        for b in bp["boxes"]:
            b.set_facecolor(colour)
            b.set_alpha(0.35)
            b.set_edgecolor(colour)
        ax.scatter(ev["x"], ev["y"], s=5, color=colour, alpha=0.35, zorder=3)
        m = res[(res["group"] == group) & (res["analysis"] == "main") & (res["model"] == "PMM")].iloc[0]
        line(ax, g, m["slope"], "black", label=f"Fitted slope {m['slope']:.2f} [{m['ci_low']:.2f}, {m['ci_high']:.2f}]")
        s0 = s0s.get(group, (None,))[0]
        if s0 is not None:
            line(ax, g, s0, "#e08a1e", ls="--", lw=2, label=f"Geometry only (S0) {s0:.2f}")
        ax.set_title(f"{name}: {len(ev)} families in {len(g)} species, median {10 ** ev['y'].median():.0f} kb  |  "
                     f"vs 0: {fmt_p(m['p_vs_0'])}, vs S0: {fmt_p(m['p_vs_S0'])}, vs 1: {fmt_p(m['p_vs_1'])}",
                     fontsize=10, loc="left")
        ax.legend(fontsize=8, frameon=False, loc="lower right")
        ax.set_ylabel("Inversion size (log scale)")
    xs = events(fam, sp, "spanning")["x"]
    axes[1].set_xlim(xs.min() - 0.03, xs.max() + 0.03)
    for ax in axes:
        ticks(ax)
    axes[1].set_xlabel("Genome size (log scale); each box is one species at its genome size")
    fig.suptitle("Inversion size vs genome size, asked separately for spanning and local inversions "
                 "(high + medium families)", fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_scatter(fam, sp, res, s0s, out):
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    for group, (name, colour) in GROUPS.items():
        ev = events(fam, sp, group)
        g = ev.groupby("species_id").agg(x=("x", "first"), y=("y", "mean"), n=("y", "size"))
        ax.scatter(g["x"], g["y"], s=12 + 6 * g["n"], color=colour, alpha=0.6, edgecolors="white", linewidths=0.5)
        m = res[(res["group"] == group) & (res["analysis"] == "main") & (res["model"] == "PMM")].iloc[0]
        line(ax, g, m["slope"], colour, lw=3,
             label=f"{name}: slope {m['slope']:.2f} [{m['ci_low']:.2f}, {m['ci_high']:.2f}], vs 0 {fmt_p(m['p_vs_0'])}")
        s0 = s0s.get(group, (None,))[0]
        if s0 is not None:
            line(ax, g, s0, colour, ls="--", lw=1.6,
                 label=f"   geometry only {s0:.2f}; observed vs geometry {fmt_p(m['p_vs_S0'])}")
    ticks(ax)
    ax.set_xlabel("Genome size (log scale)")
    ax.set_ylabel("Species mean inversion size (geometric, log scale)")
    ax.set_title("Two kinds of inversions, two relationships with genome size\n"
                 "(circle area = number of families; dashed = what geometry alone produces)", fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.11))
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--null", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    fam = pd.read_csv(a.loci / "families.tsv", sep="\t")
    sp = pd.read_csv(a.loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    tree = parse_newick_tree(a.tree.read_text().strip())
    s0s = read_s0(a.null)

    res, lopo = analyse(fam, sp, tree, s0s)
    res.to_csv(a.out / "split_results.tsv", sep="\t", index=False)
    lopo.to_csv(a.out / "split_leave_one_phylum_out.tsv", sep="\t", index=False)
    plot_boxes(fam, sp, res, s0s, a.out / "split_boxes.png")
    plot_scatter(fam, sp, res, s0s, a.out / "split_scatter.png")

    show = res.copy()
    for c in ["p_vs_0", "p_vs_1", "p_vs_S0"]:
        show[c] = show[c].map(lambda v: f"{v:.2g}" if pd.notna(v) else "")
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        for group in GROUPS:
            print(f"\n===== {group.upper()} =====")
            print(show[show["group"] == group].drop(columns="group").round(3).to_string(index=False))
            lp = lopo[lopo["group"] == group]
            if len(lp):
                print("Leave one phylum out (mixed model):")
                print(lp[["left_out", "species_left_out", "n_species", "slope", "ci_low", "ci_high",
                          "p_vs_0", "p_vs_S0"]].round(4).to_string(index=False))
    print(f"\nFigures: {a.out / 'split_boxes.png'}, {a.out / 'split_scatter.png'}")


if __name__ == "__main__":
    main()
