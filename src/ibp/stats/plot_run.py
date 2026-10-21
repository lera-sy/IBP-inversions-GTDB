"""
Figures for ONE run (no team comparison). Axes are on log scales but labelled
in Mb and kb. Every species sits at its real genome size (no even spacing).

1. boxes.png          one box per species of its inversion family sizes, at the
                      species' genome size, coloured by phylum; fitted slope with
                      its 95% range (fan), geometry-only line (S0), proportional
                      line (S1).
2. scatter.png        one circle per species (mean family size, circle area =
                      number of families) + every family as a faint dot; same lines.
3. position.png       the same scatter split into spanning (crosses origin or
                      terminus) and local (one arm), each with its own S0.
4. spanning_share.png per species, the share of families that span the origin or
                      terminus vs genome size, against the geometry-only share.
5. tiers.png          families per tier and position class.
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

FIT = "#1f5fa8"
S0C = "#e08a1e"
GENOME = "median_genome_length_bp"
MB_TICKS = [1, 1.5, 2, 3, 4, 5, 7, 10]
KB_TICKS = [3, 10, 30, 100, 300, 1000, 3000]


def label_axes(ax, ylabel="Inversion family size (kb)"):
    lo, hi = ax.get_xlim()
    ticks = [t for t in MB_TICKS if lo <= np.log10(t) <= hi]
    ax.set_xticks(np.log10(ticks))
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    lo, hi = ax.get_ylim()
    ticks = [t for t in KB_TICKS if lo <= np.log10(t) <= hi]
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels([f"{t:g}" for t in ticks])
    ax.set_xlabel("Genome size (Mb, log scale)")
    ax.set_ylabel(ylabel + ", log scale")
    ax.grid(alpha=0.25, lw=0.6)


def load(loci_dir, model_dir):
    fam = pd.read_csv(loci_dir / "families.tsv", sep="\t")
    sp = pd.read_csv(loci_dir / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    res = pd.read_csv(model_dir / "model_results.tsv", sep="\t")
    null = pd.read_csv(model_dir / "null_slopes.tsv", sep="\t")
    null = null[(null["size_distribution"] == "empirical_pooled_sizes") & (null["fit"] == "events")]
    hm = fam[fam["tier"].isin(["high", "medium"])].join(
        sp.set_index("species_id")[[GENOME, "phylum"]], on="species_id", how="inner")
    hm = hm.assign(x=np.log10(hm[GENOME] / 1e6), y=np.log10(hm["size_bp"] / 1000))
    return fam, sp, hm, res, null.set_index("subset")


def colours(phyla):
    phyla = sorted(phyla)
    return dict(zip(phyla, plt.cm.tab10.colors))


def species_means(ev):
    return ev.groupby("species_id").agg(x=("x", "first"), y=("y", "mean"), n=("y", "size"),
                                        phylum=("phylum", "first"))


def draw_lines(ax, g, fit, s0):
    """Fitted slope with its 95% range as a fan, S0 and S1, all through the weighted mean."""
    xs = np.linspace(g["x"].min() - 0.02, g["x"].max() + 0.02, 60)
    xm, ym = np.average(g["x"], weights=g["n"]), np.average(g["y"], weights=g["n"])
    lo = ym + fit["ci_low"] * (xs - xm)
    hi = ym + fit["ci_high"] * (xs - xm)
    ax.fill_between(xs, np.minimum(lo, hi), np.maximum(lo, hi), color=FIT, alpha=0.12, lw=0)
    ax.plot(xs, ym + fit["slope"] * (xs - xm), color=FIT, lw=2.5,
            label=f"Fitted: slope {fit['slope']:.2f} [{fit['ci_low']:.2f}, {fit['ci_high']:.2f}]")
    if s0 is not None:
        ax.plot(xs, ym + s0["slope_mean"] * (xs - xm), color=S0C, lw=2.2, ls="--",
                label=f"Geometry only (S0): {s0['slope_mean']:.2f}")
    ax.plot(xs, ym + 1.0 * (xs - xm), color="black", lw=1, ls=":", label="Proportional (slope 1)")


def pick(res, analysis):
    return res[(res["analysis"] == analysis) & (res["model"] == "PMM")].iloc[0]


def plot_boxes(hm, res, null, out):
    cols = colours(hm["phylum"].dropna().unique())
    g = species_means(hm).sort_values("x")
    fig, ax = plt.subplots(figsize=(13, 6.3))
    data = [hm.loc[hm["species_id"] == s, "y"].to_numpy() for s in g.index]
    width = 0.012
    bp = ax.boxplot(data, positions=g["x"].to_numpy(), widths=width, patch_artist=True,
                    showfliers=False, manage_ticks=False, medianprops={"color": "black", "lw": 1.2},
                    whiskerprops={"lw": 0.8}, capprops={"lw": 0.8})
    for patch, ph in zip(bp["boxes"], g["phylum"]):
        patch.set_facecolor(cols.get(ph, "grey"))
        patch.set_alpha(0.55)
        patch.set_linewidth(0.6)
    ax.scatter(hm["x"], hm["y"], s=5, color="black", alpha=0.25, zorder=3)
    draw_lines(ax, g, pick(res, "main"), null.loc["all"] if "all" in null.index else None)
    pad = 0.03
    ax.set_xlim(g["x"].min() - pad, g["x"].max() + pad)
    label_axes(ax)
    handles, _ = ax.get_legend_handles_labels()
    handles += [Patch(facecolor=cols[p], alpha=0.55, label=f"{p} ({(g['phylum'] == p).sum()} species)")
                for p in cols]
    ax.legend(handles=handles, fontsize=8, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.13))
    ax.set_title(f"Inversion size per species ({len(g)} species, {len(hm)} high + medium families); "
                 "each box at its species' genome size", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_scatter(ax, ev, fit, s0, cols, title):
    g = species_means(ev)
    ax.scatter(ev["x"], ev["y"], s=6, color="grey", alpha=0.25, zorder=1)
    for ph, gg in g.groupby("phylum"):
        ax.scatter(gg["x"], gg["y"], s=15 + 8 * gg["n"], facecolors=cols.get(ph, "grey"),
                   edgecolors="black", linewidths=0.5, alpha=0.75, zorder=3, label=f"{ph} ({len(gg)})")
    draw_lines(ax, g, fit, s0)
    label_axes(ax)
    ax.set_title(title, fontsize=10)
    return g


def plot_main_scatter(hm, res, null, out):
    cols = colours(hm["phylum"].dropna().unique())
    fig, ax = plt.subplots(figsize=(8, 5.5))
    plot_scatter(ax, hm, pick(res, "main"), null.loc["all"] if "all" in null.index else None, cols,
                 "Mean inversion size per species (circle area = number of families)")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_position(hm, res, null, out):
    cols = colours(hm["phylum"].dropna().unique())
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, cls, analysis, name in [(axes[0], "spanning", "expl_spanning", "Spanning origin or terminus"),
                                    (axes[1], "local", "expl_local", "Local (one replication arm)")]:
        ev = hm[hm["position_class"] == cls]
        plot_scatter(ax, ev, pick(res, analysis), null.loc[cls] if cls in null.index else None, cols,
                     f"{name}: {len(ev)} families, {ev['species_id'].nunique()} species")
        ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_spanning_share(hm, null, out):
    cols = colours(hm["phylum"].dropna().unique())
    g = hm.groupby("species_id").agg(x=("x", "first"), n=("y", "size"), phylum=("phylum", "first"),
                                     share=("position_class", lambda v: (v == "spanning").mean()))
    fig, ax = plt.subplots(figsize=(9.5, 5))
    for ph, gg in g.groupby("phylum"):
        ax.scatter(gg["x"], gg["share"] * 100, s=15 + 8 * gg["n"], facecolors=cols.get(ph, "grey"),
                   edgecolors="black", linewidths=0.5, alpha=0.75, label=f"{ph} ({len(gg)})")
    overall = (hm["position_class"] == "spanning").mean()
    ax.axhline(overall * 100, color=FIT, lw=2, label=f"All families: {overall:.0%}")
    if "spanning_fraction" in null.index:
        s = null.loc["spanning_fraction"]
        ax.axhspan(s["slope_2.5"] * 100, s["slope_97.5"] * 100, color=S0C, alpha=0.25, lw=0)
        ax.axhline(s["slope_mean"] * 100, color=S0C, lw=2, ls="--",
                   label=f"Geometry only: {s['slope_mean']:.0%}")
    ax.set_ylim(-3, 103)
    lo, hi = g["x"].min() - 0.03, g["x"].max() + 0.03
    ax.set_xlim(lo, hi)
    ticks = [t for t in MB_TICKS if lo <= np.log10(t) <= hi]
    ax.set_xticks(np.log10(ticks))
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_xlabel("Genome size (Mb, log scale)")
    ax.set_ylabel("Families spanning origin or terminus (%)")
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_title("Share of spanning inversions per species (circle area = number of families)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_tiers(fam, out):
    t = (fam.groupby(["tier", "position_class"]).size().unstack(fill_value=0)
         .reindex(["high", "medium", "low"]).fillna(0))
    fig, ax = plt.subplots(figsize=(5.5, 4.4))
    left = np.zeros(len(t))
    for cls, colour in [("spanning", "#4c72b0"), ("local", "#c44e52")]:
        if cls in t:
            ax.barh(t.index, t[cls], left=left, color=colour, label=cls)
            left += t[cls].to_numpy()
    for i, total in enumerate(left):
        ax.text(total + 5, i, f"{int(total)}", va="center", fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Inversion families")
    ax.set_title("Families per confidence tier\n(main model uses high + medium)", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
    ax.set_xlim(0, left.max() * 1.15)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--loci", type=Path, required=True)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    fam, sp, hm, res, null = load(a.loci, a.model)
    plot_boxes(hm, res, null, a.out / "boxes.png")
    plot_main_scatter(hm, res, null, a.out / "scatter.png")
    plot_position(hm, res, null, a.out / "position.png")
    plot_spanning_share(hm, null, a.out / "spanning_share.png")
    plot_tiers(fam, a.out / "tiers.png")
    print(f"Figures in {a.out}: boxes.png, scatter.png, position.png, spanning_share.png, tiers.png")


if __name__ == "__main__":
    main()
