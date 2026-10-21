#!/usr/bin/env python3
"""Local inversions vs inversions spanning the origin/terminus, against genome size.

Figure 1 (position_boxplots.png): two stacked rows, same species order (by
genome size) and same size axis. Top: inversions spanning the origin or
terminus; bottom: local inversions (both breakpoints on one replication arm).
Each box = one species' high + medium events; line = the model fit.
Figure 2 (position_slopes.png): species means of both types against genome
size with both fitted lines, the comparison in one panel.

Usage (repository root, ibp-stats environment):
  python scripts/steps/analysis/plot_position_classes.py \
      --model-dir results/final_model_pmm --model PMM_phylogenetic_mixed
"""

import argparse
import re
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLOR = {"spanning": "#eb6834", "local": "#2a78d6"}          # validated categorical pair
FILL = {"spanning": "#fbd9c9", "local": "#cde2fb"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3de"
LABEL = {"spanning": "Spanning origin or terminus", "local": "Local (one replication arm)"}
ANALYSIS = {"spanning": "exploratory_spanning", "local": "exploratory_local"}


def pretty(slug):
    words = re.sub(r"_\d+_research50.*$", "", slug).split("_")
    return f"{words[0][0].upper()}. {' '.join(words[1:])}"


def fit_row(fits, cls, model):
    row = fits[(fits["analysis"] == ANALYSIS[cls]) & (fits["genome_size"] == "genome_size_bp")
               & (fits["model"] == model)]
    if row.empty:
        raise SystemExit(f"No {model} fit for {ANALYSIS[cls]} in pgls_results.tsv; rerun final_inversion_model.py")
    return row.iloc[0]


def fit_text(f):
    return (f"slope {f['slope']:.2f} [{f['slope_ci_low']:.2f}, {f['slope_ci_high']:.2f}], "
            f"p = {f['p_slope_vs_0']:.3g}")


def style(ax):
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED)


def size_axis(ax):
    ax.set_yscale("log")
    ax.set_yticks([1e4, 3e4, 1e5, 3e5, 1e6, 3e6])
    ax.set_yticklabels(["10 kb", "30 kb", "100 kb", "300 kb", "1 Mb", "3 Mb"])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model-dir", type=Path, required=True)
    p.add_argument("--model", default="PMM_phylogenetic_mixed")
    p.add_argument("--out", type=Path, help="output folder (default: the model folder)")
    args = p.parse_args(argv)
    out = args.out or args.model_dir
    out.mkdir(parents=True, exist_ok=True)

    events = pd.read_csv(args.model_dir / "events.tsv", sep="\t")
    if "position_class" not in events:
        raise SystemExit("events.tsv has no position_class column; rerun final_inversion_model.py")
    table = pd.read_csv(args.model_dir / "model_species_table.tsv", sep="\t")
    sizes = pd.read_csv(args.model_dir / "species_gtdb_mapping.tsv", sep="\t")
    fits = pd.read_csv(args.model_dir / "pgls_results.tsv", sep="\t")
    fit = {cls: fit_row(fits, cls, args.model) for cls in ("spanning", "local")}

    rows = table[(table["analysis"] == "main_high_medium") & (table["genome_size"] == "genome_size_bp")]
    in_model = {s.strip() for joined in rows["ncbi_species"] for s in joined.split(";")}
    ev = events[events["tier"].isin(["high", "medium"]) & events["species"].isin(in_model)]
    ev = ev.merge(sizes[["species", "genome_size_bp"]], on="species")
    order = (ev.drop_duplicates("species").sort_values("genome_size_bp")
             [["species", "genome_size_bp"]].reset_index(drop=True))
    pos = {s: i for i, s in enumerate(order["species"])}
    n = len(order)
    xs = np.arange(n)
    lg = np.log10(order["genome_size_bp"].to_numpy())

    # ---- figure 1: stacked box plots ---------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(max(12, 0.24 * n + 2), 10), sharex=True, sharey=True)
    rng = np.random.default_rng(0)
    for ax, cls in zip(axes, ("spanning", "local")):
        sub = ev[ev["position_class"] == cls]
        data, where = [], []
        for s, g in sub.groupby("species"):
            data.append(g["size_bp"].to_numpy())
            where.append(pos[s])
            ax.scatter(pos[s] + rng.uniform(-0.22, 0.22, len(g)), g["size_bp"], s=6,
                       color=COLOR[cls], alpha=0.45, linewidths=0, zorder=2)
        ax.boxplot(data, positions=where, widths=0.62, showfliers=False, patch_artist=True, zorder=1,
                   boxprops={"facecolor": FILL[cls], "edgecolor": COLOR[cls], "linewidth": 0.9},
                   medianprops={"color": COLOR[cls], "linewidth": 1.4},
                   whiskerprops={"color": COLOR[cls], "linewidth": 0.8},
                   capprops={"color": COLOR[cls], "linewidth": 0.8})
        f = fit[cls]
        ax.plot(xs, 10 ** (f["intercept"] + f["slope"] * lg), color=INK, linewidth=2.2, zorder=5)
        ax.set_title(f"{LABEL[cls]}: {len(sub)} inversions in {sub['species'].nunique()} species, "
                     f"median {sub['size_bp'].median() / 1e3:.0f} kb   |   {fit_text(f)}",
                     loc="left", fontsize=11, color=INK)
        size_axis(ax)
        ax.set_ylabel("Inversion size", color=INK)
        style(ax)
    axes[1].set_xticks(xs)
    axes[1].set_xticklabels([f"{pretty(s)}  {g / 1e6:.1f} Mb" for s, g in
                             zip(order["species"], order["genome_size_bp"])],
                            rotation=90, fontsize=7.5, color=INK, fontstyle="italic")
    axes[1].set_xlim(-0.8, n - 0.2)
    axes[1].set_xlabel("Species, ordered by genome size", color=INK, labelpad=8)
    fig.text(0.01, 0.985, "Large inversions spanning the origin or terminus grow with genome size; "
             "local ones do not", ha="left", va="top", fontsize=14, color=INK)
    fig.text(0.01, 0.955, f"Black line: {args.model} fit of log10 size on log10 genome size; "
             "high + medium confidence events", ha="left", va="top", fontsize=9, color=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out / "position_boxplots.png", dpi=220)

    # ---- figure 2: both lines in one panel ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5.8))
    grid = np.linspace(lg.min(), lg.max(), 50)
    handles = []
    for cls in ("spanning", "local"):
        sub = ev[ev["position_class"] == cls]
        m = sub.groupby("species").agg(y=("log10_size", "mean"), n=("log10_size", "size"),
                                       g=("genome_size_bp", "first"))
        ax.scatter(m["g"], 10 ** m["y"], s=10 + 3 * np.sqrt(m["n"]), color=COLOR[cls], alpha=0.6,
                   linewidths=0, zorder=2)
        f = fit[cls]
        ax.plot(10 ** grid, 10 ** (f["intercept"] + f["slope"] * grid), color=COLOR[cls], linewidth=2.4, zorder=3)
        handles.append(Line2D([], [], color=COLOR[cls], marker="o", linewidth=2.4,
                              label=f"{LABEL[cls]}: {fit_text(f)}"))
    ax.set_xscale("log")
    ax.set_xticks([1e6, 2e6, 4e6, 8e6])
    ax.set_xticklabels(["1 Mb", "2 Mb", "4 Mb", "8 Mb"])
    size_axis(ax)
    ax.set_xlabel("Genome size", color=INK)
    ax.set_ylabel("Species mean inversion size (geometric)", color=INK)
    ax.legend(handles=handles, frameon=False, fontsize=8.5, loc="upper left",
              bbox_to_anchor=(0, -0.14), ncol=1)
    ax.set_title("Two kinds of inversions, two relationships with genome size", loc="left",
                 fontsize=12, color=INK)
    style(ax)
    fig.tight_layout()
    fig.savefig(out / "position_slopes.png", dpi=220)
    print(f"Written {out / 'position_boxplots.png'} and {out / 'position_slopes.png'}")


if __name__ == "__main__":
    main()
