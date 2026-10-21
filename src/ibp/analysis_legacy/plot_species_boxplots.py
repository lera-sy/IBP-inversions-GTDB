#!/usr/bin/env python3
"""One box plot of inversion sizes per species, ordered by genome size, with the PGLS line.

Each box: sizes of the species' high + medium confidence inversion events (log scale).
Dots: the individual events (blue = found by both callers, orange = one caller).
Diamond: the species mean of log10(size), the value the PGLS model uses.
Red line: PGLS prediction at each species' genome size.

Usage (from the repository root, in ibp-stats):
  python scripts/steps/analysis/plot_species_boxplots.py \
      --model-dir results/final_model_all_species \
      --out results/final_model_all_species/species_boxplots.png
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
from matplotlib.patches import Patch

HIGH, MEDIUM = "#2a78d6", "#eb6834"      # validated categorical slots 1-2
BOX_FILL, BOX_EDGE = "#cde2fb", "#184f95"
LINE = "#c0392b"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3de"


def pretty(slug):
    """escherichia_coli_562_research50_m993... -> E. coli"""
    words = re.sub(r"_\d+_research50.*$", "", slug).split("_")
    return f"{words[0][0].upper()}. {' '.join(words[1:])}"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model-dir", type=Path, required=True,
                   help="Output folder of final_inversion_model.py")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--analysis", default="main_high_medium")
    p.add_argument("--model", default="PGLS_pagel_lambda")
    p.add_argument("--tiers", default="high,medium")
    args = p.parse_args(argv)

    events = pd.read_csv(args.model_dir / "events.tsv", sep="\t")
    table = pd.read_csv(args.model_dir / "model_species_table.tsv", sep="\t")
    fits = pd.read_csv(args.model_dir / "pgls_results.tsv", sep="\t")
    sizes = pd.read_csv(args.model_dir / "species_gtdb_mapping.tsv", sep="\t")
    fit = fits[(fits["analysis"] == args.analysis) & (fits["genome_size"] == "genome_size_bp")
               & (fits["model"] == args.model)].iloc[0]

    # species that entered this model (NCBI species pooled on one GTDB tip are listed together)
    rows = table[(table["analysis"] == args.analysis) & (table["genome_size"] == "genome_size_bp")]
    in_model = {s.strip() for joined in rows["ncbi_species"] for s in joined.split(";")}
    tiers = args.tiers.split(",")
    ev = events[events["tier"].isin(tiers) & events["species"].isin(in_model)]
    sp = (sizes[sizes["species"].isin(ev["species"].unique())]
          .sort_values("genome_size_bp").reset_index(drop=True))
    n = len(sp)

    fig, ax = plt.subplots(figsize=(max(12, 0.24 * n + 2), 7.5))
    fig.patch.set_facecolor("white")
    rng = np.random.default_rng(0)
    data, means = [], []
    for i, row in sp.iterrows():
        e = ev[ev["species"] == row["species"]]
        data.append(e["size_bp"].to_numpy())
        means.append(10 ** e["log10_size"].mean())
        for tier, color in (("medium", MEDIUM), ("high", HIGH)):
            vals = e.loc[e["tier"] == tier, "size_bp"].to_numpy()
            ax.scatter(i + rng.uniform(-0.22, 0.22, len(vals)), vals, s=7, color=color,
                       alpha=0.55, linewidths=0, zorder=2)

    ax.boxplot(data, positions=range(n), widths=0.62, showfliers=False, patch_artist=True,
               boxprops={"facecolor": BOX_FILL, "edgecolor": BOX_EDGE, "linewidth": 0.9},
               medianprops={"color": BOX_EDGE, "linewidth": 1.4},
               whiskerprops={"color": BOX_EDGE, "linewidth": 0.8},
               capprops={"color": BOX_EDGE, "linewidth": 0.8}, zorder=1)
    ax.scatter(range(n), means, marker="D", s=22, color=INK, zorder=4,
               edgecolors="white", linewidths=0.8)

    predicted = 10 ** (fit["intercept"] + fit["slope"] * np.log10(sp["genome_size_bp"]))
    ax.plot(range(n), predicted, color=LINE, linewidth=2.2, zorder=5)

    ax.set_yscale("log")
    ax.set_yticks([1e4, 3e4, 1e5, 3e5, 1e6, 3e6])
    ax.set_yticklabels(["10 kb", "30 kb", "100 kb", "300 kb", "1 Mb", "3 Mb"])
    ax.set_ylabel("Inversion size", color=INK, fontsize=11)
    ax.set_xticks(range(n))
    ax.set_xticklabels([f"{pretty(s)}  {g / 1e6:.1f} Mb" for s, g in zip(sp["species"], sp["genome_size_bp"])],
                       rotation=90, fontsize=7.5, color=INK, fontstyle="italic")
    ax.set_xlim(-0.8, n - 0.2)
    ax.set_xlabel("Species, ordered by genome size", color=INK, fontsize=11, labelpad=8)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED)

    lam = "" if pd.isna(fit["lambda"]) else f", λ = {fit['lambda']:.2f}"
    handles = [
        Patch(facecolor=BOX_FILL, edgecolor=BOX_EDGE, label="Inversion sizes per species"),
        Line2D([], [], marker="o", linestyle="", color=HIGH, markersize=5, label="Event found by both callers"),
        Line2D([], [], marker="o", linestyle="", color=MEDIUM, markersize=5, label="Event found by one caller"),
        Line2D([], [], marker="D", linestyle="", color=INK, markersize=5, label="Species mean (used in the model)"),
        Line2D([], [], color=LINE, linewidth=2.2,
               label=f"PGLS fit: slope {fit['slope']:.2f} "
                     f"[{fit['slope_ci_low']:.2f}, {fit['slope_ci_high']:.2f}]{lam}"),
    ]
    ax.legend(handles=handles, loc="upper left", ncol=3, frameon=False, fontsize=9,
              bbox_to_anchor=(0, 1.13))
    tips = len(rows)
    ax.set_title(f"Inversion size per species: {n} bacterial species"
                 + (f" ({tips} GTDB species in the model)" if tips != n else ""),
                 loc="left", fontsize=13, color=INK, pad=48)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=250)
    print(f"Written {args.out} ({n} species, {len(ev)} events)")


if __name__ == "__main__":
    main()
