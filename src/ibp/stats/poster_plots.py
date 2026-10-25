"""
Poster figures (rerun as the main pipeline).

  question.png       "Do larger genomes allow larger inversions?" All high + medium
                     families against genome size, the most each species' chromosome
                     allows (half its length, the shorter-arc cap), the fitted slope and
                     the slope from room alone.
  split_boxes.png    One box per species at its real genome size (not evenly spaced),
                     spanning and local panels, fitted and room-only lines.
  split_scatter.png  Species mean sizes of both groups on one plot, fitted and room-only
                     lines, empirical p-values against both nulls.

Slopes come from the mixed models (split_results.tsv, model_results.tsv); null slopes
and empirical p-values from class_tests (null_tests.tsv), so every number on the
figures matches the results document.
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GENOME = "median_genome_length_bp"
COL = {"spanning": "#e8693b", "local": "#2f7ad6"}
NAME = {"spanning": "Spanning origin or terminus", "local": "Local (one replication arm)"}
MB = [1, 1.5, 2, 3, 4, 5, 7]
KB = [3, 10, 30, 100, 300, 1000, 3000]
plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12})


def ticks(ax):
    lo, hi = ax.get_xlim()
    t = [v for v in MB if lo <= np.log10(v) <= hi]
    ax.set_xticks(np.log10(t))
    ax.set_xticklabels([f"{v:g} Mb" for v in t])
    lo, hi = ax.get_ylim()
    t = [v for v in KB if lo <= np.log10(v) <= hi]
    ax.set_yticks(np.log10(t))
    ax.set_yticklabels([f"{v / 1000:g} Mb" if v >= 1000 else f"{v:g} kb" for v in t])
    ax.grid(alpha=0.25, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fmt_p(p):
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.2g}"


def null_p(p, reps=200):
    """Empirical p; at the floor (all replicates less extreme) say so."""
    return f"p ≤ {2 / (reps + 1):.2f}" if p <= 2 / (reps + 1) + 1e-9 else f"p = {p:.2f}"


def load(loci):
    sp = pd.read_csv(loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    fam = pd.read_csv(loci / "families.tsv", sep="\t")
    fam = fam[fam["tier"].isin(["high", "medium"]) & fam["species_id"].isin(sp["species_id"])
              & (fam["size_bp"] > 0)]
    fam = fam.join(sp.set_index("species_id")[[GENOME, "reference_length_bp"]], on="species_id")
    fam = fam.assign(x=np.log10(fam[GENOME] / 1e6), y=np.log10(fam["size_bp"] / 1000))
    return fam, sp


def species_means(ev):
    return ev.groupby("species_id").agg(x=("x", "first"), y=("y", "mean"), n=("y", "size"))


def line(ax, g, slope, **kw):
    xs = np.linspace(g["x"].min() - 0.02, g["x"].max() + 0.02, 50)
    xm, ym = np.average(g["x"], weights=g["n"]), np.average(g["y"], weights=g["n"])
    ax.plot(xs, ym + slope * (xs - xm), **kw)


def question(fam, sp, main, nulls, out):
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    ax.scatter(fam["x"], fam["y"], s=10, color="#9a9a9a", alpha=0.45, lw=0, label="Inversion family (828)")
    cap = sp.assign(x=np.log10(sp[GENOME] / 1e6), y=np.log10(sp["reference_length_bp"] / 2 / 1000))
    ax.scatter(cap["x"], cap["y"], marker="_", s=120, color="black", lw=1.6,
               label="Largest possible inversion (half the chromosome)")
    g = species_means(fam)
    room = nulls.loc[("shuffle", "slope_all"), "null_mean"]
    line(ax, g, main["slope"], color="#1f5fb4", lw=3,
         label=f"Observed trend: slope {main['slope']:.2f} [{main['ci_low']:.2f}, {main['ci_high']:.2f}]")
    line(ax, g, room, color="#e08a1e", lw=2.2, ls="--", label=f"Expected from extra room alone: {room:.2f}")
    ax.set_ylim(np.log10(5), np.log10(5000))
    ax.set_xlim(g["x"].min() - 0.03, g["x"].max() + 0.03)
    ticks(ax)
    ax.set_xlabel("Genome size (log scale)")
    ax.set_ylabel("Inversion size (log scale)")
    ax.set_title("Bigger genomes have more room, and inversions grow with it", loc="left")
    ax.legend(fontsize=9.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def split_boxes(fam, res, nulls, out):
    fig, axes = plt.subplots(2, 1, figsize=(13, 9.5), sharex=True, sharey=True)
    for ax, cls in zip(axes, ["spanning", "local"]):
        ev = fam[fam["position_class"] == cls]
        g = species_means(ev).sort_values("x")
        data = [ev.loc[ev["species_id"] == s, "y"].to_numpy() for s in g.index]
        bp = ax.boxplot(data, positions=g["x"].to_numpy(), widths=0.009, patch_artist=True, showfliers=False,
                        manage_ticks=False, medianprops={"color": "black", "lw": 1.3},
                        whiskerprops={"color": COL[cls], "lw": 0.9}, capprops={"color": COL[cls], "lw": 0.9})
        for b in bp["boxes"]:
            b.set_facecolor(COL[cls])
            b.set_alpha(0.4)
            b.set_edgecolor(COL[cls])
        rng = np.random.default_rng(1)
        ax.scatter(ev["x"] + rng.uniform(-0.003, 0.003, len(ev)), ev["y"], s=7, color=COL[cls], alpha=0.5,
                   lw=0, zorder=3)
        m = res[(res["group"] == cls) & (res["analysis"] == "main") & (res["model"] == "PMM")].iloc[0]
        room = nulls.loc[("shuffle", f"slope_{cls}"), "null_mean"]
        pr = nulls.loc[("shuffle", f"slope_{cls}"), "p_empirical_two_sided"]
        line(ax, g, m["slope"], color="black", lw=2.6,
             label=f"Observed trend: slope {m['slope']:.2f} [{m['ci_low']:.2f}, {m['ci_high']:.2f}], vs 0 {fmt_p(m['p_vs_0'])}")
        line(ax, g, room, color="#e08a1e", lw=2, ls="--",
             label=f"Expected from extra room alone: {room:.2f}; observed vs room {null_p(pr)}")
        ax.set_title(f"{NAME[cls]}: {len(ev)} families in {len(g)} species, median {10 ** ev['y'].median():.0f} kb",
                     loc="left")
        ax.legend(fontsize=9.5, frameon=True, framealpha=0.9, edgecolor="none",
                  loc="lower right" if cls == "spanning" else "upper left")
        ax.set_ylabel("Inversion size (log scale)")
    xs = fam["x"]
    axes[1].set_xlim(xs.min() - 0.03, xs.max() + 0.03)
    axes[1].set_ylim(np.log10(5), np.log10(5000))
    for ax in axes:
        ticks(ax)
    axes[1].set_xlabel("Genome size (log scale); each box is one species, placed at its genome size")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def split_scatter(fam, res, nulls, diff, out):
    fig, ax = plt.subplots(figsize=(8.5, 6.8))
    for cls in ["spanning", "local"]:
        g = species_means(fam[fam["position_class"] == cls])
        ax.scatter(g["x"], g["y"], s=14 + 7 * g["n"], color=COL[cls], alpha=0.6, edgecolors="white", lw=0.5)
        m = res[(res["group"] == cls) & (res["analysis"] == "main") & (res["model"] == "PMM")].iloc[0]
        room = nulls.loc[("shuffle", f"slope_{cls}"), "null_mean"]
        pr = nulls.loc[("shuffle", f"slope_{cls}"), "p_empirical_two_sided"]
        line(ax, g, m["slope"], color=COL[cls], lw=3.2,
             label=f"{NAME[cls]}: slope {m['slope']:.2f} [{m['ci_low']:.2f}, {m['ci_high']:.2f}], vs 0 {fmt_p(m['p_vs_0'])}")
        line(ax, g, room, color=COL[cls], lw=1.6, ls="--",
             label=f"    expected from room alone {room:.2f}; observed vs room {null_p(pr)}")
    ax.set_ylim(np.log10(5), np.log10(4000))
    ticks(ax)
    ax.set_xlabel("Genome size (log scale)")
    ax.set_ylabel("Species mean inversion size (log scale)")
    ax.set_title(f"Two kinds of inversions, two relationships with genome size\n"
                 f"slope difference {diff['estimate']:.2f} [{diff['ci_low']:.2f}, {diff['ci_high']:.2f}], "
                 f"{fmt_p(diff['p_vs_0'])}  (circle area = number of families)", fontsize=12)
    ax.legend(fontsize=9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True, help="results/model/<cohort>")
    ap.add_argument("--split", type=Path, required=True, help="results/figures/<cohort>/split")
    ap.add_argument("--class-tests", type=Path, required=True, help="results/class_tests/<cohort>")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    fam, sp = load(a.loci)
    mr = pd.read_csv(a.model / "model_results.tsv", sep="\t")
    main_fit = mr[(mr["analysis"] == "main") & (mr["model"] == "PMM")].iloc[0]
    res = pd.read_csv(a.split / "split_results.tsv", sep="\t")
    nulls = pd.read_csv(a.class_tests / "null_tests.tsv", sep="\t").set_index(["null", "quantity"])
    inter = pd.read_csv(a.class_tests / "interaction_test.tsv", sep="\t")
    diff = inter[inter["quantity"].str.startswith("slope difference")].iloc[0]
    question(fam, sp, main_fit, nulls, a.out / "question.png")
    split_boxes(fam, res, nulls, a.out / "split_boxes.png")
    split_scatter(fam, res, nulls, diff, a.out / "split_scatter.png")
    print(f"Written to {a.out}: question.png, split_boxes.png, split_scatter.png")


if __name__ == "__main__":
    main()
