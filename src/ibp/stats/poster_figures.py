"""
Final poster figures, drawn at the size they are placed on the poster (inches),
in the poster palette (teal / terracotta / slate, from the pipeline figure).

  fig_question.png      pooled sizes vs genome size, cap, observed vs room-only slope
  fig_split_scatter.png species means, spanning vs local, observed vs room-only
  fig_split_boxes.png   one box per species at its genome size, spanning | local
  fig_classes.png       chromosome diagram: spanning vs local inversion
  qr_repo.png           QR code to the repository

Numbers come from the same result tables as poster_plots.py.
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Arc, Circle  # noqa: E402

from ibp.stats.poster_plots import load, species_means  # noqa: E402

# Okabe-Ito colour-blind safe pair: vermillion (spanning) and blue (local)
VERM, BLUE, NAVY, GREY, LIGHT = "#D55E00", "#0072B2", "#1F3A5F", "#7D8C99", "#C5CFD8"
TEAL, TERRA, SLATE = BLUE, VERM, NAVY
COL = {"spanning": VERM, "local": BLUE}
MB = [1, 1.5, 2, 3, 5, 7]
KB = [10, 30, 100, 300, 1000, 3000]
plt.rcParams.update({"font.family": "Carlito", "font.size": 20, "axes.edgecolor": SLATE,
                     "axes.labelcolor": SLATE, "xtick.color": SLATE, "ytick.color": SLATE,
                     "axes.linewidth": 1.4, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
                     "xtick.labelsize": 18, "ytick.labelsize": 18})


def axes_style(ax, ylabel="Inversion size", xlabel="Genome size", ylabel_pos=1.03):
    lo, hi = ax.get_xlim()
    t = [v for v in MB if lo <= np.log10(v) <= hi]
    ax.set_xticks(np.log10(t))
    ax.set_xticklabels([f"{v:g} Mb" for v in t])
    lo, hi = ax.get_ylim()
    t = [v for v in KB if lo <= np.log10(v) <= hi]
    ax.set_yticks(np.log10(t))
    ax.set_yticklabels([f"{v / 1000:g} Mb" if v >= 1000 else f"{v:g} kb" for v in t])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#E3E9EB", lw=1)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel if "log" in xlabel else xlabel + " (log scale)", fontsize=19)
    if ylabel:  # horizontal, above the axis
        ax.text(0, ylabel_pos, ylabel + " (log scale)", transform=ax.transAxes, fontsize=19, color=SLATE,
                ha="left", va="bottom")


def fit_line(g, slope, x0, x1):
    xs = np.linspace(x0, x1, 50)
    xm, ym = np.average(g["x"], weights=g["n"]), np.average(g["y"], weights=g["n"])
    return xs, ym + slope * (xs - xm)


def end_label(ax, xs, ys, text, color, dy=0.0, weight="bold", size=19):
    ax.text(xs[-1] + 0.012, ys[-1] + dy, text, color=color, fontsize=size, fontweight=weight, va="center",
            ha="left")


def fig_question(fam, sp, main, room, out, size=(11.8, 7.6)):
    fig, ax = plt.subplots(figsize=size)
    ax.scatter(fam["x"], fam["y"], s=26, color=LIGHT, lw=0, zorder=2)
    cap = sp.assign(x=np.log10(sp["median_genome_length_bp"] / 1e6),
                    y=np.log10(sp["reference_length_bp"] / 2 / 1000),
                    ratio=sp["reference_length_bp"] / sp["median_genome_length_bp"])
    ax.scatter(cap["x"], cap["y"], marker="_", s=260, color=SLATE, lw=2.4, zorder=3)
    for _, r in cap[cap["ratio"] < 0.8].iterrows():
        ax.text(r["x"] + 0.012, r["y"] - 0.02, "*", ha="left", va="center", fontsize=24, color=SLATE, fontweight="bold")
    g = species_means(fam)
    x0, x1 = g["x"].min() - 0.01, g["x"].max() + 0.01
    xs, ys = fit_line(g, main["slope"], x0, x1)
    ax.plot(xs, ys, color=NAVY, lw=5, zorder=4, solid_capstyle="round")
    xs2, ys2 = fit_line(g, room, x0, x1)
    ax.plot(xs2, ys2, color=GREY, lw=3, ls=(0, (5, 3)), zorder=4)
    end_label(ax, xs, ys, f"Observed {main['slope']:.2f}", NAVY, dy=0.13)
    end_label(ax, xs2, ys2, f"Room alone {room:.2f}", GREY, dy=-0.13, weight="normal")
    ax.text(np.log10(1.0), np.log10(4200), "\u2014 Largest possible inversion\n    (half the chromosome; * two chromosomes)",
            color=SLATE, fontsize=15, va="top")
    ax.set_ylim(np.log10(5), np.log10(4500))
    ax.set_xlim(g["x"].min() - 0.03, g["x"].max() + 0.27)
    axes_style(ax)
    fig.tight_layout()
    fig.savefig(out, dpi=220, transparent=True)
    plt.close(fig)


def fig_split_scatter(fam, res, rooms, out, size=(10.4, 7.6)):
    fig, ax = plt.subplots(figsize=size)
    for cls in ["spanning", "local"]:
        g = species_means(fam[fam["position_class"] == cls])
        ax.scatter(g["x"], g["y"], s=40 + 14 * g["n"], color=COL[cls], alpha=0.8, edgecolors="white", lw=1,
                   zorder=3)
        m = res[(res["group"] == cls) & (res["analysis"] == "main") & (res["model"] == "PMM")].iloc[0]
        x0, x1 = g["x"].min() - 0.01, g["x"].max() + 0.01
        xs, ys = fit_line(g, m["slope"], x0, x1)
        xs2, ys2 = fit_line(g, rooms[cls], x0, x1)
        ax.plot(xs, ys, color=COL[cls], lw=5.5, zorder=4, solid_capstyle="round")
        ax.plot(xs2, ys2, color=COL[cls], lw=2.4, ls=(0, (5, 3)), zorder=4)
        name = "Spanning" if cls == "spanning" else "Local"
        up = ys[-1] >= ys2[-1]
        end_label(ax, xs, ys, f"{name} {m['slope']:+.2f}", COL[cls], dy=0.11 if up else -0.11)
        end_label(ax, xs2, ys2, f"room {rooms[cls]:+.2f}", COL[cls], dy=-0.11 if up else 0.11,
                  weight="normal", size=17)
    ax.text(0.995, 0.015, "circle area = number of families", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=14, color=GREY)
    ax.set_ylim(np.log10(5), np.log10(4000))
    ax.set_xlim(fam["x"].min() - 0.03, fam["x"].max() + 0.30)
    axes_style(ax, ylabel="Mean inversion size per species")
    fig.tight_layout()
    fig.savefig(out, dpi=220, transparent=True)
    plt.close(fig)


def fig_split_boxes(fam, res, rooms, out, size=(10.4, 7.6)):
    fig, axes = plt.subplots(2, 1, figsize=size, sharex=True, sharey=True)
    for ax, cls in zip(axes, ["spanning", "local"]):
        ev = fam[fam["position_class"] == cls]
        g = species_means(ev).sort_values("x")
        data = [ev.loc[ev["species_id"] == s, "y"].to_numpy() for s in g.index]
        bp = ax.boxplot(data, positions=g["x"].to_numpy(), widths=0.011, patch_artist=True, showfliers=False,
                        manage_ticks=False, medianprops={"color": NAVY, "lw": 1.4},
                        whiskerprops={"color": COL[cls], "lw": 1.0}, capprops={"color": COL[cls], "lw": 1.0})
        for b in bp["boxes"]:
            b.set_facecolor(COL[cls])
            b.set_alpha(0.45)
            b.set_edgecolor(COL[cls])
        rng = np.random.default_rng(1)
        ax.scatter(ev["x"] + rng.uniform(-0.003, 0.003, len(ev)), ev["y"], s=8, color=COL[cls], alpha=0.6, lw=0,
                   zorder=3)
        m = res[(res["group"] == cls) & (res["analysis"] == "main") & (res["model"] == "PMM")].iloc[0]
        x0, x1 = fam["x"].min() - 0.01, fam["x"].max() + 0.01
        xs, ys = fit_line(g, m["slope"], x0, x1)
        xs2, ys2 = fit_line(g, rooms[cls], x0, x1)
        ax.plot(xs, ys, color=NAVY, lw=3.5, zorder=5, solid_capstyle="round")
        ax.plot(xs2, ys2, color=NAVY, lw=2, ls=(0, (5, 3)), zorder=5)
        up = ys[-1] >= ys2[-1]
        end_label(ax, xs, ys, f"observed {m['slope']:+.2f}", NAVY, dy=0.2 if up else -0.2, size=16)
        end_label(ax, xs2, ys2, f"room {rooms[cls]:+.2f}", NAVY, dy=-0.2 if up else 0.2, weight="normal",
                  size=16)
        name = "Spanning" if cls == "spanning" else "Local"
        ax.text(0.01, 0.97, f"{name}: {len(ev)} families, median {10 ** ev['y'].median():.0f} kb",
                transform=ax.transAxes, ha="left", va="top", fontsize=17, color=COL[cls], fontweight="bold")
        ax.set_xlim(fam["x"].min() - 0.03, fam["x"].max() + 0.30)
    axes[1].set_ylim(np.log10(5), np.log10(8000))
    axes_style(axes[0], xlabel=None, ylabel="Inversion size", ylabel_pos=1.02)
    axes_style(axes[1], ylabel=None, xlabel="Genome size, log scale (box = one species)")
    for ax in axes:
        ax.tick_params(labelsize=16)
    fig.tight_layout(h_pad=0.6)
    fig.savefig(out, dpi=220, transparent=True)
    plt.close(fig)


def fig_classes(out, size=(5.6, 4.6)):
    fig, ax = plt.subplots(figsize=size)
    ax.set_aspect("equal")
    ax.axis("off")
    R = 1.0
    ax.add_patch(Circle((0, 0), R, fill=False, color=GREY, lw=7))
    # ori at top (90 deg), ter at bottom (270 deg)
    for ang, lab, dy in [(90, "ori", -0.22), (270, "ter", 0.22)]:
        x, y = np.cos(np.radians(ang)) * R, np.sin(np.radians(ang)) * R
        ax.plot([x], [y], "o", color=SLATE, ms=11, zorder=5)
        ax.text(x, y + dy, lab, ha="center", va="center", fontsize=20, color=SLATE, fontweight="bold")
    ax.plot([0, 0], [-R + 0.38, R - 0.38], color=GREY, lw=1.4, ls=(0, (4, 3)))
    # spanning: symmetric around ori (60..120 deg), outside the circle
    ax.add_patch(Arc((0, 0), 2.36, 2.36, theta1=55, theta2=125, color=TERRA, lw=9, capstyle="butt"))
    ax.text(-0.95, 1.12, "spanning", ha="right", va="center", fontsize=20, color=TERRA, fontweight="bold")
    # local: on the right arm only (-25..15 deg)
    ax.add_patch(Arc((0, 0), 2.36, 2.36, theta1=-28, theta2=8, color=TEAL, lw=9, capstyle="butt"))
    ax.text(1.32, -0.12, "local", ha="left", va="center", fontsize=20, color=TEAL, fontweight="bold")
    ax.text(-0.5, 0, "arm 1", ha="center", va="center", fontsize=16, color=GREY)
    ax.text(0.5, 0, "arm 2", ha="center", va="center", fontsize=16, color=GREY)
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-1.45, 1.45)
    fig.tight_layout(pad=0.1)
    fig.savefig(out, dpi=220, transparent=True)
    plt.close(fig)


def qr(url, out):
    import qrcode
    q = qrcode.QRCode(border=1, box_size=20, error_correction=qrcode.constants.ERROR_CORRECT_M)
    q.add_data(url)
    q.make(fit=True)
    q.make_image(fill_color=SLATE, back_color="white").save(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--split", type=Path, required=True)
    ap.add_argument("--class-tests", type=Path, required=True)
    ap.add_argument("--repo-url", default="https://github.com/lera-sy/IBP-inversions-GTDB")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    fam, sp = load(a.loci)
    mr = pd.read_csv(a.model / "model_results.tsv", sep="\t")
    main_fit = mr[(mr["analysis"] == "main") & (mr["model"] == "PMM")].iloc[0]
    res = pd.read_csv(a.split / "split_results.tsv", sep="\t")
    nulls = pd.read_csv(a.class_tests / "null_tests.tsv", sep="\t").set_index(["null", "quantity"])
    rooms = {c: nulls.loc[("shuffle", f"slope_{c}"), "null_mean"] for c in ["spanning", "local"]}
    fig_question(fam, sp, main_fit, nulls.loc[("shuffle", "slope_all"), "null_mean"], a.out / "fig_question.png")
    fig_split_scatter(fam, res, rooms, a.out / "fig_split_scatter.png")
    fig_split_boxes(fam, res, rooms, a.out / "fig_split_boxes.png")
    fig_classes(a.out / "fig_classes.png")
    qr(a.repo_url, a.out / "qr_repo.png")
    print(f"Written to {a.out}")


if __name__ == "__main__":
    main()
