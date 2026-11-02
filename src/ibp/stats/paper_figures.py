"""
Paper figures and the numbers behind them, for both datasets.

  core      species with >= 50 genomes (50 per species)       e.g. selected_50_v1
  expanded  species with >= 20 genomes (up to 50 per species) e.g. selected_min20_v1

Figures are drawn at Bioinformatics single-column width (8.5 cm), Arial-metric font:
  fig1_design.png      position classes on the chromosome + workflow with counts
  fig2_core.png        core dataset: pooled sizes and the spanning/local split
  figS1_expanded.png   the same for the expanded dataset
  fig3_slopes.png      forest plot of slopes, both datasets, with room-alone nulls
  fig4_position.png    symmetry around ori/ter and spanning share vs genome size

Also written: paper_numbers.tsv (extra fits computed here: species not in the core
set, core species within the expanded run, spanning share vs genome size).

Inputs are the result folders written by the pipeline (results/<step>/<cohort>) and
the species table and tree of each cohort.
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import font_manager
from scipy import optimize, stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Arc, Circle, FancyBboxPatch  # noqa: E402

from ibp.analysis_legacy.final_inversion_model import gls, lambda_matrix  # noqa: E402
from ibp.stats.class_tests import interaction_test  # noqa: E402
from ibp.stats.robustness import GENOME, HM, Tree, pmm, summarise  # noqa: E402

VERM, BLUE, INK, GREY, LIGHT, BAND = "#D55E00", "#0072B2", "#1F2A36", "#6B7785", "#C9D1D9", "#E4E8EC"
COL = {"spanning": VERM, "local": BLUE, "all": INK}
W = 3.35  # single column, inches
MB = [1, 2, 3, 5, 8]
KB = [10, 30, 100, 300, 1000, 3000]
MARK = {"core": dict(marker="o", ms=4.2, mfc="fill", label="Core (≥50 genomes)"),
        "expanded": dict(marker="D", ms=3.6, mfc="white", label="Expanded (≥20 genomes)")}


def set_style():
    names = {f.name for f in font_manager.fontManager.ttflist}
    family = next((f for f in ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"] if f in names), "sans-serif")
    plt.rcParams.update({"font.family": family, "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                         "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                         "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                         "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                         "legend.fontsize": 6.5, "legend.frameon": False, "savefig.dpi": 600})


# ------------------------------------------------------------------ data

class Data:
    def __init__(self, root, cohort, species_table, tree):
        r = Path(root) / "results"
        self.loci = r / "loci" / cohort
        sp = pd.read_csv(self.loci / "species.tsv", sep="\t")
        self.sp_all = sp
        self.sp = sp[sp["excluded"].fillna("") == ""]
        fam = pd.read_csv(self.loci / "families.tsv", sep="\t")
        fam = fam[fam["tier"].isin(HM) & fam["species_id"].isin(self.sp["species_id"]) & (fam["size_bp"] > 0)]
        fam = fam.join(self.sp.set_index("species_id")[[GENOME, "reference_length_bp"]], on="species_id")
        self.fam = fam.assign(x=np.log10(fam[GENOME] / 1e6), y=np.log10(fam["size_bp"] / 1000))
        self.tax = pd.read_csv(species_table, sep="\t")
        self.tree_path = tree
        self.tree = Tree(tree, self.sp["species_id"].tolist())
        rd = lambda step, f: pd.read_csv(r / step / cohort / f, sep="\t")  # noqa: E731
        self.model = rd("model", "model_results.tsv")
        self.split = pd.read_csv(r / "figures" / cohort / "split" / "split_results.tsv", sep="\t")
        self.nulls = rd("class_tests", "null_tests.tsv").set_index(["null", "quantity"])
        self.inter = rd("class_tests", "interaction_test.tsv")
        self.symmetry = rd("class_tests", "symmetry.tsv")
        self.symsplit = rd("class_tests", "symmetric_split.tsv")
        self.share = rd("class_tests", "spanning_share.tsv")
        self.rob = rd("robustness", "robustness.tsv")
        extra = r / "extra" / cohort
        if not (extra / "ori_ter.tsv").is_file():
            extra = r / "ori_ter" / cohort
        self.oriter = pd.read_csv(extra / "ori_ter.tsv", sep="\t") if (extra / "ori_ter.tsv").is_file() else None

    def rob_row(self, check, cls, setting=None):
        d = self.rob[(self.rob["check"] == check) & (self.rob["class"] == cls)]
        if setting is not None:
            d = d[np.isclose(d["setting"].astype(float), setting)]
        return d.iloc[0]

    def split_row(self, cls, analysis):
        d = self.split
        return d[(d["group"] == cls) & (d["analysis"] == analysis) & (d["model"] == "PMM")].iloc[0]

    def room(self, q):
        r = self.nulls.loc[("shuffle", q)]
        return r["null_mean"], r["null_2.5"], r["null_97.5"]


def est(slope, lo, hi):
    return float(slope), float(lo), float(hi)


def pgls_slope(y, x, C):
    X = np.column_stack([np.ones(len(y)), x])
    lam = float(optimize.minimize_scalar(lambda v: -gls(y, X, lambda_matrix(C, v))[3], bounds=(0, 1),
                                         method="bounded").x)
    for edge in (0.0, 1.0):
        if gls(y, X, lambda_matrix(C, edge))[3] > gls(y, X, lambda_matrix(C, lam))[3]:
            lam = edge
    b, cov, rss, _ = gls(y, X, lambda_matrix(C, lam))
    n = len(y)
    se = float(np.sqrt(cov[1, 1] * rss / (n - 2)))
    t = stats.t.ppf(0.975, n - 2)
    return {"slope": b[1], "ci_low": b[1] - t * se, "ci_high": b[1] + t * se,
            "p": 2 * stats.t.sf(abs(b[1] / se), n - 2), "lambda": lam, "n_species": n, "intercept": b[0]}


def share_by_species(d, min_families=3):
    g = d.fam.groupby("species_id").agg(n=("y", "size"), share=("position_class", lambda s: (s == "spanning").mean()),
                                        x=("x", "first"))
    return g[(g["n"] >= min_families) & g.index.isin(d.tree.pos)]


def extra_numbers(core, exp):
    rows = []
    for name, d in [("core", core), ("expanded", exp)]:
        g = share_by_species(d)
        r = pgls_slope(g["share"].to_numpy(), g["x"].to_numpy(), d.tree.sub(g.index))
        rows.append({"dataset": name, "analysis": "spanning share vs log10 genome size (PGLS, species >= 3 families)",
                     **{k: r[k] for k in ["slope", "ci_low", "ci_high", "p", "lambda", "n_species"]}})
    core_ids = set(core.tax["species_id"])
    for label, keep in [("species not in the core set", lambda s: s not in core_ids),
                        ("core species, expanded run", lambda s: s in core_ids)]:
        sp = exp.sp[exp.sp["species_id"].map(keep)]
        tree = Tree(exp.tree_path, sp["species_id"].tolist())
        fam = exp.fam[exp.fam["species_id"].isin(sp["species_id"])]
        for cls in ["all", "spanning", "local"]:
            ev = fam if cls == "all" else fam[fam["position_class"] == cls]
            r = pmm(summarise(ev, sp), tree)
            rows.append({"dataset": "expanded", "analysis": f"{label}: {cls}", "slope": r["slope"],
                         "ci_low": r["ci_low"], "ci_high": r["ci_high"], "p": r["p_vs_0"],
                         "n_species": r["n_species"], "n_families": r["n_families"]})
        it, meta, _, _ = interaction_test(fam.assign(**{"class": fam["position_class"]}), sp, tree)
        dif = it[it["quantity"] == "slope difference spanning - local"].iloc[0]
        rows.append({"dataset": "expanded", "analysis": f"{label}: spanning - local (joint model)",
                     "slope": dif["estimate"], "ci_low": dif["ci_low"], "ci_high": dif["ci_high"],
                     "p": dif["p_vs_0"], "p_LRT": dif["p_LRT"], "n_species": meta["n_species"]})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ helpers

def log_axes(ax, xlabel="Genome size (Mb, log scale)", ylabel="Inversion size (log scale)", ysize=True):
    lo, hi = ax.get_xlim()
    t = [v for v in MB if lo <= np.log10(v) <= hi]
    ax.set_xticks(np.log10(t))
    ax.set_xticklabels([f"{v:g}" for v in t])
    if ysize:
        lo, hi = ax.get_ylim()
        t = [v for v in KB if lo <= np.log10(v) <= hi]
        ax.set_yticks(np.log10(t))
        ax.set_yticklabels([f"{v / 1000:g} Mb" if v >= 1000 else f"{v:g} kb" for v in t])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#EEF1F4", lw=0.5)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def species_means(ev):
    return ev.groupby("species_id").agg(x=("x", "first"), y=("y", "mean"), n=("y", "size"))


def fit_line(g, slope, x0, x1):
    xs = np.linspace(x0, x1, 50)
    xm, ym = np.average(g["x"], weights=g["n"]), np.average(g["y"], weights=g["n"])
    return xs, ym + slope * (xs - xm)


def panel_label(ax, s, x=-0.16, y=1.03):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="bottom", ha="left", color=INK)


# ------------------------------------------------------------------ figure 1

def fig1(core, exp, out):
    fig = plt.figure(figsize=(W, 4.45))
    ax = fig.add_axes([0.0, 0.60, 1.0, 0.40])
    ax.set_aspect("equal")
    ax.axis("off")
    R = 1.0
    ax.add_patch(Circle((0, 0), R, fill=False, color=LIGHT, lw=3.2))
    for ang, lab, dy in [(90, "ori", -0.2), (270, "ter", 0.2)]:
        x, y = np.cos(np.radians(ang)) * R, np.sin(np.radians(ang)) * R
        ax.plot([x], [y], "o", color=INK, ms=4, zorder=5)
        ax.text(x, y + dy, lab, ha="center", va="center", fontsize=7, color=INK, fontweight="bold")
    ax.plot([0, 0], [-R + 0.32, R - 0.32], color=GREY, lw=0.6, ls=(0, (3, 2)))
    ax.text(-0.45, 0.0, "replication\narm 1", ha="center", va="center", fontsize=6, color=GREY)
    ax.text(0.45, 0.0, "replication\narm 2", ha="center", va="center", fontsize=6, color=GREY)
    r2 = 2.3
    # symmetric spanning, around ori
    ax.add_patch(Arc((0, 0), r2, r2, theta1=58, theta2=122, color=VERM, lw=3.4, capstyle="butt"))
    ax.text(0, 1.38, "spanning, symmetric (d1 ≈ d2)", ha="center", va="bottom", fontsize=6.3, color=VERM,
            fontweight="bold")
    for a, lab in [(74, "d1"), (106, "d2")]:
        ax.text(np.cos(np.radians(a)) * 0.8, np.sin(np.radians(a)) * 0.8, lab, fontsize=5.8, color=VERM,
                ha="center", va="center")
    # asymmetric spanning, around ter
    ax.add_patch(Arc((0, 0), r2, r2, theta1=232, theta2=282, color=VERM, lw=3.4, capstyle="butt", alpha=0.55))
    ax.text(0.05, -1.4, "spanning, asymmetric (d1 ≠ d2)", ha="center", va="top", fontsize=6.3, color=VERM)
    # local
    ax.add_patch(Arc((0, 0), r2, r2, theta1=-22, theta2=10, color=BLUE, lw=3.4, capstyle="butt"))
    ax.text(1.24, -0.1, "local\n(one arm)", ha="left", va="center", fontsize=6.3, color=BLUE, fontweight="bold")
    fig.text(0.02, 0.86, "Size = shorter arc\nbetween the two\nbreakpoints (at most\nhalf the chromosome)",
             ha="left", va="top", fontsize=5.8, color=GREY)
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-1.62, 1.6)
    fig.text(0.02, 0.995, "a", fontsize=8.5, fontweight="bold", va="top")

    bx = fig.add_axes([0.0, 0.0, 1.0, 0.58])
    bx.axis("off")
    bx.set_xlim(0, 1)
    bx.set_ylim(0, 1)
    bx.text(0.02, 0.995, "b", fontsize=8.5, fontweight="bold", va="top")
    n = lambda v: f"{int(v):,}"  # noqa: E731
    cs, es = core.sp, exp.sp
    steps = [
        ("GTDB r232 genomes", "RefSeq, complete, not MAG; CheckM2 ≥95% / ≤5%;\nlong-read de novo; one per BioSample",
         "37,451 genomes"),
        ("Species and sampling", "core: species ≥50 genomes · expanded: ≥20;\nat most 50 genomes per species (seed 42)",
         f"{len(core.sp_all)} / {len(exp.sp_all)} species\n{n(core.tax['n_selected'].sum())} / "
         f"{n(exp.tax['n_selected'].sum())} genomes"),
        ("Chromosome preparation", "largest chromosome rotated to dnaA (dnaapler);\nreference = Mash medoid, circular only",
         f"{len(cs)} / {len(es)} species kept"),
        ("Two inversion callers", "synteny-block caller (2 kb windows, minimap2)\nand minimap2 + SyRI; shorter arc, ≥6 kb",
         f"{n(cs['calls_kept'].sum())} / {n(es['calls_kept'].sum())} calls"),
        ("Loci, families, tiers", "overlapping calls → loci; same breakpoints\n(±5 kb) → one family per species; 3 evidence lines",
         f"{n(len(core.fam))} / {n(len(exp.fam))} high +\nmedium families"),
        ("Position and models", "spanning vs local (ori = 0, ter = L/2);\nphylogenetic mixed model + geometry nulls",
         f"{core.fam['species_id'].nunique()} / {exp.fam['species_id'].nunique()} species\nin the models"),
    ]
    top, h, gap = 0.93, 0.128, 0.03
    for i, (title, body, count) in enumerate(steps):
        y = top - i * (h + gap)
        bx.add_patch(FancyBboxPatch((0.03, y - h), 0.94, h, boxstyle="round,pad=0,rounding_size=0.015",
                                    fc="#F4F6F8", ec=LIGHT, lw=0.6))
        bx.text(0.055, y - 0.018, f"{i + 1}  {title}", fontsize=6.6, fontweight="bold", va="top", color=INK)
        bx.text(0.055, y - 0.05, body, fontsize=5.7, va="top", color=INK, linespacing=1.15)
        bx.text(0.95, y - h / 2, count, fontsize=5.9, va="center", ha="right", color=INK, linespacing=1.15,
                fontweight="bold")
        if i < len(steps) - 1:
            bx.annotate("", xy=(0.5, y - h - gap + 0.003), xytext=(0.5, y - h - 0.003),
                        arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.6, mutation_scale=5))
    bx.text(0.97, 0.985, "counts: core / expanded", fontsize=5.8, color=GREY, ha="right", va="top")
    fig.savefig(out, dpi=600)
    plt.close(fig)


# ------------------------------------------------------------------ figure 2 / S1

def fig2(d, out, title):
    fam, sp = d.fam, d.sp
    fig, axes = plt.subplots(1, 2, figsize=(2 * W + 0.3, 3.1), sharey=True)
    ax = axes[0]
    ax.scatter(fam["x"], fam["y"], s=3, color=LIGHT, lw=0, zorder=2, label=f"Inversion family (n = {len(fam):,})")
    cap = sp.assign(x=np.log10(sp[GENOME] / 1e6), y=np.log10(sp["reference_length_bp"] / 2 / 1000))
    ax.scatter(cap["x"], cap["y"], marker="_", s=14, color=INK, lw=0.8, zorder=3,
               label="Largest possible size (half the chromosome)")
    g = species_means(fam)
    main = d.model[(d.model["analysis"] == "main") & (d.model["model"] == "PMM")].iloc[0]
    room = d.room("slope_all")[0]
    x0, x1 = g["x"].min(), g["x"].max()
    xs, ys = fit_line(g, main["slope"], x0, x1)
    ax.plot(xs, ys, color=INK, lw=1.6, zorder=4, label=f"Fitted slope {main['slope']:.2f}")
    xs, ys = fit_line(g, room, x0, x1)
    ax.plot(xs, ys, color=GREY, lw=1.1, ls=(0, (4, 2)), zorder=4, label=f"Room-alone null {room:.2f}")
    ax.set_ylim(np.log10(0.9), np.log10(6500))
    ax.set_xlim(fam["x"].min() - 0.03, fam["x"].max() + 0.03)
    log_axes(ax)
    ax.legend(loc="lower right", handlelength=1.6, borderaxespad=0.2, markerscale=1.5, frameon=True,
              facecolor="white", edgecolor="none", framealpha=0.92)
    ax.set_title(f"All inversions, {title}", loc="left", fontsize=7, color=INK)
    panel_label(ax, "a")

    ax = axes[1]
    for cls in ["spanning", "local"]:
        ev = fam[fam["position_class"] == cls]
        g = species_means(ev)
        ax.scatter(g["x"], g["y"], s=3 + 1.2 * g["n"], color=COL[cls], alpha=0.75, edgecolors="white", lw=0.3,
                   zorder=3)
        r = d.rob_row("as reported (terminus at half)", cls)
        room, lo, hi = d.room(f"slope_{cls}")
        x0, x1 = g["x"].min(), g["x"].max()
        xs, ys = fit_line(g, r["slope"], x0, x1)
        name = "Spanning" if cls == "spanning" else "Local"
        ax.plot(xs, ys, color=COL[cls], lw=1.8, zorder=4, solid_capstyle="round",
                label=f"{name} {r['slope']:+.2f} [{r['ci_low']:.2f}, {r['ci_high']:.2f}]")
        xs, ys = fit_line(g, room, x0, x1)
        ax.plot(xs, ys, color=COL[cls], lw=0.9, ls=(0, (4, 2)), zorder=4, label=f"{name}, room alone {room:+.2f}")
    ax.set_ylim(np.log10(0.9), np.log10(6500))
    ax.set_xlim(fam["x"].min() - 0.03, fam["x"].max() + 0.03)
    log_axes(ax, ylabel=None)
    ax.tick_params(labelleft=True)
    ax.legend(loc="lower right", ncol=1, handlelength=1.6, borderaxespad=0.2, frameon=True, facecolor="white",
              edgecolor="none", framealpha=0.92)
    ax.set_title("Split by position: mean size per species (circle area = families)", loc="left", fontsize=7,
                 color=INK)
    panel_label(ax, "b", x=-0.1)
    fig.tight_layout(w_pad=1.6)
    fig.savefig(out, dpi=600)
    plt.close(fig)


# ------------------------------------------------------------------ figure 3

def forest_rows(core, exp):
    rows = []

    def add(label, cls, getter, null=None, header=False):
        rows.append({"label": label, "cls": cls, "header": header,
                     "core": getter(core), "expanded": getter(exp),
                     "null_core": null(core) if null else None, "null_exp": null(exp) if null else None})

    def mm(d):
        r = d.model[(d.model["analysis"] == "main") & (d.model["model"] == "PMM")].iloc[0]
        return est(r["slope"], r["ci_low"], r["ci_high"])

    def rob(check, cls, setting=None):
        return lambda d: est(*d.rob_row(check, cls, setting)[["slope", "ci_low", "ci_high"]])

    def spl(cls, analysis):
        return lambda d: est(*d.split_row(cls, analysis)[["slope", "ci_low", "ci_high"]])

    def inter(d):
        r = d.inter[d.inter["quantity"] == "slope difference spanning - local"].iloc[0]
        return est(r["estimate"], r["ci_low"], r["ci_high"])

    def sym(which):
        return lambda d: est(*d.symsplit[d.symsplit["group"].str.contains(which)].iloc[0][["slope", "ci_low", "ci_high"]])

    def ot(group):
        def f(d):
            o = d.oriter[(d.oriter["terminus"] == "terminus at half") & (d.oriter["group"] == group)].iloc[0]
            return est(o["slope"], o["ci_low"], o["ci_high"])
        return f

    room = lambda q: (lambda d: d.room(q))  # noqa: E731
    rows.append({"label": "Main tests", "header": True})
    add("All inversions (pooled)", "all", mm, room("slope_all"))
    add("Spanning", "spanning", rob("as reported (terminus at half)", "spanning"), room("slope_spanning"))
    add("Local", "local", rob("as reported (terminus at half)", "local"), room("slope_local"))
    add("Spanning − local", "all", inter, room("slope_difference"))
    rows.append({"label": "Spanning: checks", "header": True})
    add("High tier only", "spanning", spl("spanning", "high_only"))
    add("All tiers (low included)", "spanning", spl("spanning", "all_tiers"))
    add("Chromosome length as predictor", "spanning", spl("spanning", "chromosome_length"))
    add("Terminus from GC skew", "spanning", rob("terminus from GC skew", "spanning"))
    add("Near-ori/ter families removed (5%)", "spanning", rob("margin around ori/ter", "spanning", 0.05))
    add("One species per genus", "spanning", rob("one species per genus (most families)", "spanning"))
    add("Largest 10% per species removed", "spanning",
        rob("largest 10% per species removed (mixed model)", "spanning"))
    add("10th percentile per species", "spanning", rob("percentile per species (PGLS)", "spanning", 0.1))
    add("Median per species", "spanning", rob("percentile per species (PGLS)", "spanning", 0.5))
    add("90th percentile per species", "spanning", rob("percentile per species (PGLS)", "spanning", 0.9))
    rows.append({"label": "Spanning: subgroups", "header": True})
    add("Symmetric (A < 0.5)", "spanning", sym("symmetric \\("))
    add("Asymmetric (A ≥ 0.5)", "spanning", sym("asymmetric"))
    add("Crossing the origin", "spanning", ot("ori"))
    add("Crossing the terminus", "spanning", ot("ter"))
    rows.append({"label": "Local: checks", "header": True})
    add("High tier only", "local", spl("local", "high_only"))
    add("Terminus from GC skew", "local", rob("terminus from GC skew", "local"))
    add("Largest 10% per species removed", "local", rob("largest 10% per species removed (mixed model)", "local"))
    add("Median per species", "local", rob("percentile per species (PGLS)", "local", 0.5))
    return rows


def fig3(core, exp, out):
    rows = forest_rows(core, exp)
    n = len(rows)
    fig, ax = plt.subplots(figsize=(W, 0.168 * n + 0.75))
    yy = np.arange(n)[::-1]
    off = {"core": 0.17, "expanded": -0.17}
    for y, r in zip(yy, rows):
        if r.get("header"):
            ax.text(-0.02, y, r["label"], transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=6.6, fontweight="bold", color=INK)
            continue
        ax.text(-0.02, y, r["label"], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.2,
                color=INK)
        c = COL[r["cls"]]
        for ds, nk in [("core", "null_core"), ("expanded", "null_exp")]:
            if r[nk] is not None:
                m, lo, hi = r[nk]
                ax.plot([lo, hi], [y + off[ds]] * 2, color=BAND, lw=3.2, solid_capstyle="butt", zorder=1)
                ax.plot([m], [y + off[ds]], marker="|", color=GREY, ms=4.5, mew=0.8, zorder=2)
            s, lo, hi = r[ds]
            ax.plot([lo, hi], [y + off[ds]] * 2, color=c, lw=0.9, zorder=3)
            mk = MARK[ds]
            ax.plot([s], [y + off[ds]], marker=mk["marker"], ms=mk["ms"], color=c,
                    mfc=c if mk["mfc"] == "fill" else "white", mew=0.8, zorder=4, ls="none")
    ax.axvline(0, color=INK, lw=0.6)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (3, 2)))
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_yticks([])
    ax.set_xlim(-1.3, 1.6)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1, 1.5])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("Slope [95% CI]")
    ax.text(1.0, n - 0.35, "proportional", fontsize=5.6, color=GREY, ha="center", va="bottom")
    h = [plt.Line2D([], [], marker="o", color=INK, mfc=INK, ls="-", lw=0.9, ms=4, label=MARK["core"]["label"]),
         plt.Line2D([], [], marker="D", color=INK, mfc="white", ls="-", lw=0.9, ms=3.4,
                    label=MARK["expanded"]["label"]),
         plt.Line2D([], [], color=BAND, lw=3.2, label="Room-alone null, 95% range")]
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, handlelength=1.4,
               columnspacing=0.8, fontsize=5.8)
    fig.subplots_adjust(left=0.45, right=0.98, top=0.985, bottom=0.115)
    fig.savefig(out, dpi=600)
    plt.close(fig)


# ------------------------------------------------------------------ figure 4

def fig4(core, exp, out, numbers):
    fig, axes = plt.subplots(2, 1, figsize=(W, 5.0))
    ax = axes[0]
    centres = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    for ds, d, dx in [("core", core, -0.02), ("expanded", exp, 0.02)]:
        s = d.symmetry[(d.symmetry["terminus"] == "terminus at half") & d.symmetry["asymmetry_bin"].str.match(r"0\.\d-")]
        mk = MARK[ds]
        ax.errorbar(centres + dx, s["obs_over_exp"], yerr=[s["obs_over_exp"] - s["ci_low_species"],
                                                           s["ci_high_species"] - s["obs_over_exp"]],
                    fmt=mk["marker"], ms=mk["ms"], color=VERM, mfc=VERM if mk["mfc"] == "fill" else "white",
                    mew=0.8, elinewidth=0.8, capsize=0, label=mk["label"])
    ax.axhline(1, color=GREY, lw=0.6, ls=(0, (3, 2)))
    ax.set_yscale("log")
    ax.set_yticks([0.25, 0.5, 1, 2, 4])
    ax.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.minorticks_off()
    ax.set_xticks(centres)
    ax.set_xticklabels(["0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1"])
    ax.set_xlabel("Asymmetry A = |d1 − d2| / (d1 + d2) of spanning families")
    ax.set_ylabel("Observed / expected (random placement)")
    ax.text(0.98, 1.08, "chance", fontsize=5.8, color=GREY, ha="right", va="bottom",
            transform=ax.get_yaxis_transform())
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="upper right", handlelength=1.2)
    ax.set_title("Spanning inversions are symmetric around ori or ter", loc="left", fontsize=7)
    panel_label(ax, "a")

    ax = axes[1]
    core_ids = set(core.tax["species_id"])
    g = share_by_species(exp)
    inc = g.index.isin(core_ids)
    ax.scatter(g.loc[inc, "x"], g.loc[inc, "share"], s=2 + 0.8 * g.loc[inc, "n"], color=VERM, alpha=0.7, lw=0,
               label="Species in the core set")
    ax.scatter(g.loc[~inc, "x"], g.loc[~inc, "share"], s=2 + 0.8 * g.loc[~inc, "n"], facecolors="white",
               edgecolors=VERM, lw=0.6, alpha=0.9, label="Species added in the expanded set")
    for ds, d, ls in [("core", core, "-"), ("expanded", exp, (0, (4, 2)))]:
        gg = share_by_species(d)
        r = numbers[(numbers["dataset"] == ds) & numbers["analysis"].str.startswith("spanning share")].iloc[0]
        xs = np.linspace(gg["x"].min(), gg["x"].max(), 20)
        b0 = gg["share"].mean() - r["slope"] * gg["x"].mean()
        ax.plot(xs, b0 + r["slope"] * xs, color=INK, lw=1.3, ls=ls,
                label=f"{'Core' if ds == 'core' else 'Expanded'}: +{r['slope'] * np.log10(2) * 100:.0f} points per "
                      f"doubling (P = {r['p']:.2g})")
    ax.set_ylim(-0.03, 1.03)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlim(g["x"].min() - 0.03, g["x"].max() + 0.03)
    log_axes(ax, ylabel="Share of families that span", ysize=False)
    ax.grid(axis="y", color="#EEF1F4", lw=0.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=2, handlelength=1.8, fontsize=5.6,
              columnspacing=0.8)
    ax.set_title("Bigger genomes have a larger share of spanning inversions", loc="left", fontsize=7)
    panel_label(ax, "b")
    fig.tight_layout(h_pad=1.2)
    fig.savefig(out, dpi=600)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--core-root", type=Path, required=True)
    ap.add_argument("--core-cohort", default="selected_50_v1")
    ap.add_argument("--core-species", type=Path, required=True)
    ap.add_argument("--core-tree", type=Path, required=True)
    ap.add_argument("--expanded-root", type=Path, required=True)
    ap.add_argument("--expanded-cohort", default="selected_min20_v1")
    ap.add_argument("--expanded-species", type=Path, required=True)
    ap.add_argument("--expanded-tree", type=Path, required=True)
    ap.add_argument("--extra-core", type=Path, help="folder with ori_ter.tsv for the core set, if elsewhere")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    set_style()
    core = Data(a.core_root, a.core_cohort, a.core_species, a.core_tree)
    exp = Data(a.expanded_root, a.expanded_cohort, a.expanded_species, a.expanded_tree)
    if core.oriter is None and a.extra_core:
        core.oriter = pd.read_csv(a.extra_core / "ori_ter.tsv", sep="\t")
    numbers = extra_numbers(core, exp)
    numbers.to_csv(a.out / "paper_numbers.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(numbers.round(4).to_string(index=False))
    fig1(core, exp, a.out / "fig1_design.png")
    fig2(core, a.out / "fig2_core.png", "core dataset")
    fig2(exp, a.out / "figS1_expanded.png", "expanded dataset")
    fig3(core, exp, a.out / "fig3_slopes.png")
    fig4(core, exp, a.out / "fig4_position.png", numbers)
    print(f"Written to {a.out}")


if __name__ == "__main__":
    main()
