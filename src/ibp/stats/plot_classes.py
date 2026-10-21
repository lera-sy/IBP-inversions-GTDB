"""
Spanning vs local inversions: is each slope significant, against which reference,
and do the two classes differ?

For each class (and all inversions together) the main mixed model (PMM) and the
species-mean model (PGLS) are shown with their 95% CI, against three references:
  0    no relationship (the test the team used)
  S0   geometry only (shuffled sizes, same cap and floor), its 95% range as a band
  1    proportional scaling
p-values for each reference are written next to the estimate.

Class difference: spanning slope minus local slope, with a normal-approximation
test (standard errors from the CIs; the two classes come from the same species,
so this is approximate). Compared with the difference geometry alone gives
(S0 spanning - S0 local).
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COL = {"spanning": "#4c72b0", "local": "#c44e52", "all": "#555555"}
S0C = "#e08a1e"


def fmt_p(p):
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.2g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    res = pd.read_csv(a.model / "model_results.tsv", sep="\t")
    null = pd.read_csv(a.model / "null_slopes.tsv", sep="\t")
    null = null[(null["size_distribution"] == "empirical_pooled_sizes") & (null["fit"] == "events")].set_index("subset")

    rows = []
    for cls, analysis in [("all", "main"), ("spanning", "expl_spanning"), ("local", "expl_local")]:
        for model in ["PMM", "PGLS"]:
            r = res[(res["analysis"] == analysis) & (res["model"] == model)]
            if r.empty:
                continue
            r = r.iloc[0]
            df = r["n_species"] - 2
            se = (r["ci_high"] - r["ci_low"]) / (2 * stats.t.ppf(0.975, df))
            rows.append({"class": cls, "model": model, "n_species": int(r["n_species"]),
                         "n_families": r["n_families"], "slope": r["slope"], "ci_low": r["ci_low"],
                         "ci_high": r["ci_high"], "se": se, "S0": null.loc[cls, "slope_mean"],
                         "S0_low": null.loc[cls, "slope_2.5"], "S0_high": null.loc[cls, "slope_97.5"],
                         "p_vs_0": r["p_vs_0"], "p_vs_S0": r["p_vs_S0"], "p_vs_1": r["p_vs_1"]})
    t = pd.DataFrame(rows)

    diff = []
    for model in ["PMM", "PGLS"]:
        s = t[(t["class"] == "spanning") & (t["model"] == model)]
        lo = t[(t["class"] == "local") & (t["model"] == model)]
        if s.empty or lo.empty:
            continue
        s, lo = s.iloc[0], lo.iloc[0]
        d, se = s["slope"] - lo["slope"], np.hypot(s["se"], lo["se"])
        d0 = s["S0"] - lo["S0"]
        diff.append({"model": model, "difference": d, "se": se, "ci_low": d - 1.96 * se, "ci_high": d + 1.96 * se,
                     "p_vs_0": 2 * stats.norm.sf(abs(d / se)), "geometry_difference": d0,
                     "p_vs_geometry": 2 * stats.norm.sf(abs((d - d0) / se))})
    diff = pd.DataFrame(diff)

    # figure
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    order = [("all", "PMM"), ("all", "PGLS"), ("spanning", "PMM"), ("spanning", "PGLS"),
             ("local", "PMM"), ("local", "PGLS")]
    labels = {"all": "All inversions", "spanning": "Spanning ori/ter", "local": "Local (one arm)"}
    y = 0
    yticks, ylabels = [], []
    for cls, model in order:
        r = t[(t["class"] == cls) & (t["model"] == model)]
        if r.empty:
            continue
        r = r.iloc[0]
        ax.fill_betweenx([y - 0.35, y + 0.35], r["S0_low"], r["S0_high"], color=S0C, alpha=0.22, lw=0)
        ax.plot([r["S0"]] * 2, [y - 0.35, y + 0.35], color=S0C, lw=1.6)
        ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=COL[cls], lw=2.5)
        ax.plot(r["slope"], y, "o" if model == "PMM" else "s", color=COL[cls], ms=8)
        ax.text(1.62, y, f"{r['slope']:.2f}   vs 0: {fmt_p(r['p_vs_0'])}   vs S0: {fmt_p(r['p_vs_S0'])}",
                va="center", fontsize=8)
        yticks.append(y)
        ylabels.append(f"{labels[cls]} - {'mixed model' if model == 'PMM' else 'species means'}")
        y -= 1 if model == "PMM" else 1.6
    ax.axvline(0, color="black", lw=0.9)
    ax.axvline(1, color="black", lw=0.9, ls=":")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_xlim(-1.3, 1.6)
    ax.set_xlabel("Slope: log10 inversion size on log10 genome size (95% CI)")
    ax.plot([], [], color=S0C, lw=6, alpha=0.4, label="Geometry only (S0, 95% range)")
    ax.plot([], [], ":", color="black", label="Proportional (1)")
    ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.12), ncol=2)
    if len(diff):
        d = diff[diff["model"] == "PMM"].iloc[0]
        ax.set_title(f"Spanning vs local: difference {d['difference']:.2f} [{d['ci_low']:.2f}, {d['ci_high']:.2f}], "
                     f"{fmt_p(d['p_vs_0'])}; geometry alone gives {d['geometry_difference']:.2f}", fontsize=10)
    fig.subplots_adjust(left=0.27, right=0.72, bottom=0.2)
    fig.savefig(a.out / "classes_significance.png", dpi=200)
    plt.close(fig)

    t.to_csv(a.out / "classes_significance.tsv", sep="\t", index=False)
    diff.to_csv(a.out / "classes_difference.tsv", sep="\t", index=False)
    show = t.drop(columns=["se"]).copy()
    for c in ["p_vs_0", "p_vs_S0", "p_vs_1"]:
        show[c] = show[c].map(lambda v: f"{v:.2g}")
    with pd.option_context("display.width", 200):
        print(show.round(3).to_string(index=False))
        print("\nSpanning minus local:")
        print(diff.round(4).to_string(index=False))
    print(f"\nFigure: {a.out / 'classes_significance.png'}")


if __name__ == "__main__":
    main()
