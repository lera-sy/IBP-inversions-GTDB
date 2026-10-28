"""
Is the size distribution normal, and does it matter? (lab-head question, 28 Oct)

The mixed model does not need inversion sizes to be normal; it assumes that
log10 size is normal AROUND the fitted line: (a) within-species deviations and
(b) species-level deviations. Reported:

1. Shape of log10 size: pooled, spanning, local (skewness, excess kurtosis,
   Shapiro-Wilk, number of modes is visible in the figure).
2. Model residuals per group: within-species (y - species mean) and species-level
   (species mean - fitted line), Shapiro-Wilk + skew/kurtosis, QQ plots.
3. Checks that do not assume normality:
   - Spearman rank correlation of species median size with genome size;
   - species bootstrap (resample species with replacement, refit the mixed model)
     percentile CI of the slope;
   (the empirical room-alone null and the percentile slopes elsewhere are also
   distribution-free with respect to the size distribution).
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

from ibp.analysis_legacy.final_inversion_model import fit_pmm
from ibp.stats.robustness import GENOME, HM, Tree, summarise

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COL = {"all": "#7D8C99", "spanning": "#D55E00", "local": "#0072B2"}


def shape(y):
    y = np.asarray(y, float)
    sh = stats.shapiro(y if len(y) <= 5000 else np.random.default_rng(1).choice(y, 5000, replace=False))
    return {"n": len(y), "mean": y.mean(), "sd": y.std(ddof=1), "skew": stats.skew(y),
            "excess_kurtosis": stats.kurtosis(y), "shapiro_W": sh.statistic, "shapiro_p": sh.pvalue}


def pmm_fit(d, tree, idx=None):
    C = tree.sub(d.index) if idx is None else tree.C[np.ix_(idx, idx)]
    return fit_pmm(d["ybar"].to_numpy(), d["n"].to_numpy(float), float(d["ssw"].sum()), d["x"].to_numpy(), C,
                   [0.0], phylogenetic=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(a.loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    fam = pd.read_csv(a.loci / "families.tsv", sep="\t")
    fam = fam[fam["tier"].isin(HM) & fam["species_id"].isin(sp["species_id"]) & (fam["size_bp"] > 0)].copy()
    fam["y"] = np.log10(fam["size_bp"] / 1000)
    fam["x"] = np.log10(fam["species_id"].map(sp.set_index("species_id")[GENOME]).astype(float) / 1e6)
    tree = Tree(a.tree, sp["species_id"].tolist())
    rng = np.random.default_rng(a.seed)

    shapes, resid_rows, np_rows, resid = [], [], [], {}
    for grp in ["all", "spanning", "local"]:
        ev = fam if grp == "all" else fam[fam["position_class"] == grp]
        ev = ev[ev["species_id"].isin(tree.pos)]
        shapes.append({"group": grp, "what": "log10 size (raw)", **shape(ev["y"])})
        d = summarise(ev, sp)
        d = d[d.index.isin(tree.pos)]
        r = pmm_fit(d, tree)
        within = (ev["y"] - ev["species_id"].map(d["ybar"])).to_numpy()
        within = within[ev["species_id"].map(d["n"]).to_numpy() > 1]
        between = (d["ybar"] - (r["intercept"] + r["slope"] * d["x"])).to_numpy()
        resid[grp] = (within, between)
        resid_rows.append({"group": grp, "what": "within-species residuals", **shape(within)})
        resid_rows.append({"group": grp, "what": "species-level residuals", **shape(between)})
        # distribution-free checks
        med = ev.groupby("species_id")["y"].median()
        xs = d.loc[med.index, "x"]
        rho, p_rho = stats.spearmanr(xs, med)
        ids = list(d.index)
        pos = np.array([tree.pos[s] for s in ids])
        boots = []
        for _ in range(a.boot):
            pick = rng.integers(0, len(ids), len(ids))
            db = d.iloc[pick]
            try:
                boots.append(pmm_fit(db, tree, idx=pos[pick])["slope"])
            except Exception:  # noqa: BLE001 - a degenerate resample is skipped
                continue
        boots = np.array(boots)
        np_rows.append({"group": grp, "model_slope": r["slope"], "model_ci": f"{r['slope_ci_low']:.2f} to {r['slope_ci_high']:.2f}",
                        "bootstrap_ci": f"{np.quantile(boots, 0.025):.2f} to {np.quantile(boots, 0.975):.2f}",
                        "bootstrap_share_above_0": float((boots > 0).mean()), "boot_reps": len(boots),
                        "spearman_rho_species_median": rho, "spearman_p": p_rho, "n_species": len(d)})

    shapes, resid_t, np_t = pd.DataFrame(shapes), pd.DataFrame(resid_rows), pd.DataFrame(np_rows)
    pd.concat([shapes, resid_t]).to_csv(a.out / "normality.tsv", sep="\t", index=False)
    np_t.to_csv(a.out / "distribution_free.tsv", sep="\t", index=False)

    # figure
    plt.rcParams.update({"font.size": 11})
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    bins = np.linspace(fam["y"].min() - 0.1, fam["y"].max() + 0.1, 40)
    ax = axes[0, 0]
    for grp in ["spanning", "local"]:
        ax.hist(fam.loc[fam["position_class"] == grp, "y"], bins=bins, color=COL[grp], alpha=0.6, label=grp)
    ax.set_title("Pooled log10 size: two groups, two peaks", loc="left", fontsize=11)
    ax.legend(frameon=False)
    for ax, grp in zip(axes[0, 1:], ["spanning", "local"]):
        y = fam.loc[fam["position_class"] == grp, "y"]
        ax.hist(y, bins=bins, color=COL[grp], alpha=0.7, density=True)
        xs = np.linspace(bins[0], bins[-1], 200)
        ax.plot(xs, stats.norm.pdf(xs, y.mean(), y.std()), color="black", lw=1.5)
        s = shapes[(shapes["group"] == grp)].iloc[0]
        ax.set_title(f"{grp}: skew {s['skew']:.2f}, excess kurtosis {s['excess_kurtosis']:.2f}", loc="left",
                     fontsize=11)
    for ax in axes[0]:
        ax.set_xlabel("log10 inversion size (kb)")
    for ax, grp in zip(axes[1], ["spanning", "local", "all"]):
        within, between = resid[grp]
        for vals, mk, lab in [(within, "o", "within species"), (between, "s", "species level")]:
            (osm, osr), _ = stats.probplot(vals / vals.std(), dist="norm")
            ax.scatter(osm, osr, s=10 if mk == "o" else 26, marker=mk, alpha=0.6,
                       color=COL[grp] if mk == "o" else "black", label=lab)
        ax.plot([-3.5, 3.5], [-3.5, 3.5], color="grey", lw=1, ls="--")
        ax.set_title(f"QQ plot of model residuals: {grp}", loc="left", fontsize=11)
        ax.set_xlabel("normal quantile")
        ax.set_ylabel("standardised residual")
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(a.out / "normality.png", dpi=180)
    plt.close(fig)

    pd.set_option("display.width", 200)
    print("== Shape of log10 size and of the model residuals")
    print(pd.concat([shapes, resid_t])[["group", "what", "n", "skew", "excess_kurtosis", "shapiro_W", "shapiro_p"]]
          .round(3).to_string(index=False))
    print("\n== Checks that do not assume normality")
    print(np_t.round(4).to_string(index=False))
    print(f"\nWritten to {a.out}: normality.tsv, distribution_free.tsv, normality.png")


if __name__ == "__main__":
    main()
