"""
Figures comparing this rerun (GTDB, loci, null model) with the team run.

1. slopes.png        slope estimates side by side (team vs rerun) for the main
                     model, PGLS, chromosome length, spanning and local, with
                     the rerun's geometry-only band (S0) and proportional
                     scaling (S1 = 1).
2. species_boxes.png every species as a box of its family sizes, ordered by
                     genome size and coloured by phylum, with the fitted slope
                     and the S0 line (the supervisor's per-species graph).
3. position.png      share of spanning vs local families: observed vs the
                     geometry null.
4. raw_backing.png   share of loci backed by the raw alignment per evidence
                     class, next to the team's Part 6 numbers.
5. summary.tsv       the numbers behind 1, 3 and 4.

Team numbers come from config/team_run/team_results.tsv (edit there).
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TEAM = "#8c8c8c"
RERUN = "#1f5fa8"
S0_COLOR = "#e6a23c"


def load_rerun(model_dir):
    res = pd.read_csv(model_dir / "model_results.tsv", sep="\t")
    null = pd.read_csv(model_dir / "null_slopes.tsv", sep="\t")
    null = null[(null["size_distribution"] == "empirical_pooled_sizes") & (null["fit"] == "events")]
    return res, null.set_index("subset")


def slope_rows(res, team):
    def pick(analysis, model):
        r = res[(res["analysis"] == analysis) & (res["model"] == model)].iloc[0]
        return r["slope"], r["ci_low"], r["ci_high"], int(r["n_species"])
    t = team.set_index("measure")
    spec = [("Main model (mixed, tree)", "PMM_main", ("main", "PMM"), "all"),
            ("Species means (PGLS)", "PGLS_main", ("main", "PGLS"), "all"),
            ("Chromosome length", "PMM_chromosome_length", ("sens_chromosome_length", "PMM"), "all"),
            ("Spanning ori/ter", "spanning", ("expl_spanning", "PMM"), "spanning"),
            ("Local (one arm)", "local", ("expl_local", "PMM"), "local")]
    rows = []
    for label, team_key, ours, subset in spec:
        tr = t.loc[team_key] if team_key in t.index else None
        est, lo, hi, n = pick(*ours)
        rows.append({"label": label, "subset": subset,
                     "team": tr["estimate"] if tr is not None else np.nan,
                     "team_lo": tr["ci_low"] if tr is not None else np.nan,
                     "team_hi": tr["ci_high"] if tr is not None else np.nan,
                     "team_n": tr["n_species"] if tr is not None else np.nan,
                     "rerun": est, "rerun_lo": lo, "rerun_hi": hi, "rerun_n": n})
    return pd.DataFrame(rows)


def plot_slopes(rows, null, out):
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    y = np.arange(len(rows))[::-1]
    for yi, (_, r) in zip(y, rows.iterrows()):
        if r["subset"] in null.index:
            s = null.loc[r["subset"]]
            ax.fill_betweenx([yi - 0.38, yi + 0.38], s["slope_2.5"], s["slope_97.5"], color=S0_COLOR, alpha=0.25, lw=0)
            ax.plot([s["slope_mean"]] * 2, [yi - 0.38, yi + 0.38], color=S0_COLOR, lw=1.5)
        if pd.notna(r["team"]):
            if pd.notna(r["team_lo"]):
                ax.plot([r["team_lo"], r["team_hi"]], [yi + 0.15] * 2, color=TEAM, lw=2)
            ax.plot(r["team"], yi + 0.15, "o", color=TEAM, ms=7)
        ax.plot([r["rerun_lo"], r["rerun_hi"]], [yi - 0.15] * 2, color=RERUN, lw=2)
        ax.plot(r["rerun"], yi - 0.15, "o", color=RERUN, ms=7)
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(1, color="black", lw=0.8, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels(rows["label"])
    ax.set_xlabel("Slope: log10 inversion size on log10 genome size")
    ax.plot([], [], "o-", color=TEAM, label="Team run (NCBI species, events)")
    ax.plot([], [], "o-", color=RERUN, label="Rerun (GTDB species, loci, 95% CI)")
    ax.fill_between([], [], color=S0_COLOR, alpha=0.4, label="Geometry only (S0, 95% range)")
    ax.plot([], [], ":", color="black", label="Proportional (S1 = 1)")
    ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=2)
    ax.set_title("Team run vs rerun, against what geometry alone produces", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_species_boxes(fam, sp, res, null, out):
    hm = fam[fam["tier"].isin(["high", "medium"])].copy()
    hm = hm.join(sp.set_index("species_id")[["median_genome_length_bp", "phylum", "scientific_name"]],
                 on="species_id", how="inner")
    hm["y"] = np.log10(hm["size_bp"] / 1000)
    order = (hm.groupby("species_id")["median_genome_length_bp"].first().sort_values().index.tolist())
    phyla = sorted(hm["phylum"].dropna().unique())
    colors = dict(zip(phyla, plt.cm.tab10.colors))
    fig, ax = plt.subplots(figsize=(12, 5))
    xs = np.log10(hm.groupby("species_id")["median_genome_length_bp"].first() / 1e6)
    data = [hm.loc[hm["species_id"] == s, "y"].to_numpy() for s in order]
    pos = xs.loc[order].to_numpy()
    width = (pos.max() - pos.min()) / len(order) * 0.9
    bp = ax.boxplot(data, positions=pos, widths=width, patch_artist=True, showfliers=False,
                    manage_ticks=False, medianprops={"color": "black", "lw": 1.2})
    for patch, s in zip(bp["boxes"], order):
        patch.set_facecolor(colors.get(hm.loc[hm["species_id"] == s, "phylum"].iloc[0], "grey"))
        patch.set_alpha(0.6)
    ax.scatter(np.log10(hm["median_genome_length_bp"] / 1e6), hm["y"], s=4, color="black", alpha=0.25)
    main = res[(res["analysis"] == "main") & (res["model"] == "PMM")].iloc[0]
    xline = np.linspace(pos.min(), pos.max(), 50)
    g = hm.groupby("species_id").agg(x=("median_genome_length_bp", lambda v: np.log10(v.iloc[0] / 1e6)),
                                     y=("y", "mean"), n=("y", "size"))
    xm, ym = np.average(g["x"], weights=g["n"]), np.average(g["y"], weights=g["n"])
    ax.plot(xline, ym + main["slope"] * (xline - xm), color=RERUN, lw=2.5,
            label=f"Fitted (mixed model) {main['slope']:.2f}")
    s0 = null.loc["all", "slope_mean"]
    ax.plot(xline, ym + s0 * (xline - xm), color=S0_COLOR, lw=2.5, ls="--", label=f"Geometry only (S0) {s0:.2f}")
    from matplotlib.patches import Patch
    handles, _ = ax.get_legend_handles_labels()
    handles += [Patch(facecolor=colors[ph], alpha=0.6,
                      label=f"{ph} ({hm.loc[hm['phylum'] == ph, 'species_id'].nunique()})") for ph in phyla]
    pad = (pos.max() - pos.min()) * 0.03
    ax.set_xlim(pos.min() - pad, pos.max() + pad)
    ax.set_xlabel("Genome size, log10 Mb (one box per species)")
    ax.set_ylabel("Inversion family size, log10 kb")
    ax.set_title(f"Inversion size per species ({len(order)} GTDB species, high + medium families)", fontsize=11)
    ax.legend(handles=handles, fontsize=8, frameon=False, ncol=2, loc="lower left")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_position(fam, null, out):
    hm = fam[fam["tier"].isin(["high", "medium"])]
    obs = (hm["position_class"] == "spanning").mean()
    exp = null.loc["spanning_fraction"]
    fig, ax = plt.subplots(figsize=(4.2, 4))
    ax.bar([0], [obs * 100], color=RERUN, width=0.6)
    ax.bar([1], [exp["slope_mean"] * 100], color=S0_COLOR, width=0.6,
           yerr=[[100 * (exp["slope_mean"] - exp["slope_2.5"])], [100 * (exp["slope_97.5"] - exp["slope_mean"])]],
           capsize=6)
    ax.axhline(50, color="black", lw=0.8, ls=":")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Observed", "Geometry only"])
    ax.set_ylabel("Families spanning origin or terminus, %")
    ax.set_ylim(0, 100)
    ax.text(0, obs * 100 + 2, f"{obs:.0%}", ha="center")
    ax.text(1, exp["slope_97.5"] * 100 + 2, f"{exp['slope_mean']:.0%}", ha="center")
    ax.set_title(f"Excess of spanning inversions\n({len(hm)} high + medium families)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return obs, exp["slope_mean"]


def plot_raw(loci, team, out):
    if "raw_support" not in loci:
        return pd.DataFrame()
    loc = loci.dropna(subset=["raw_support"]).copy()
    loc["raw_support"] = loc["raw_support"].astype(bool)
    two = loc["n_leen"].gt(0) & loc["n_syri"].gt(0)
    classes = [("Leen + SyRI\nmatching", two & loc["breakpoints_agree"]),
               ("Leen + SyRI\noverlap only", two & ~loc["breakpoints_agree"]),
               ("Leen only", loc["n_leen"].gt(0) & loc["n_syri"].eq(0)),
               ("SyRI only", loc["n_leen"].eq(0))]
    rows = [{"class": c.replace("\n", " "), "loci": int(m.sum()),
             "raw_backed": float(loc.loc[m, "raw_support"].mean()) if m.any() else np.nan} for c, m in classes]
    rows = pd.DataFrame(rows)
    t = team.set_index("measure")
    team_vals = [t.loc[k, "estimate"] if k in t.index else np.nan
                 for k in ["raw_backing_high", "raw_backing_medium", "raw_backing_low"]]
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(rows))
    ax.bar(x, rows["raw_backed"] * 100, color=RERUN, width=0.6, label="Rerun: loci by evidence")
    for xi, (_, r) in zip(x, rows.iterrows()):
        ax.text(xi, r["raw_backed"] * 100 + 2, f"{r['raw_backed']:.0%}\nn={r['loci']}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([c for c, _ in classes], fontsize=8)
    for val, lab, ls in zip(team_vals, ["team high", "team medium", "team low"], ["-", "--", ":"]):
        if pd.notna(val):
            ax.axhline(val * 100, color=TEAM, lw=1.2, ls=ls, label=f"{lab} tier: {val:.0%}")
    ax.set_ylim(0, 115)
    ax.set_ylabel("Backed by the raw alignment, %")
    ax.set_title("Raw-alignment check per evidence class", fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--loci", type=Path, required=True)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--team", type=Path, default=Path("config/team_run/team_results.tsv"))
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    fam = pd.read_csv(a.loci / "families.tsv", sep="\t")
    loci = pd.read_csv(a.loci / "loci.tsv", sep="\t")
    sp = pd.read_csv(a.loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    team = pd.read_csv(a.team, sep="\t")
    res, null = load_rerun(a.model)

    rows = slope_rows(res, team)
    plot_slopes(rows, null, a.out / "slopes.png")
    plot_species_boxes(fam, sp, res, null, a.out / "species_boxes.png")
    obs, exp = plot_position(fam, null, a.out / "position.png")
    raw = plot_raw(loci, team, a.out / "raw_backing.png")

    with open(a.out / "summary.tsv", "w") as h:
        rows.to_csv(h, sep="\t", index=False)
        h.write(f"\nspanning_share_observed\t{obs:.3f}\nspanning_share_geometry\t{exp:.3f}\n\n")
        if len(raw):
            raw.to_csv(h, sep="\t", index=False)
    with pd.option_context("display.width", 200):
        print(rows.round(3).to_string(index=False))
        print(f"\nSpanning share: observed {obs:.0%}, geometry only {exp:.0%}")
        if len(raw):
            print("\n" + raw.round(3).to_string(index=False))
    print(f"\nFigures in {a.out}: slopes.png, species_boxes.png, position.png, raw_backing.png")


if __name__ == "__main__":
    main()
