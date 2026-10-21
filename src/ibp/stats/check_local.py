"""
Checks on the negative slope of LOCAL inversions (one replication arm).

1. Same model on the high tier only.
2. Leave one phylum out.
3. Leave one species out: which species move the slope most.
4. Only species with >= 3 local families.
5. Counts: local and spanning families per species vs genome size
   (do large genomes just have MORE small local inversions?).
6. Figure: every local family (dots) and species means (circles sized by
   number of families), coloured by phylum, with the fitted line and the S0
   line through the same mean.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ibp.analysis_legacy.final_inversion_model import parse_newick_tree
from ibp.stats.model import fit_all, read_s0, species_data

GENOME = "median_genome_length_bp"


def pmm(d, tree, s0):
    r = fit_all(d, tree, 3, s0)["PMM"]
    return {"n_species": r["n_species"], "n_families": r["n_events"],
            "slope": round(r["slope"], 3), "ci_low": round(r["slope_ci_low"], 3),
            "ci_high": round(r["slope_ci_high"], 3),
            "p_vs_0": r["p_slope_vs_0"], "p_vs_S0": r.get(f"p_slope_vs_{s0:g}") if s0 is not None else np.nan}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--loci", type=Path, required=True)
    p.add_argument("--tree", type=Path, required=True)
    p.add_argument("--null", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    fam = pd.read_csv(a.loci / "families.tsv", sep="\t")
    sp = pd.read_csv(a.loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    tree = parse_newick_tree(a.tree.read_text().strip())
    s0 = read_s0(a.null).get("local", (None,))[0]
    hm = ["high", "medium"]

    rows = []
    d = species_data(fam, sp, hm, "local", "size_bp", GENOME)
    rows.append({"check": "local, high+medium (as reported)", **pmm(d, tree, s0)})
    rows.append({"check": "local, high only",
                 **pmm(species_data(fam, sp, ["high"], "local", "size_bp", GENOME), tree, s0)})
    rows.append({"check": "local, species with >= 3 local families", **pmm(d[d["n"] >= 3], tree, s0)})
    for ph in sorted(d["phylum"].dropna().unique()):
        if (d["phylum"] != ph).sum() >= 10:
            rows.append({"check": f"local, without {ph} ({int((d['phylum'] == ph).sum())} sp.)",
                         **pmm(d[d["phylum"] != ph], tree, s0)})
    checks = pd.DataFrame(rows)

    loso = []
    for s in d.index:
        r = pmm(d.drop(index=s), tree, s0)
        loso.append({"species_left_out": s, "n_local": int(d.at[s, "n"]),
                     "genome_Mb": round(10 ** d.at[s, "x"], 2), "slope_without": r["slope"]})
    loso = pd.DataFrame(loso)
    base = checks.iloc[0]["slope"]
    loso["change"] = (loso["slope_without"] - base).round(3)
    loso = loso.reindex(loso["change"].abs().sort_values(ascending=False).index)

    cnt = (fam[fam["tier"].isin(hm)].groupby(["species_id", "position_class"]).size()
           .unstack(fill_value=0).reindex(sp["species_id"], fill_value=0))
    cnt = cnt.join(sp.set_index("species_id")[[c for c in [GENOME, "phylum", "scientific_name"] if c in sp]])
    cnt["genome_Mb"] = cnt[GENOME] / 1e6
    cnt["local_share"] = (cnt.get("local", 0) / (cnt.get("local", 0) + cnt.get("spanning", 0))).round(2)
    lx = np.log10(cnt["genome_Mb"])
    corr = {}
    for col in ["local", "spanning", "local_share"]:
        ok = cnt[col].notna()
        corr[col] = round(float(np.corrcoef(lx[ok], cnt.loc[ok, col])[0, 1]), 2)

    size_by_bin = (fam[fam["tier"].isin(hm) & (fam["position_class"] == "local")]
                   .join(sp.set_index("species_id")[GENOME], on="species_id")
                   .assign(genome_bin=lambda t: pd.cut(t[GENOME] / 1e6, [0, 2, 3, 4, 5, 6, 10]))
                   .groupby("genome_bin", observed=True)["size_bp"]
                   .agg(families="size", median_kb=lambda v: round(v.median() / 1000, 1),
                        q90_kb=lambda v: round(v.quantile(0.9) / 1000, 1)))

    checks.to_csv(a.out / "local_checks.tsv", sep="\t", index=False)
    loso.to_csv(a.out / "local_leave_one_species_out.tsv", sep="\t", index=False)
    cnt.to_csv(a.out / "families_per_species_by_position.tsv", sep="\t")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        show = checks.copy()
        for c in ["p_vs_0", "p_vs_S0"]:
            show[c] = show[c].map(lambda v: f"{v:.2g}" if pd.notna(v) else "")
        print(f"S0 for local inversions: {s0}\n")
        print(show.to_string(index=False))
        print("\nSpecies that move the local slope most (leave one out):")
        print(loso.head(8).to_string(index=False))
        print("\nLocal family sizes by genome size bin:")
        print(size_by_bin.to_string())
        print(f"\nCorrelation with log genome size: local families {corr['local']}, "
              f"spanning families {corr['spanning']}, local share {corr['local_share']}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ev = (fam[fam["tier"].isin(hm) & (fam["position_class"] == "local")]
              .join(sp.set_index("species_id")[[GENOME, "phylum"]], on="species_id"))
        fig, ax = plt.subplots(figsize=(7, 5))
        phyla = sorted(d["phylum"].dropna().unique())
        colors = dict(zip(phyla, plt.cm.tab10.colors))
        for ph in phyla:
            e = ev[ev["phylum"] == ph]
            ax.scatter(np.log10(e[GENOME] / 1e6), np.log10(e["size_bp"] / 1000), s=8, alpha=0.35,
                       color=colors[ph])
            m = d[d["phylum"] == ph]
            ax.scatter(m["x"], m["ybar"], s=12 + 6 * m["n"], facecolors="none", edgecolors=colors[ph],
                       linewidths=1.2, label=f"{ph} ({len(m)})")
        xs = np.linspace(d["x"].min(), d["x"].max(), 50)
        xm, ym = np.average(d["x"], weights=d["n"]), np.average(d["ybar"], weights=d["n"])
        ax.plot(xs, ym + base * (xs - xm), color="black", lw=2, label=f"fitted slope {base:.2f}")
        if s0 is not None:
            ax.plot(xs, ym + s0 * (xs - xm), color="grey", lw=2, ls="--",
                    label=f"geometry only (S0) {s0:.2f}")
        ax.set_xlabel("Genome size, log10 Mb")
        ax.set_ylabel("Local inversion size, log10 kb")
        ax.set_title("Local inversions (one replication arm), high + medium")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(a.out / "local_inversions.png", dpi=150)
        print(f"\nFigure: {a.out / 'local_inversions.png'}")
    except ImportError:
        print("matplotlib not available; no figure.")


if __name__ == "__main__":
    main()
