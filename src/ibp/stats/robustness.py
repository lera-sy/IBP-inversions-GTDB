"""
Robustness of the spanning-vs-local result (reviewer feedback, 24 Oct 2026).

All checks reuse existing outputs (families, species table, tree, reference
FASTAs); nothing is called again. Each check refits the main model
(phylogenetic mixed model, high + medium families, log size on log genome
size) separately for spanning and local inversions.

1. Terminus position (the spanning/local split depends on it)
   a. shift:    assumed terminus moved by -15, -10, -5, +5, +10, +15 % of the chromosome
   b. margin:   families with a breakpoint within 2.5 / 5 / 10 % of the
                chromosome from the origin or terminus are left out
   c. inferred: terminus taken from each reference's cumulative GC skew
                (origin is at 0 after dnaapler; terminus = maximum of the
                cumulative skew). A dif-site search is not done.
2. Random-breakpoint null: per species the same number of families on the
   real chromosome length, both breakpoints uniform on the circle, shorter
   arc, 6 kb floor, classified the same way, same model. Gives what slope and
   what spanning share geometry alone produces when sizes have NO fixed
   distribution (sizes scale with the chromosome by construction, so the
   overall slope is near 1). Together with the shuffle null (real sizes, room
   only) it brackets the result.
3. Shared ancestry across species: one species per genus (and per GTDB family),
   (a) the species with most families, (b) 200 random picks.
4. Whole distribution (supervisor / reviewers: not just a few outliers):
   per species the 10/25/50/75/90th percentile of log size regressed on genome
   size (PGLS), and the main model after removing each species' largest 10%.
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from ibp.analysis_legacy.final_inversion_model import fit_pgls, fit_pmm, parse_newick_tree, tip_covariance

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GENOME = "median_genome_length_bp"
HM = ["high", "medium"]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

class Tree:
    """Tree covariance for all species once; subsets are sliced, not recomputed."""

    def __init__(self, path, species_ids):
        parent, length, tips = parse_newick_tree(Path(path).read_text().strip())
        self.ids = [s for s in species_ids if s in tips]
        self.C = tip_covariance(parent, length, [tips[s] for s in self.ids])
        self.pos = {s: i for i, s in enumerate(self.ids)}

    def sub(self, ids):
        idx = [self.pos[s] for s in ids]
        return self.C[np.ix_(idx, idx)]


def classify(start, end, length, ter):
    """Spanning if the two breakpoints are on different arms; arms = [0, ter) and [ter, L)."""
    a, b = np.asarray(start, float) % length, np.asarray(end, float) % length
    return np.where((a < ter) != (b < ter), "spanning", "local")


def summarise(ev, sp):
    """Per species: number of families, mean and within SS of log10 size (kb), log10 genome size (Mb)."""
    ev = ev.assign(y=np.log10(ev["size_bp"].astype(float) / 1000))
    g = ev.groupby("species_id")["y"]
    d = pd.DataFrame({"n": g.size(), "ybar": g.mean(),
                      "ssw": g.apply(lambda v: float(((v - v.mean()) ** 2).sum()))})
    d = d.join(sp.set_index("species_id")[[GENOME]], how="inner")
    d["x"] = np.log10(d[GENOME].astype(float) / 1e6)
    return d


def pmm(d, tree, nulls=(0.0, 1.0)):
    d = d[d.index.isin(tree.pos)]
    if len(d) < 5:
        return None
    r = fit_pmm(d["ybar"].to_numpy(), d["n"].to_numpy(float), float(d["ssw"].sum()),
                d["x"].to_numpy(), tree.sub(d.index), list(nulls), phylogenetic=True)
    out = {"n_species": int(r["n_species"]), "n_families": int(r["n_events"]), "slope": r["slope"],
           "ci_low": r["slope_ci_low"], "ci_high": r["slope_ci_high"]}
    for s in nulls:
        out[f"p_vs_{s:g}"] = r[f"p_slope_vs_{s:g}"]
    return out


def fit_classes(ev, sp, tree, label, **extra):
    rows = []
    for cls in ["spanning", "local"]:
        r = pmm(summarise(ev[ev["class"] == cls], sp), tree)
        if r:
            rows.append({"check": label, "class": cls, **extra, **r,
                         "share_spanning": float((ev["class"] == "spanning").mean())})
    return rows


# ---------------------------------------------------------------------------
# 1. Terminus position
# ---------------------------------------------------------------------------

def gc_skew_terminus(fasta, window=10_000):
    seq = "".join(line.strip() for line in open(fasta) if not line.startswith(">")).upper()
    n = len(seq)
    skew = []
    for i in range(0, n, window):
        w = seq[i:i + window]
        g, c = w.count("G"), w.count("C")
        skew.append((g - c) / (g + c) if g + c else 0.0)
    cum = np.cumsum(skew)
    return int(np.argmax(cum)) * window + window // 2, n


def inferred_termini(prepared_root, cohort, species_ids):
    rows = []
    for s in species_ids:
        pkg = Path(prepared_root) / s / cohort
        g = pkg / "genomes.tsv"
        if not g.is_file():
            continue
        t = pd.read_csv(g, sep="\t", dtype=str)
        ref = t[t["is_reference"].str.lower() == "true"]
        if ref.empty:
            continue
        ter, n = gc_skew_terminus(pkg / ref["fasta_path"].iloc[0])
        rows.append({"species_id": s, "terminus_bp": ter, "length_bp": n, "terminus_fraction": ter / n})
    return pd.DataFrame(rows)


def terminus_checks(fam, sp, tree, ter_table):
    L = fam["species_id"].map(sp.set_index("species_id")["reference_length_bp"]).astype(float)
    rows = []
    base = fam.assign(**{"class": classify(fam["start"], fam["end"], L, L / 2)})
    rows += fit_classes(base, sp, tree, "as reported (terminus at half)", setting=0)
    for shift in [-0.15, -0.10, -0.05, 0.05, 0.10, 0.15]:
        ev = fam.assign(**{"class": classify(fam["start"], fam["end"], L, L * (0.5 + shift))})
        rows += fit_classes(ev, sp, tree, "terminus shifted", setting=shift)
    a, b = fam["start"] % L, fam["end"] % L
    for m in [0.025, 0.05, 0.10]:
        near = lambda p: (np.minimum(p, L - p) < m * L) | (np.abs(p - L / 2) < m * L)  # noqa: E731
        keep = ~(near(a) | near(b))
        rows += fit_classes(base[keep], sp, tree, "margin around ori/ter", setting=m)
    if ter_table is not None and len(ter_table):
        T = fam["species_id"].map(ter_table.set_index("species_id")["terminus_bp"]).astype(float)
        ok = T.notna()
        ev = fam[ok].assign(**{"class": classify(fam.loc[ok, "start"], fam.loc[ok, "end"], L[ok], T[ok])})
        rows += fit_classes(ev, sp, tree, "terminus from GC skew", setting=np.nan)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 2. Random-breakpoint null
# ---------------------------------------------------------------------------

def breakpoint_null(fam, sp, tree, replicates, floor, seed):
    rng = np.random.default_rng(seed)
    counts = fam.groupby("species_id").size()
    L = sp.set_index("species_id")["reference_length_bp"].astype(float)
    sim_rows = []
    for rep in range(replicates):
        parts = []
        for s, k in counts.items():
            length, got = L[s], []
            while sum(len(x) for x in got) < k:
                a, b = rng.uniform(0, length, (2, 4 * k))
                start, end = np.minimum(a, b), np.maximum(a, b)
                span = end - start
                flip = span > length / 2            # shorter arc: the other way round the circle
                start, end = np.where(flip, end, start), np.where(flip, start + length, end)
                size = np.minimum(span, length - span)
                keep = size >= floor
                got.append(np.column_stack([start[keep], end[keep], size[keep]]))
            arr = np.vstack(got)[:k]
            parts.append(pd.DataFrame({"species_id": s, "start": arr[:, 0], "end": arr[:, 1],
                                       "size_bp": arr[:, 2]}))
        ev = pd.concat(parts, ignore_index=True)
        Ls = ev["species_id"].map(L)
        ev["class"] = classify(ev["start"], ev["end"], Ls, Ls / 2)
        r_all = pmm(summarise(ev, sp), tree)
        row = {"replicate": rep, "share_spanning": float((ev["class"] == "spanning").mean()),
               "slope_all": r_all["slope"] if r_all else np.nan}
        for cls in ["spanning", "local"]:
            r = pmm(summarise(ev[ev["class"] == cls], sp), tree)
            row[f"slope_{cls}"] = r["slope"] if r else np.nan
        sim_rows.append(row)
    sims = pd.DataFrame(sim_rows)
    summary = []
    for col in ["slope_all", "slope_spanning", "slope_local", "share_spanning"]:
        v = sims[col].dropna()
        summary.append({"quantity": col, "mean": v.mean(), "low_2.5": v.quantile(0.025),
                        "high_97.5": v.quantile(0.975), "replicates": len(v)})
    return sims, pd.DataFrame(summary)


# ---------------------------------------------------------------------------
# 3. One species per genus / family
# ---------------------------------------------------------------------------

def thinning(fam, sp, tree, taxa, picks, seed):
    rng = np.random.default_rng(seed)
    rows, dist, n_groups = [], [], {}
    n_fam = fam.groupby("species_id").size()
    for rank in ["genus", "family"]:
        if rank not in taxa:
            continue
        groups = taxa.loc[taxa.index.isin(n_fam.index), rank].dropna()
        # (a) deterministic: species with most families per group
        best = groups.groupby(groups).apply(lambda g: n_fam.loc[g.index].idxmax()).tolist()
        rows += fit_classes(fam[fam["species_id"].isin(best)], sp, tree, f"one species per {rank} (most families)",
                            setting=len(best))
        # (b) random picks
        members = {k: v.index.tolist() for k, v in groups.groupby(groups)}
        n_groups[rank] = len(members)
        for i in range(picks):
            chosen = [m[rng.integers(len(m))] for m in members.values()]
            for r in fit_classes(fam[fam["species_id"].isin(chosen)], sp, tree, f"one species per {rank} (random)",
                                 setting=i):
                dist.append({"rank": rank, **r})
    dist = pd.DataFrame(dist)
    if len(dist):
        for (rank, cls), g in dist.groupby(["rank", "class"]):
            rows.append({"check": f"one species per {rank} (random, {picks} picks)", "class": cls,
                         "setting": n_groups[rank], "n_species": int(g["n_species"].median()),
                         "slope": g["slope"].median(), "ci_low": g["slope"].quantile(0.025),
                         "ci_high": g["slope"].quantile(0.975), "p_vs_0": g["p_vs_0"].median(),
                         "share_p_below_0.05": float((g["p_vs_0"] < 0.05).mean())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. Whole distribution
# ---------------------------------------------------------------------------

def distribution(fam, sp, tree, min_families=5):
    rows = []
    for cls in ["spanning", "local"]:
        ev = fam[fam["class"] == cls].assign(y=lambda t: np.log10(t["size_bp"] / 1000))
        g = ev.groupby("species_id")["y"]
        ok = g.size()[g.size() >= min_families].index
        x = np.log10(sp.set_index("species_id").loc[ok, GENOME].astype(float) / 1e6)
        ids = [s for s in ok if s in tree.pos]
        for q in [0.10, 0.25, 0.50, 0.75, 0.90]:
            y = g.quantile(q).loc[ids].to_numpy()
            r = fit_pgls(y, x.loc[ids].to_numpy(), tree.sub(ids), [0.0, 1.0])
            rows.append({"check": "percentile per species (PGLS)", "class": cls, "setting": q,
                         "n_species": len(ids), "slope": r["slope"], "ci_low": r["slope_ci_low"],
                         "ci_high": r["slope_ci_high"], "p_vs_0": r["p_slope_vs_0"], "p_vs_1": r["p_slope_vs_1"]})
        # largest 10% of each species removed
        cut = ev.groupby("species_id")["size_bp"].transform(lambda v: v.quantile(0.9))
        trimmed = ev[(ev["size_bp"] <= cut) | (ev.groupby("species_id")["size_bp"].transform("size") < 3)]
        r = pmm(summarise(trimmed, sp), tree)
        if r:
            rows.append({"check": "largest 10% per species removed (mixed model)", "class": cls, "setting": 0.9, **r})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure + main
# ---------------------------------------------------------------------------

def forest(table, null_summary, out):
    t = table.dropna(subset=["slope"]).copy()
    t = t[~t["check"].str.contains("random\\)$|picks", regex=True) | t["check"].str.contains("picks")]
    labels = []
    for _, r in t.iterrows():
        s = r["setting"]
        if r["check"] == "terminus shifted":
            lab = f"terminus shifted {s:+.0%}"
        elif r["check"] == "margin around ori/ter":
            lab = f"margin {s:.1%} around ori/ter"
        elif r["check"] == "percentile per species (PGLS)":
            lab = f"{int(round(s * 100))}th percentile per species"
        else:
            lab = r["check"]
        labels.append(lab)
    t["label"] = labels
    order = list(dict.fromkeys(t["label"]))
    fig, axes = plt.subplots(1, 2, figsize=(13, 0.32 * len(order) + 2), sharey=True)
    nulls = null_summary.set_index("quantity") if null_summary is not None else None
    for ax, cls, colour in [(axes[0], "spanning", "#e8693b"), (axes[1], "local", "#2f7ad6")]:
        sub = t[t["class"] == cls].set_index("label")
        for i, lab in enumerate(order):
            y = len(order) - 1 - i
            if lab in sub.index:
                r = sub.loc[lab]
                r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
                ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=colour, lw=2)
                ax.plot(r["slope"], y, "o", color=colour, ms=5)
        ax.axvline(0, color="black", lw=0.9)
        ax.axvline(1, color="black", lw=0.9, ls=":")
        if nulls is not None and f"slope_{cls}" in nulls.index:
            n = nulls.loc[f"slope_{cls}"]
            ax.axvspan(n["low_2.5"], n["high_97.5"], color="grey", alpha=0.18, lw=0)
            ax.axvline(n["mean"], color="grey", lw=1.5, ls="--")
        ax.set_title(f"{cls.capitalize()} inversions", fontsize=11)
        ax.set_xlabel("Slope (95% CI)")
        ax.grid(axis="x", alpha=0.25)
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels(order[::-1], fontsize=8)
    fig.suptitle("Robustness of the spanning vs local result  (grey band: random-breakpoint null)",
                 fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--taxonomy", type=Path, required=True, help="data/processed/species.tsv (genus, family)")
    ap.add_argument("--prepared-root", type=Path, help="for the GC-skew terminus; skipped if absent")
    ap.add_argument("--cohort", default="selected_50_v1")
    ap.add_argument("--replicates", type=int, default=100)
    ap.add_argument("--picks", type=int, default=200)
    ap.add_argument("--floor-bp", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    sp = pd.read_csv(a.loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    fam = pd.read_csv(a.loci / "families.tsv", sep="\t")
    fam = fam[fam["tier"].isin(HM) & fam["species_id"].isin(sp["species_id"]) & (fam["size_bp"] > 0)].copy()
    fam["class"] = fam["position_class"]
    tree = Tree(a.tree, sp["species_id"].tolist())
    taxa = pd.read_csv(a.taxonomy, sep="\t", dtype=str).set_index("species_id")

    print("== 1. Terminus position")
    ter = None
    if a.prepared_root and a.prepared_root.is_dir():
        ter = inferred_termini(a.prepared_root, a.cohort, sp["species_id"].tolist())
        ter.to_csv(a.out / "termini_gc_skew.tsv", sep="\t", index=False)
        f = ter["terminus_fraction"]
        print(f"GC-skew terminus for {len(ter)} species: median at {f.median():.0%} of the chromosome "
              f"(middle 90%: {f.quantile(0.05):.0%}-{f.quantile(0.95):.0%})")
    t1 = terminus_checks(fam, sp, tree, ter)

    print("== 2. Random-breakpoint null")
    sims, null_summary = breakpoint_null(fam, sp, tree, a.replicates, a.floor_bp, a.seed)
    sims.to_csv(a.out / "breakpoint_null_replicates.tsv", sep="\t", index=False)
    null_summary.to_csv(a.out / "breakpoint_null_summary.tsv", sep="\t", index=False)

    print("== 3. One species per genus / family")
    t3 = thinning(fam, sp, tree, taxa, a.picks, a.seed)

    print("== 4. Whole distribution")
    t4 = distribution(fam, sp, tree)

    table = pd.concat([t1, t3, t4], ignore_index=True)
    table.to_csv(a.out / "robustness.tsv", sep="\t", index=False)
    forest(table, null_summary, a.out / "robustness.png")

    show = table.copy()
    for c in ["p_vs_0", "p_vs_1"]:
        if c in show:
            show[c] = show[c].map(lambda v: f"{v:.2g}" if pd.notna(v) else "")
    cols = [c for c in ["check", "setting", "class", "n_species", "n_families", "slope", "ci_low", "ci_high",
                        "p_vs_0", "p_vs_1", "share_spanning", "share_p_below_0.05"] if c in show]
    with pd.option_context("display.width", 230, "display.max_columns", 20, "display.max_rows", 200):
        print(show[cols].round(3).to_string(index=False))
        print("\nRandom-breakpoint null (no biology, sizes scale with the chromosome):")
        print(null_summary.round(3).to_string(index=False))
    print(f"\nWritten to {a.out}: robustness.tsv, robustness.png, breakpoint_null_summary.tsv"
          + (", termini_gc_skew.tsv" if ter is not None else ""))


if __name__ == "__main__":
    main()
