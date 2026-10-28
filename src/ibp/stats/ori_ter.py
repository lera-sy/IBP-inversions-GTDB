"""
Spanning inversions split by what they cross: the origin (ori) or the terminus (ter).

Biological reason to split: the two ends of the replication arms are different
places. Around ori sit the highly expressed, dosage-sensitive genes (rRNA,
ribosomal proteins; more copies during fast growth). Around ter sit replication
fork termination, chromosome-dimer resolution (dif/XerCD/FtsK) and a gene-poor,
AT-rich, horizontally acquired region. If the size scaling and the symmetry are
a property of replication-arm balance itself, both groups should behave alike;
if one end drives it, the groups differ.

Per group (ori-crossing, ter-crossing, local), high + medium families:
  - number, species, median size;
  - frequency vs chance given each family's own size (random placement of the
    same arc), observed / expected with species-resampled CI;
  - symmetry: share with asymmetry A < 0.2 (A = |d1 - d2| / (d1 + d2)) vs chance;
  - size slope (phylogenetic mixed model) and the room-alone null slope
    (real sizes shuffled across species, same rules, same split), empirical p;
  - ori vs ter: one joint model, slope difference with CI and likelihood ratio.
Done with the terminus at half and, if given, the GC-skew terminus.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ibp.stats.class_tests import fit_joint, simulate
from ibp.stats.robustness import HM, Tree, pmm, summarise

GROUPS = ["ori", "ter", "local"]


def axis_geometry(start, end, L, T):
    start, end, L, T = (np.asarray(v, float) for v in (start, end, L, T))
    s = (end - start) % L
    o0, oT = (0 - start) % L, (T - start) % L
    in0, inT = (o0 > 0) & (o0 < s), (oT > 0) & (oT < s)
    grp = np.where(in0 & ~inT, "ori", np.where(inT & ~in0, "ter", np.where(in0 & inT, "both", "local")))
    o = np.where(in0, o0, oT)
    A = np.where(in0 ^ inT, np.abs(2 * o - s) / np.where(s > 0, s, 1), np.nan)
    return s, grp, A


def expected(s, L, T, rng, draws=400):
    n = len(s)
    st = rng.uniform(0, 1, (n, draws)) * L[:, None]
    _, g, A = axis_geometry(st, (st + s[:, None]) % L[:, None], L[:, None], T[:, None])
    return {k: (g == k).mean(1) for k in ["ori", "ter"]}, {k: ((g == k) & (A < 0.2)).mean(1) for k in ["ori", "ter"]}


def boot_ratio(df, num, den, rng, B=2000):
    g = df.groupby("species_id")[[num, den]].sum()
    a, b = g[num].to_numpy(), g[den].to_numpy()
    idx = rng.integers(0, len(g), (B, len(g)))
    r = a[idx].sum(1) / b[idx].sum(1)
    return np.quantile(r, 0.025), np.quantile(r, 0.975)


def joint(ev, sp, tree):
    rows = []
    for c, grp in enumerate(["ori", "ter"]):
        d = summarise(ev[ev["grp"] == grp], sp)
        d = d[d.index.isin(tree.pos)]
        for s_, r in d.iterrows():
            rows.append({"species_id": s_, "c": c, "n": r["n"], "ybar": r["ybar"], "ssw": r["ssw"], "x": r["x"]})
    cl = pd.DataFrame(rows)
    full, red = fit_joint(cl, tree), fit_joint(cl, tree, common_slope=True)
    b, V = full["beta"], full["cov"]
    d = b[2] - b[3]
    se = float(np.sqrt(V[2, 2] + V[3, 3] - 2 * V[2, 3]))
    t = stats.t.ppf(0.975, full["n_species"] - 4)
    lr = 2 * (full["loglik"] - red["loglik"])
    return {"slope_ori_joint": b[2], "slope_ter_joint": b[3], "difference": d, "ci_low": d - t * se,
            "ci_high": d + t * se, "p_wald": 2 * stats.t.sf(abs(d / se), full["n_species"] - 4),
            "p_LRT": float(stats.chi2.sf(max(lr, 0), 1)), "n_species": full["n_species"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--termini", type=Path)
    ap.add_argument("--replicates", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    sp = pd.read_csv(a.loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    fam = pd.read_csv(a.loci / "families.tsv", sep="\t")
    fam = fam[fam["tier"].isin(HM) & fam["species_id"].isin(sp["species_id"]) & (fam["size_bp"] > 0)].reset_index(drop=True)
    tree = Tree(a.tree, sp["species_id"].tolist())
    Lmap = sp.set_index("species_id")["reference_length_bp"].astype(float)
    L = fam["species_id"].map(Lmap).to_numpy()
    settings = [("terminus at half", L / 2)]
    if a.termini and Path(a.termini).is_file():
        ter = pd.read_csv(a.termini, sep="\t").set_index("species_id")["terminus_bp"]
        T = fam["species_id"].map(ter).astype(float).to_numpy()
        settings.append(("terminus from GC skew", np.where(np.isnan(T), L / 2, T)))

    rows, joints = [], []
    for label, T in settings:
        s, grp, A = axis_geometry(fam["start"], fam["end"], L, T)
        ev = fam.assign(grp=grp, A=A)
        e_grp, e_sym = expected(s, L, T, rng)
        per = pd.DataFrame({"species_id": fam["species_id"]})
        for k in ["ori", "ter"]:
            per[f"o_{k}"], per[f"e_{k}"] = (grp == k).astype(float), e_grp[k]
            per[f"os_{k}"], per[f"es_{k}"] = ((grp == k) & (A < 0.2)).astype(float), e_sym[k]
        for g in GROUPS:
            sub = ev[ev["grp"] == g]
            r = pmm(summarise(sub, sp), tree)
            row = {"terminus": label, "group": g, "families": len(sub), "species": sub["species_id"].nunique(),
                   "median_kb": float(np.median(sub["size_bp"]) / 1000) if len(sub) else np.nan}
            if g != "local":
                o, e = per[f"o_{g}"].sum(), per[f"e_{g}"].sum()
                lo, hi = boot_ratio(per, f"o_{g}", f"e_{g}", rng)
                os_, es_ = per[f"os_{g}"].sum(), per[f"es_{g}"].sum()
                slo, shi = boot_ratio(per, f"os_{g}", f"es_{g}", rng)
                row.update({"obs_over_exp": o / e, "oe_ci": f"{lo:.2f}-{hi:.2f}", "share_symmetric_A<0.2": float((sub["A"] < 0.2).mean()),
                            "symmetric_obs_over_exp": os_ / es_, "sym_ci": f"{slo:.2f}-{shi:.2f}"})
            if r:
                row.update({"slope": r["slope"], "ci_low": r["ci_low"], "ci_high": r["ci_high"], "p_vs_0": r["p_vs_0"]})
            rows.append(row)
        joints.append({"terminus": label, **joint(ev, sp, tree)})

    # room-alone null for the three groups (terminus at half, as in the main analysis)
    counts = fam.groupby("species_id").size()
    pooled = fam["size_bp"].astype(float).to_numpy()
    reps = []
    for i in range(a.replicates):
        sim = simulate(counts, Lmap, pooled, "shuffle", rng, 6000, 5000)
        Ls = sim["species_id"].map(Lmap).to_numpy()
        _, g, _ = axis_geometry(sim["start"], sim["end"], Ls, Ls / 2)
        sim = sim.assign(grp=g)
        reps.append({k: (pmm(summarise(sim[sim["grp"] == k], sp), tree) or {"slope": np.nan})["slope"] for k in GROUPS})
        if (i + 1) % 50 == 0:
            print(f"   null {i + 1}/{a.replicates}", flush=True)
    reps = pd.DataFrame(reps)
    res = pd.DataFrame(rows)
    for k in GROUPS:
        v = reps[k].dropna().to_numpy()
        m = (res["terminus"] == "terminus at half") & (res["group"] == k)
        obs = res.loc[m, "slope"].iloc[0]
        p = min(1.0, 2 * min((1 + (v <= obs).sum()) / (len(v) + 1), (1 + (v >= obs).sum()) / (len(v) + 1)))
        res.loc[m, "room_alone"] = v.mean()
        res.loc[m, "room_range"] = f"{np.quantile(v, 0.025):.2f}-{np.quantile(v, 0.975):.2f}"
        res.loc[m, "p_vs_room"] = p
    res.to_csv(a.out / "ori_ter.tsv", sep="\t", index=False)
    pd.DataFrame(joints).to_csv(a.out / "ori_vs_ter_test.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(res.round(3).to_string(index=False))
    print("\nori vs ter, one joint model:")
    print(pd.DataFrame(joints).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
