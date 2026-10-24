"""
Spanning vs local: formal tests (reviewer feedback round 2, 24 Oct 2026).

1. Interaction test. One phylogenetic mixed model fitted on both classes at
   once:  log size ~ class + class x log genome size, with
     - a tree-shared species effect per class (2x2 covariance: each class its
       own variance plus a correlation between the classes),
     - a species-unique effect per class (also 2x2),
     - within-species variance per class.
   Reported: slope per class, slope difference (spanning - local) with 95% CI
   (Wald), and a likelihood-ratio test against a model with one common slope.

2. Empirical p-values against two nulls, with the SAME model and the SAME
   classification as the observed data. Per species: the real number of
   families, real chromosome length, real genome size; breakpoints at random;
   shorter-arc rule; 6 kb floor; events with both breakpoints within the
   family tolerance (5 kb) of an earlier event are merged (as families are);
   terminus at half.
     shuffle            sizes drawn from the pooled real family sizes (room only)
     random_breakpoints both breakpoints uniform (size scales with the chromosome)
   For each quantity (slope all / spanning / local, slope difference, spanning
   share): null mean, 95% range, empirical two-sided p.

3. Spanning share and symmetry.
   - Share of families that span, with uncertainty from resampling species and
     genera (families of one species are not independent), and a tree-aware
     mean of per-species shares.
   - Expected share GIVEN EACH FAMILY'S SIZE: an arc of length s placed at a
     random position spans ori or ter with probability 2s/L (larger inversions
     span more often by chance). Observed / expected > 1 means enrichment
     beyond what the sizes explain.
   - Symmetry of spanning families: A = |d1 - d2| / (d1 + d2), with d1, d2
     the distances from the crossed axis (ori or ter) to each breakpoint.
     A = 0 perfectly symmetric, A = 1 one breakpoint on the axis. Under random
     placement of the same arc A is close to uniform. Observed / expected per
     A bin shows whether the enrichment sits in the symmetric ones. Done with
     the terminus at half and, if given, at the GC-skew terminus.
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import optimize, stats

from ibp.analysis_legacy.final_inversion_model import gls, lambda_matrix
from ibp.stats.robustness import HM, Tree, classify, pmm, summarise

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CLASSES = ["spanning", "local"]
BINS = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0001])


# ---------------------------------------------------------------------------
# 1. Joint (bivariate) phylogenetic mixed model
# ---------------------------------------------------------------------------

def cells(ev, sp, tree):
    """One row per (species, class) with n, mean, within SS, x."""
    rows = []
    for c_i, cls in enumerate(CLASSES):
        d = summarise(ev[ev["class"] == cls], sp)
        d = d[d.index.isin(tree.pos)]
        for s, r in d.iterrows():
            rows.append({"species_id": s, "c": c_i, "n": r["n"], "ybar": r["ybar"], "ssw": r["ssw"], "x": r["x"]})
    return pd.DataFrame(rows)


def _cov2(log_v1, log_v2, z):
    v1, v2, r = np.exp(log_v1), np.exp(log_v2), np.tanh(z)
    c = r * np.sqrt(v1 * v2)
    return np.array([[v1, c], [c, v2]])


def fit_joint(cl, tree, common_slope=False):
    """ML fit of the bivariate PMM on cell summaries. Returns dict."""
    ids = list(dict.fromkeys(cl["species_id"]))
    idx = np.array([ids.index(s) for s in cl["species_id"]])
    C = tree.sub(ids)
    Cn = C / np.mean(np.diag(C))
    Ccell = Cn[np.ix_(idx, idx)]
    same = (idx[:, None] == idx[None, :]).astype(float)
    c = cl["c"].to_numpy()
    n = cl["n"].to_numpy(float)
    y = cl["ybar"].to_numpy()
    x = cl["x"].to_numpy()
    is_s, is_l = (c == 0).astype(float), (c == 1).astype(float)
    if common_slope:
        X = np.column_stack([is_s, is_l, x])
    else:
        X = np.column_stack([is_s, is_l, is_s * x, is_l * x])
    m = len(y)
    N = [n[c == k].sum() for k in (0, 1)]
    K = [(c == k).sum() for k in (0, 1)]
    SSW = [cl.loc[c == k, "ssw"].sum() for k in (0, 1)]

    def build(theta):
        St = _cov2(*theta[0:3])
        Ss = _cov2(*theta[3:6])
        e = np.exp(theta[6:8])
        V = St[c][:, c] * Ccell + Ss[c][:, c] * same + np.diag(e[c] / n)
        return V, e

    def nll(theta, want=False):
        V, e = build(theta)
        try:
            L = np.linalg.cholesky(V)
        except np.linalg.LinAlgError:
            return 1e12
        Xw = np.linalg.solve(L, X)
        yw = np.linalg.solve(L, y)
        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        r = yw - Xw @ beta
        ll = -0.5 * (m * np.log(2 * np.pi) + 2 * np.log(np.diag(L)).sum() + r @ r)
        for k in (0, 1):
            if N[k] > K[k]:
                ll += -0.5 * ((N[k] - K[k]) * np.log(2 * np.pi * e[k]) + SSW[k] / e[k])
        if want:
            return beta, np.linalg.inv(Xw.T @ Xw), ll
        return -ll

    v = [np.var(y[c == k]) + 1e-3 for k in (0, 1)]
    w = [SSW[k] / max(N[k] - K[k], 1) + 1e-3 for k in (0, 1)]
    base = np.log([v[0] / 2, v[1] / 2, 1, v[0] / 2, v[1] / 2, 1, w[0], w[1]])
    base[2] = base[5] = 0.0
    bounds = [(-20, 4), (-20, 4), (-3, 3)] * 2 + [(-20, 4)] * 2
    best = None
    for shift in (0.0, -2.0, 1.0):
        start = base + np.array([shift, shift, 0, shift, shift, 0, 0, 0])
        res = optimize.minimize(nll, start, method="L-BFGS-B", bounds=bounds)
        if best is None or res.fun < best.fun:
            best = res
    beta, cov, ll = nll(best.x, want=True)
    out = {"loglik": ll, "beta": beta, "cov": cov, "n_species": len(ids), "n_cells": m,
           "corr_tree": float(np.tanh(best.x[2])), "corr_species": float(np.tanh(best.x[5]))}
    return out


def interaction_test(ev, sp, tree):
    cl = cells(ev, sp, tree)
    full = fit_joint(cl, tree)
    red = fit_joint(cl, tree, common_slope=True)
    b, V = full["beta"], full["cov"]
    se = np.sqrt(np.diag(V))
    df = full["n_species"] - 4
    t = stats.t.ppf(0.975, df)
    d = b[2] - b[3]
    se_d = float(np.sqrt(V[2, 2] + V[3, 3] - 2 * V[2, 3]))
    lr = 2 * (full["loglik"] - red["loglik"])
    rows = []
    for k, cls in enumerate(CLASSES):
        s, e = b[2 + k], se[2 + k]
        rows.append({"quantity": f"slope {cls} (joint model)", "estimate": s, "ci_low": s - t * e,
                     "ci_high": s + t * e, "p_vs_0": 2 * stats.t.sf(abs(s / e), df)})
    rows.append({"quantity": "slope difference spanning - local", "estimate": d, "ci_low": d - t * se_d,
                 "ci_high": d + t * se_d, "p_vs_0": 2 * stats.t.sf(abs(d / se_d), df),
                 "p_LRT": float(stats.chi2.sf(max(lr, 0), 1))})
    rows.append({"quantity": "common slope (reduced model)", "estimate": red["beta"][2]})
    rows.append({"quantity": "between-class correlation of species effects (tree / unique)",
                 "estimate": full["corr_tree"], "ci_low": np.nan, "ci_high": full["corr_species"]})
    meta = {"n_species": full["n_species"], "n_cells": full["n_cells"]}
    return pd.DataFrame(rows), meta, d, se_d


def joint_difference(ev, sp, tree):
    cl = cells(ev, sp, tree)
    if cl["c"].nunique() < 2 or (cl["c"] == 1).sum() < 5 or (cl["c"] == 0).sum() < 5:
        return np.nan
    f = fit_joint(cl, tree)
    return float(f["beta"][2] - f["beta"][3])


# ---------------------------------------------------------------------------
# 2. Nulls with the same model and the same procedure
# ---------------------------------------------------------------------------

def circ(p, q, L):
    d = np.abs(p - q) % L
    return np.minimum(d, L - d)


def simulate(counts, L, pooled, null, rng, floor, tol):
    parts = []
    for s, k in counts.items():
        length = L[s]
        st, en, sz = [], [], []
        tries = 0
        while len(sz) < k and tries < 200:
            tries += 1
            m = 4 * k
            a = rng.uniform(0, length, m)
            if null == "random_breakpoints":
                true = rng.uniform(0, length, m)
            else:
                true = rng.choice(pooled, m)
                keep = true < length
                a, true = a[keep], true[keep]
            size = np.minimum(true, length - true)
            start = np.where(true > length / 2, (a + true) % length, a)
            end = (start + size) % length
            for i in np.flatnonzero(size >= floor):
                if len(sz) >= k:
                    break
                if st:
                    S, E = np.array(st), np.array(en)
                    dup = ((circ(S, start[i], length) <= tol) & (circ(E, end[i], length) <= tol)) | \
                          ((circ(S, end[i], length) <= tol) & (circ(E, start[i], length) <= tol))
                    if dup.any():
                        continue
                st.append(start[i])
                en.append(end[i])
                sz.append(size[i])
        parts.append(pd.DataFrame({"species_id": s, "start": st, "end": en, "size_bp": sz}))
    ev = pd.concat(parts, ignore_index=True)
    Ls = ev["species_id"].map(L)
    ev["class"] = classify(ev["start"], ev["end"], Ls, Ls / 2)
    return ev


def quantities(ev, sp, tree, with_joint=True):
    out = {"share_spanning": float((ev["class"] == "spanning").mean())}
    r = pmm(summarise(ev, sp), tree)
    out["slope_all"] = r["slope"] if r else np.nan
    for cls in CLASSES:
        r = pmm(summarise(ev[ev["class"] == cls], sp), tree)
        out[f"slope_{cls}"] = r["slope"] if r else np.nan
    out["slope_difference"] = joint_difference(ev, sp, tree) if with_joint else np.nan
    return out


def null_tests(fam, sp, tree, observed, replicates, floor, tol, seed):
    rng = np.random.default_rng(seed)
    counts = fam.groupby("species_id").size()
    L = sp.set_index("species_id")["reference_length_bp"].astype(float)
    pooled = fam["size_bp"].astype(float).to_numpy()
    reps = []
    for null in ["shuffle", "random_breakpoints"]:
        for i in range(replicates):
            ev = simulate(counts, L, pooled, null, rng, floor, tol)
            reps.append({"null": null, "replicate": i, **quantities(ev, sp, tree)})
            if (i + 1) % 25 == 0:
                print(f"   {null}: {i + 1}/{replicates}", flush=True)
    reps = pd.DataFrame(reps)
    rows = []
    for null, g in reps.groupby("null", sort=False):
        for q, obs in observed.items():
            v = g[q].dropna().to_numpy()
            lo = (1 + (v <= obs).sum()) / (len(v) + 1)
            hi = (1 + (v >= obs).sum()) / (len(v) + 1)
            rows.append({"null": null, "quantity": q, "observed": obs, "null_mean": v.mean(),
                         "null_2.5": np.quantile(v, 0.025), "null_97.5": np.quantile(v, 0.975),
                         "p_empirical_two_sided": min(1.0, 2 * min(lo, hi)),
                         "observed_is": "below" if obs < v.mean() else "above", "replicates": len(v)})
    return pd.DataFrame(rows), reps


# ---------------------------------------------------------------------------
# 3. Spanning share and symmetry
# ---------------------------------------------------------------------------

def arc_geometry(start, end, L, T):
    """Arc length, class and asymmetry for arcs start -> end (forward) with axes at 0 and T."""
    start, end, L, T = (np.asarray(v, float) for v in (start, end, L, T))
    s = (end - start) % L
    o0 = (0 - start) % L          # offset of the origin from start along the arc
    oT = (T - start) % L
    in0, inT = (o0 > 0) & (o0 < s), (oT > 0) & (oT < s)
    span = in0 ^ inT
    o = np.where(in0, o0, oT)
    A = np.where(span, np.abs(2 * o - s) / np.where(s > 0, s, 1), np.nan)
    return s, span, A


def expected_counts(s, L, T, rng, draws=400):
    """Per family: P(spanning) and P(spanning and A in bin) for the same arc at a random position."""
    n = len(s)
    start = rng.uniform(0, 1, (n, draws)) * L[:, None]
    end = (start + s[:, None]) % L[:, None]
    _, span, A = arc_geometry(start, end, L[:, None], T[:, None])
    p_span = span.mean(axis=1)
    p_bins = np.stack([(span & (A >= BINS[k]) & (A < BINS[k + 1])).mean(axis=1) for k in range(len(BINS) - 1)], 1)
    return p_span, p_bins


def pgls_mean(y, C):
    """Tree-aware mean (intercept-only GLS, Pagel's lambda by ML)."""
    X = np.ones((len(y), 1))
    res = optimize.minimize_scalar(lambda lam: -gls(y, X, lambda_matrix(C, lam))[3], bounds=(0, 1), method="bounded")
    lam = float(res.x)
    for edge in (0.0, 1.0):
        if gls(y, X, lambda_matrix(C, edge))[3] > gls(y, X, lambda_matrix(C, lam))[3]:
            lam = edge
    beta, cov, rss, _ = gls(y, X, lambda_matrix(C, lam))
    se = float(np.sqrt(cov[0, 0] * rss / (len(y) - 1)))
    t = stats.t.ppf(0.975, len(y) - 1)
    return beta[0], beta[0] - t * se, beta[0] + t * se, 2 * stats.t.sf(abs(beta[0] / se), len(y) - 1), lam


def cluster_ci(per, num, den, groups, rng, B=2000):
    """Ratio sum(num)/sum(den) with CI from resampling clusters (species or genera)."""
    g = per.groupby(groups)[[num, den]].sum()
    a, b = g[num].to_numpy(), g[den].to_numpy()
    idx = rng.integers(0, len(g), (B, len(g)))
    r = a[idx].sum(1) / b[idx].sum(1)
    return np.quantile(r, 0.025), np.quantile(r, 0.975)


def share_and_symmetry(fam, sp, tree, taxa, ter_table, seed, min_families=3):
    rng = np.random.default_rng(seed)
    L = fam["species_id"].map(sp.set_index("species_id")["reference_length_bp"]).astype(float).to_numpy()
    settings = [("terminus at half", L / 2)]
    if ter_table is not None and len(ter_table):
        T = fam["species_id"].map(ter_table.set_index("species_id")["terminus_bp"]).astype(float).to_numpy()
        settings.append(("terminus from GC skew", np.where(np.isnan(T), L / 2, T)))
    share_rows, sym_rows, per_species = [], [], {}
    genus = fam["species_id"].map(taxa["genus"]).fillna(fam["species_id"])
    for label, T in settings:
        s, span, A = arc_geometry(fam["start"], fam["end"], L, T)
        p_span, p_bins = expected_counts(s, L, T, rng)
        per = pd.DataFrame({"species_id": fam["species_id"].to_numpy(), "genus": genus.to_numpy(),
                            "one": 1.0, "obs": span.astype(float), "exp": p_span})
        for k in range(len(BINS) - 1):
            per[f"obs{k}"] = (span & (A >= BINS[k]) & (A < BINS[k + 1])).astype(float)
            per[f"exp{k}"] = p_bins[:, k]
        per["A"] = A
        per_species[label] = per

        # share and observed/expected
        obs, exp = per["obs"].sum(), per["exp"].sum()
        n = len(per)
        sp_lo, sp_hi = cluster_ci(per, "obs", "one", "species_id", rng)
        ge_lo, ge_hi = cluster_ci(per, "obs", "one", "genus", rng)
        oe_lo, oe_hi = cluster_ci(per, "obs", "exp", "species_id", rng)
        oeg_lo, oeg_hi = cluster_ci(per, "obs", "exp", "genus", rng)
        g = per.groupby("species_id").agg(n=("one", "sum"), o=("obs", "sum"), e=("exp", "sum"))
        g = g[(g["n"] >= min_families) & g.index.isin(tree.pos)]
        m, lo, hi, p, lam = pgls_mean(((g["o"] - g["e"]) / g["n"]).to_numpy(), tree.sub(g.index))
        m2, lo2, hi2, _, _ = pgls_mean((g["o"] / g["n"]).to_numpy(), tree.sub(g.index))
        share_rows += [
            {"terminus": label, "quantity": "share spanning (pooled)", "estimate": obs / n,
             "ci_species": f"{sp_lo:.3f}-{sp_hi:.3f}", "ci_genus": f"{ge_lo:.3f}-{ge_hi:.3f}"},
            {"terminus": label, "quantity": f"share spanning, tree-aware mean of species ({len(g)} sp.)",
             "estimate": m2, "ci_species": f"{lo2:.3f}-{hi2:.3f}"},
            {"terminus": label, "quantity": "expected share given each family's size", "estimate": exp / n},
            {"terminus": label, "quantity": "observed / expected (pooled)", "estimate": obs / exp,
             "ci_species": f"{oe_lo:.3f}-{oe_hi:.3f}", "ci_genus": f"{oeg_lo:.3f}-{oeg_hi:.3f}"},
            {"terminus": label, "quantity": "observed - expected share, tree-aware mean of species",
             "estimate": m, "ci_species": f"{lo:.3f}-{hi:.3f}", "p": p, "lambda": lam},
        ]

        # symmetry
        for k in range(len(BINS) - 1):
            o, e = per[f"obs{k}"].sum(), per[f"exp{k}"].sum()
            lo, hi = cluster_ci(per, f"obs{k}", f"exp{k}", "species_id", rng)
            glo, ghi = cluster_ci(per, f"obs{k}", f"exp{k}", "genus", rng)
            sym_rows.append({"terminus": label, "asymmetry_bin": f"{BINS[k]:.1f}-{min(BINS[k + 1], 1):.1f}",
                             "observed": int(o), "expected": round(e, 1), "obs_over_exp": o / e,
                             "ci_low_species": lo, "ci_high_species": hi, "ci_low_genus": glo, "ci_high_genus": ghi})
        # mean asymmetry vs its null (null mean = expected A among spanning draws, about 0.5)
        sp_ev = per[per["obs"] == 1]
        null_mean = float((p_bins * (BINS[:-1] + 0.1)).sum() / p_bins.sum())
        gm = sp_ev.groupby("species_id").agg(n=("one", "sum"), A=("A", "mean"))
        gm = gm[(gm["n"] >= min_families) & gm.index.isin(tree.pos)]
        m, lo, hi, p, lam = pgls_mean(gm["A"].to_numpy() - 0.5, tree.sub(gm.index))
        a_lo, a_hi = cluster_ci(sp_ev.assign(Aw=sp_ev["A"]), "Aw", "one", "species_id", rng)
        sym_rows.append({"terminus": label, "asymmetry_bin": "mean A (pooled)", "observed": round(sp_ev["A"].mean(), 3),
                         "expected": round(null_mean, 3), "ci_low_species": a_lo, "ci_high_species": a_hi})
        sym_rows.append({"terminus": label, "asymmetry_bin": f"mean A - 0.5, tree-aware ({len(gm)} sp.)",
                         "observed": round(m, 3), "expected": 0.0, "ci_low_species": lo, "ci_high_species": hi,
                         "p": p})
    return pd.DataFrame(share_rows), pd.DataFrame(sym_rows), per_species


def symmetric_split(per, fam, sp, tree):
    """Size scaling for symmetric (A < 0.5) vs asymmetric spanning families."""
    rows = []
    ev = fam.assign(A=per["A"].to_numpy(), spanning=per["obs"].to_numpy() == 1)
    for label, keep in [("spanning, symmetric (A < 0.5)", ev["spanning"] & (ev["A"] < 0.5)),
                        ("spanning, asymmetric (A >= 0.5)", ev["spanning"] & (ev["A"] >= 0.5))]:
        r = pmm(summarise(ev[keep], sp), tree)
        if r:
            rows.append({"group": label, **r})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def figure(reps, observed, sym, out):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    cols = {"shuffle": "#e08a1e", "random_breakpoints": "#7f7f7f"}
    names = {"shuffle": "shuffle null (room only)", "random_breakpoints": "random-breakpoint null"}
    for ax, (q, title) in zip(axes.flat[:3], [("slope_spanning", "Spanning slope"), ("slope_local", "Local slope"),
                                              ("slope_difference", "Slope difference (spanning - local)")]):
        for null, g in reps.groupby("null", sort=False):
            ax.hist(g[q].dropna(), bins=30, color=cols[null], alpha=0.55, label=names[null])
        ax.axvline(observed[q], color="black", lw=2.2, label=f"observed {observed[q]:.2f}")
        ax.axvline(0, color="black", lw=0.7, ls=":")
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlabel("slope" if "difference" not in q else "difference in slope")
        ax.legend(fontsize=7, frameon=False)
    ax = axes.flat[3]
    s = sym[(sym["terminus"] == "terminus at half") & sym["asymmetry_bin"].str.match(r"^\d")]
    x = np.arange(len(s))
    ax.bar(x, s["obs_over_exp"], color="#e8693b", alpha=0.75)
    ax.errorbar(x, s["obs_over_exp"], yerr=[s["obs_over_exp"] - s["ci_low_species"],
                                           s["ci_high_species"] - s["obs_over_exp"]], fmt="none", color="black", lw=1)
    ax.axhline(1, color="black", lw=0.9, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(s["asymmetry_bin"], fontsize=8)
    ax.set_xlabel("Asymmetry A (0 = breakpoints equally far from ori/ter)")
    ax.set_ylabel("observed / expected (same sizes, random position)")
    ax.set_title("Spanning families by symmetry (95% CI: species resampled)", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--taxonomy", type=Path, required=True)
    ap.add_argument("--termini", type=Path, help="termini_gc_skew.tsv from robustness.py")
    ap.add_argument("--replicates", type=int, default=200)
    ap.add_argument("--floor-bp", type=int, default=6000)
    ap.add_argument("--tolerance-bp", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    sp = pd.read_csv(a.loci / "species.tsv", sep="\t")
    sp = sp[sp["excluded"].fillna("") == ""]
    fam = pd.read_csv(a.loci / "families.tsv", sep="\t")
    fam = fam[fam["tier"].isin(HM) & fam["species_id"].isin(sp["species_id"]) & (fam["size_bp"] > 0)].copy()
    fam = fam.reset_index(drop=True)
    L = fam["species_id"].map(sp.set_index("species_id")["reference_length_bp"]).astype(float)
    fam["class"] = classify(fam["start"], fam["end"], L, L / 2)
    tree = Tree(a.tree, sp["species_id"].tolist())
    taxa = pd.read_csv(a.taxonomy, sep="\t", dtype=str).set_index("species_id")
    ter = pd.read_csv(a.termini, sep="\t") if a.termini and a.termini.is_file() else None
    pd.set_option("display.width", 230)
    pd.set_option("display.max_columns", 20)

    print("== 1. Interaction test: one phylogenetic mixed model, both classes")
    inter, meta, d_obs, _ = interaction_test(fam, sp, tree)
    inter.to_csv(a.out / "interaction_test.tsv", sep="\t", index=False)
    print(f"{meta['n_species']} species, {meta['n_cells']} species x class cells")
    print(inter.round(4).to_string(index=False))

    print("\n== 2. Observed vs null distributions (same model, same procedure)")
    observed = quantities(fam, sp, tree, with_joint=False)
    observed["slope_difference"] = d_obs
    nt, reps = null_tests(fam, sp, tree, observed, a.replicates, a.floor_bp, a.tolerance_bp, a.seed)
    nt.to_csv(a.out / "null_tests.tsv", sep="\t", index=False)
    reps.to_csv(a.out / "null_replicates.tsv", sep="\t", index=False)
    print(nt.round(4).to_string(index=False))

    print("\n== 3. Spanning share (clustered uncertainty) and symmetry")
    share, sym, per = share_and_symmetry(fam, sp, tree, taxa, ter, a.seed)
    share.to_csv(a.out / "spanning_share.tsv", sep="\t", index=False)
    sym.to_csv(a.out / "symmetry.tsv", sep="\t", index=False)
    print(share.round(4).to_string(index=False))
    print()
    print(sym.round(3).to_string(index=False))
    split = symmetric_split(per["terminus at half"], fam, sp, tree)
    split.to_csv(a.out / "symmetric_split.tsv", sep="\t", index=False)
    print("\nSize scaling, symmetric vs asymmetric spanning families (mixed model):")
    print(split.round(4).to_string(index=False))

    figure(reps, observed, sym, a.out / "class_tests.png")
    print(f"\nWritten to {a.out}: interaction_test.tsv, null_tests.tsv, null_replicates.tsv, "
          "spanning_share.tsv, symmetry.tsv, symmetric_split.tsv, class_tests.png")


if __name__ == "__main__":
    main()
