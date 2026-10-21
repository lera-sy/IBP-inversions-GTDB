"""
Phase 3A: geometry null model. No genomes needed.

Question: what slope of log10(inversion size) on log10(chromosome length) do
the measurement rules produce on their own, with no biology?

Two rules can create a slope:
- the shorter-arc rule: on a circle an inversion is reported as the shorter of
  its two arcs, so a reported size can never exceed half the chromosome;
- the detection floor: calls below `floor_bp` are dropped.

Two nulls are simulated:
- fixed_size   inversion size does NOT depend on chromosome length: sizes come
               from one distribution for every species (several shapes tried,
               or an empirical list of sizes). Any slope here is an artefact.
               This is S0.
- random_breakpoints  both breakpoints placed uniformly on the circle, so
               size scales with the chromosome by construction. This is the
               proportional reference (expected slope close to 1).

Each event also gets the team's position class (origin at 0 after dnaapler,
terminus at half the chromosome): `local` if both breakpoints are on one
replication arm, `spanning` if they are on opposite arms. Spanning inversions
are larger by construction, so their slope is reported separately.

Slopes are fitted two ways, as in the real analysis:
- events: every event is a data point (OLS; the mixed model's fixed slope
  without phylogeny);
- means:  one geometric mean per species (min `min_events` events).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

DISTRIBUTIONS = {
    # name: (kind, parameters); sizes in bp
    "loguniform_floor_to_5Mb": ("loguniform", {"low": None, "high": 5e6}),
    "lognormal_median_30kb": ("lognormal", {"median": 3e4, "sd_log10": 0.6}),
    "lognormal_median_200kb": ("lognormal", {"median": 2e5, "sd_log10": 0.6}),
}


def draw_sizes(kind, params, n, rng, floor_bp, empirical=None):
    if kind == "loguniform":
        low = np.log10(params["low"] or floor_bp)
        return 10 ** rng.uniform(low, np.log10(params["high"]), n)
    if kind == "lognormal":
        return 10 ** rng.normal(np.log10(params["median"]), params["sd_log10"], n)
    if kind == "empirical":
        return rng.choice(empirical, size=n, replace=True)
    raise ValueError(f"Unknown distribution kind: {kind}")


def position_class(start, end, length):
    """Team definition: arms are [0, L/2) and [L/2, L); opposite arms = spanning."""
    p1, p2 = start % length, end % length
    return np.where((p1 < length / 2) != (p2 < length / 2), "spanning", "local")


def simulate_species(length, n_events, null, rng, floor_bp, dist=None, empirical=None):
    """Return reported sizes and position classes for one chromosome."""
    sizes, starts = [], []
    have = 0
    for _ in range(1000):
        if have >= n_events:
            break
        m = max(4 * (n_events - have), 64)
        start = rng.uniform(0, length, m)
        if null == "random_breakpoints":
            true = rng.uniform(0, length, m)
        else:
            kind, params = dist
            true = draw_sizes(kind, params, m, rng, floor_bp, empirical)
            true = true[true < length]          # cannot be longer than the chromosome
            start = start[: len(true)]
        reported = np.minimum(true, length - true)  # shorter-arc rule
        # The reported arc starts where the shorter arc starts.
        flipped = true > length / 2
        arc_start = np.where(flipped, start + true, start)
        keep = reported >= floor_bp
        sizes.append(reported[keep])
        starts.append(arc_start[keep])
        have += int(keep.sum())
    if have < n_events:
        raise RuntimeError("Could not draw enough events above the floor; check the size distribution.")
    size = np.concatenate(sizes)[:n_events]
    arc_start = np.concatenate(starts)[:n_events]
    return size, position_class(arc_start, arc_start + size, length)


def ols_slope(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) == 0:
        return np.nan
    return float(np.polyfit(x, y, 1)[0])


def one_replicate(lengths, n_events, null, rng, floor_bp, min_events, dist=None, empirical=None,
                  x_values=None):
    """n_events: one number for all species, or one per species.
    x_values: per-species value regressed on (default: chromosome length)."""
    counts = np.broadcast_to(np.asarray(n_events), (len(lengths),))
    xs = lengths if x_values is None else x_values
    rows = []
    for species, (length, k, xv) in enumerate(zip(lengths, counts, xs)):
        if k < 1:
            continue
        size, cls = simulate_species(length, int(k), null, rng, floor_bp, dist, empirical)
        rows.append(pd.DataFrame({"species": species, "length": xv, "size": size, "cls": cls}))
    ev = pd.concat(rows, ignore_index=True)
    ev["x"], ev["y"] = np.log10(ev["length"]), np.log10(ev["size"] / 1000)
    out = {}
    for subset in ["all", "local", "spanning"]:
        d = ev if subset == "all" else ev[ev["cls"] == subset]
        out[(subset, "events")] = ols_slope(d["x"], d["y"])
        g = d.groupby("species").agg(x=("x", "first"), y=("y", "mean"), n=("y", "size"))
        g = g[g["n"] >= min_events]
        out[(subset, "means")] = ols_slope(g["x"], g["y"])
    out[("spanning_fraction", "events")] = float((ev["cls"] == "spanning").mean())
    return out


def run(lengths_mb=None, species_table=None, n_species=80, n_events=30, replicates=200,
        floor_bp=6000, min_events=3, seed=42, empirical_sizes=None):
    rng = np.random.default_rng(seed)
    x_values, counts = None, n_events
    if species_table is not None:
        tab = pd.read_csv(species_table, sep="\t")
        lengths = tab["chromosome_length_bp"].astype(float).to_numpy()
        if "genome_length_bp" in tab:      # regress on genome size, cap on chromosome
            x_values = tab["genome_length_bp"].astype(float).to_numpy()
        if "n_events" in tab:              # real number of events per species
            counts = tab["n_events"].astype(int).to_numpy()
        sample_lengths = lambda: lengths  # noqa: E731
        source = (f"species table ({len(lengths)} species"
                  + (", real event counts" if "n_events" in tab else "")
                  + (", x = genome size" if x_values is not None else "") + ")")
    else:
        lo, hi = min(lengths_mb), max(lengths_mb)
        sample_lengths = lambda: 10 ** rng.uniform(np.log10(lo * 1e6), np.log10(hi * 1e6), n_species)  # noqa: E731
        source = f"{n_species} species, log-uniform {lo}-{hi} Mb"

    scenarios = [("random_breakpoints", "uniform_breakpoints", None)]
    scenarios += [("fixed_size", name, spec) for name, spec in DISTRIBUTIONS.items()]
    empirical = None
    if empirical_sizes is not None:
        empirical = pd.read_csv(empirical_sizes, sep="\t")["size_bp"].astype(float).to_numpy()
        scenarios.append(("fixed_size", "empirical_pooled_sizes", ("empirical", {})))

    records = []
    for null, name, spec in scenarios:
        reps = [one_replicate(sample_lengths(), counts, null, rng, floor_bp, min_events, spec, empirical,
                              x_values) for _ in range(replicates)]
        for key in reps[0]:
            vals = np.array([r[key] for r in reps], float)
            vals = vals[~np.isnan(vals)]
            if len(vals) == 0:      # e.g. too few species with >= min_events in a subset
                stats_ = (np.nan, np.nan, np.nan)
            else:
                stats_ = (round(float(vals.mean()), 3), round(float(np.percentile(vals, 2.5)), 3),
                          round(float(np.percentile(vals, 97.5)), 3))
            records.append({"null": null, "size_distribution": name, "subset": key[0], "fit": key[1],
                            "slope_mean": stats_[0], "slope_2.5": stats_[1], "slope_97.5": stats_[2],
                            "replicates": len(vals)})
    table = pd.DataFrame(records)
    table.attrs["source"] = source
    return table


def main():
    p = argparse.ArgumentParser(description="Phase 3A geometry null model")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--lengths-mb", type=float, nargs="+", default=[1, 10])
    p.add_argument("--species-table", type=Path, help="TSV with chromosome_length_bp per species")
    p.add_argument("--empirical-sizes", type=Path, help="TSV with size_bp (pooled real sizes)")
    p.add_argument("--n-species", type=int, default=80)
    p.add_argument("--n-events", type=int, default=30)
    p.add_argument("--replicates", type=int, default=200)
    p.add_argument("--floor-bp", type=int, default=6000)
    p.add_argument("--min-events", type=int, default=3)
    p.add_argument("--seed", type=int, default=42)
    a = p.parse_args()
    table = run(a.lengths_mb, a.species_table, a.n_species, a.n_events, a.replicates,
                a.floor_bp, a.min_events, a.seed, a.empirical_sizes)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(a.out, sep="\t", index=False)
    print(f"Lengths: {table.attrs['source']}; floor {a.floor_bp} bp; {a.replicates} replicates\n")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
