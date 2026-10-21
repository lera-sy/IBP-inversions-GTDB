# Results: selected_50_v1 (21 Oct 2026)

79 GTDB species (min 50 genomes), 7 phyla, 1-7 Mb. 676 high+medium inversion
families (counted once per species). Tiers provisional (no raw-alignment check yet).
Slope = log10(family size, kb) on log10(genome size, Mb).

## Main question: does inversion size scale with genome size beyond geometry?

| Model | Species | Slope [95% CI] | S0 (geometry) | p vs S0 | p vs 0 | p vs 1 |
|---|---|---|---|---|---|---|
| PMM (main) | 77 | 0.54 [0.04, 1.03] | 0.53 [0.29, 0.82] | 0.98 | 0.033 | 0.067 |
| MM (no tree) | 77 | 0.32 [-0.09, 0.72] | 0.53 | 0.29 | 0.13 | 0.001 |
| PGLS (means, >= 3) | 60 | 0.46 [-0.05, 0.97] | 0.53 | 0.77 | 0.079 | 0.037 |

- The slope equals what the shorter-arc cap and 6 kb floor produce with no biology (S0).
  The positive slope vs 0 (p = 0.033, as in the team run: 0.49) is fully accounted for by geometry.
- Same in every sensitivity check (high only, all tiers, largest call, chromosome length:
  p vs S0 0.39-1.0) and leaving out any phylum (slopes 0.48-0.65, p vs S0 >= 0.70).
- The CI is wide (about +-0.5): a moderate biological effect cannot be excluded either.

## Position on the chromosome

| Class | Families | Observed share | Expected share (null) | Slope (PMM) | S0 | p vs S0 |
|---|---|---|---|---|---|---|
| Spanning origin/terminus | 419 | 62% | 35% [32, 39] | 0.49 [0.25, 0.73] | 0.77 | 0.025 |
| Local (one arm) | 257 | 38% | 65% | -0.58 [-1.07, -0.09] | 0.62 | 1e-5 |

- Clear excess of inversions spanning the origin or terminus (62% vs 35% expected).
- Spanning inversions grow LESS with genome size than geometry predicts.
- Local inversions get SMALLER in larger genomes (negative slope, far below S0).

### Checks on the local slope (`hpc/check_local.sh`)

| Check | Species | Slope [95% CI] | p vs 0 | p vs S0 (0.62) |
|---|---|---|---|---|
| As reported (high + medium) | 54 | -0.58 [-1.07, -0.09] | 0.022 | 1e-5 |
| High tier only | 46 | -0.35 [-0.89, 0.20] | 0.20 | 0.0009 |
| Species with >= 3 local families | 31 | -0.62 [-1.16, -0.08] | 0.026 | 6e-5 |
| Without Bacillota | 41 | -0.60 | 0.021 | 2e-5 |
| Without Campylobacterota | 50 | -0.56 | 0.05 | 1e-4 |
| Without Pseudomonadota | 19 | -0.32 [-1.24, 0.60] | 0.60 | 0.048 |

- No single species drives it (leaving any one out moves the slope by at most 0.08).
- Median local size drops from 139 kb (genomes <= 2 Mb) to 42-61 kb (3-6 Mb), but the 90th
  percentile stays at 700-900 kb: large genomes still have large local inversions, they also have
  many more small ones (4-5 Mb bin: 100 families, median 52 kb).
- Robust: local inversions are much smaller than geometry predicts (p vs S0 < 0.05 in every check).
- Less certain: whether the slope is truly negative or flat (high tier only and without
  Pseudomonadota include 0). Pseudomonadota are 35 of the 54 species.
- Hypothesis consistent with both position results: selection on replichore structure. A local
  inversion moves genes to the other strand relative to replication, costlier the more genes it
  moves, so large local inversions are disfavoured; spanning inversions that are symmetric around
  the origin or terminus keep both arms balanced and are favoured.

## Caveats
- Tiers provisional until the raw-alignment check (Phase 6).
- Dominated by Pseudomonadota (45 species) and Bacillota (24); 5 other phyla have <= 4 species.
- S0 treats the pooled observed sizes as exchangeable between species.
