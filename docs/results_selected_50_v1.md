# Results: selected_50_v1 (final tiers, 21 Oct 2026)

79 GTDB species (min 50 genomes each), 7 phyla, genomes 1-7 Mb.
Slope = log10(family size, kb) on log10(genome size, Mb); each inversion family counted once per species.
Figures: `results/compare/selected_50_v1/` (made by `hpc/finish_50.sh`).

## Evidence and tiers

| Locus evidence | Loci | Backed by raw alignment |
|---|---|---|
| Leen + SyRI, matching breakpoints | 1,171 | 91% |
| Leen + SyRI, overlapping only | 491 | 81% |
| Leen only | 1,183 | 30% |
| SyRI only | 99 | 85% (same alignment as the check) |

Families: 379 high, 449 medium, 414 low; 828 high + medium in 77 species.
The raw-alignment check separates the classes the same way as the team's Part 6 (95 / 50 / 8%):
loci found by both callers are almost always backed; Leen-only loci mostly are not.

## Main question: does inversion size scale with genome size beyond geometry?

| Model | Species | Rerun slope [95% CI] | Team | S0 (geometry only) | p vs S0 | p vs 0 | p vs 1 |
|---|---|---|---|---|---|---|---|
| Mixed model with tree (main) | 77 | 0.46 [0.04, 0.89] | 0.49 [0.10, 0.88] | 0.52 [0.28, 0.78] | 0.78 | 0.032 | 0.014 |
| Mixed model, no tree | 77 | 0.32 [-0.05, 0.69] | | 0.52 | 0.27 | 0.09 | 0.0004 |
| PGLS, species means (>= 3) | 62 | 0.33 [-0.14, 0.80] | 0.41 [-0.07, 0.89] | 0.52 | 0.41 | 0.17 | 0.006 |
| Chromosome length | 77 | 0.45 [0.01, 0.90] | 0.38 | 0.52 | 0.76 | 0.046 | 0.017 |

- The rerun reproduces the team's slope (0.46 vs 0.49) with GTDB species, loci and new tiers.
- That slope is what the shorter-arc cap and the 6 kb floor produce with no biology (S0 = 0.52,
  p = 0.78). The "significant vs 0" result is explained by geometry.
- Same in every sensitivity check (p vs S0 0.26-0.92) and leaving out any phylum (0.41-0.74, p vs S0 >= 0.41).
- The slope is below proportional scaling (p = 0.014), but so is S0; the interval (about +-0.4)
  cannot exclude a moderate biological effect either.

## Position on the chromosome

| | Observed | Geometry only | p |
|---|---|---|---|
| Share of families spanning origin/terminus | 59% | 34% [31, 37] | far outside the null |
| Spanning slope | 0.51 [0.29, 0.72] (team 0.59) | 0.77 | 0.019 (PGLS 0.007) |
| Local slope | -0.42 [-0.87, 0.04] (team 0.09) | 0.60 | 3e-5 |

- Strong excess of inversions spanning the origin or terminus.
- Spanning inversions grow LESS with genome size than geometry predicts. The team's 0.59 is
  reproduced (0.51), but it is below, not above, what geometry alone gives.
- Local inversions do not grow with genome size (slope <= 0, far below S0). The sign is not
  robust: without Pseudomonadota (39 of 60 species) the slope is +0.12 [-0.73, 0.96].

## Caveats
- 69 of 77 species are Pseudomonadota or Bacillota; the min-20 run (11 phyla) is the next test.
- S0 treats pooled observed sizes as exchangeable between species.
- Leen-only loci without raw support (824) are low tier and excluded from the main model.
