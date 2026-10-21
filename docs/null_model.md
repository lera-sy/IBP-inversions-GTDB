# Null model notes

## 3A. Geometry null (21 Oct 2026)

`src/ibp/nullmodel/geometry.py`, rule `null_geometry`. 80 species, chromosomes
log-uniform 1-10 Mb, 30 events each, 6 kb floor, shorter-arc rule, 100 replicates.
Slope of log10(size, kb) on log10(chromosome length), mean [2.5-97.5%].

| Null | Sizes | All | Local | Spanning | Spanning fraction |
|---|---|---|---|---|---|
| Random breakpoints | scale with chromosome | 0.98 [0.93, 1.03] | 0.96 | 1.00 | 0.50 |
| Fixed size | log-uniform 6 kb-5 Mb | 0.46 [0.36, 0.57] | 0.44 | 0.97 | 0.22 |
| Fixed size | lognormal, median 200 kb | 0.19 [0.14, 0.26] | 0.30 | 0.61 | 0.22 |
| Fixed size | lognormal, median 30 kb | 0.02 [-0.05, 0.09] | 0.09 | 0.25 | 0.06 |

What this means
- With no biological link at all, the cap and floor alone give slopes from about 0 to 0.46,
  depending on how heavy the tail of the size distribution is.
- Spanning inversions get a steep slope by construction (0.25-0.97): a large inversion is
  more likely to cross the origin or terminus, and the cap is tighter on small chromosomes.
- So the team's results (all 0.49, spanning 0.59, local 0.09) are within the range geometry
  alone can produce. They are not evidence of biology until tested against S0 from the
  real size distribution.

Next: rerun with `--species-table` (real chromosome lengths) and `--empirical-sizes`
(pooled real sizes, shuffled across species) once Phase 5 calls exist. That empirical
fixed-size slope is S0.
