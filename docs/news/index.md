# Changelog

## r4pde 0.2.3

### New features

- [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md)
  gains `env_smooth` and `ge_term` (with `k_ge`). `env_smooth = TRUE`
  replaces the environment random intercept with a random environment
  factor-smooth, `s(time, env, bs = "fs")`, which adjusts for
  environments that differ in epidemic timing or shape and for
  treatments that are unbalanced across environments.
  `ge_term = "intercept"` or `"smooth"` adds a genotype-by-environment
  random intercept or factor-smooth per treatment x environment curve.
  Defaults keep the previous model. `k_env` now sets the basis dimension
  of the environment factor-smooth.
- Treatment mean curves average environment factor-smooth and G x E
  deviations over environments instead of dropping them, because
  factor-smooths are not centred. The same construction is shared by
  [`functional_rate()`](https://emdelponte.github.io/r4pde/reference/functional_rate.md)
  and
  [`diagnose_curves()`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md),
  so fitted rates match the curves.

### Bug fixes

- [`functional_instability()`](https://emdelponte.github.io/r4pde/reference/functional_instability.md)
  now predicts genotype-by-environment curves with the G x E term,
  averages out plot-level unit effects, and no longer fails when
  [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md)
  was fitted with `block` or `unit`. It warns when the model has no G x
  E term, because nFI then mainly reflects each genotype’s disease
  level.
- [`functional_distances()`](https://emdelponte.github.io/r4pde/reference/functional_distances.md)
  (curve level) and `functional_resistance(group_method = "bootstrap")`
  supply all model variables when predicting curve-level trajectories;
  the bootstrap in
  [`functional_resistance()`](https://emdelponte.github.io/r4pde/reference/functional_resistance.md)
  previously failed when an environment was modelled.
- [`functional_distances()`](https://emdelponte.github.io/r4pde/reference/functional_distances.md):
  `perm_unit` is now honoured. Labels are permuted between permutation
  units (for example genotypes) and mapped back to their curves, instead
  of being shuffled between curves; this matters when `test_factor` is a
  genotype-level property such as a resistance class. `perm_strata` is
  applied at the unit level when strata are constant within units;
  otherwise whole units are permuted without strata, with a warning.
- [`functional_distances()`](https://emdelponte.github.io/r4pde/reference/functional_distances.md)
  gains `seed` (default 1) and, like `simulate_dsp_data(seed = 123)`, no
  longer resets the user’s random stream: the previous `.Random.seed` is
  restored afterwards. The bootstrap seed (`boot_seed`) is handled the
  same way.
- `functional_resistance(instability = )` accepts the list returned by
  `functional_instability(return_curves = TRUE)` (it read a non-existent
  `nfi_table` element) and joins on the treatment column instead of
  guessing it.
- `functional_pca_clusters(choose_k = "silhouette" | "elbow")` is
  implemented; it previously used k = 2 silently. The criterion table is
  attached as attribute `"k_selection"`.
- [`reconstruct_curves()`](https://emdelponte.github.io/r4pde/reference/reconstruct_curves.md)
  multiplies back by the column scaling when
  `functional_pca(scale = TRUE)` was used, and returns numeric times.
- [`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md): the
  summary mean shape index now averages the per-focus
  `SIF = rows / cols` (it averaged `cols / rows`).
- [`functional_contrast()`](https://emdelponte.github.io/r4pde/reference/functional_contrast.md)
  averages replicated observations explicitly, with a message, instead
  of letting [`approx()`](https://rdrr.io/r/stats/approxfun.html)
  average ties silently, and reports when constant extrapolation is
  used.
- [`functional_summary()`](https://emdelponte.github.io/r4pde/reference/functional_summary.md):
  `persistence` is the total time above `threshold`, summed over
  sampling intervals with the profile above the threshold at both ends;
  gaps below the threshold are no longer counted. Results for profiles
  without gaps are unchanged.
- [`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md):
  `results_table` is a data frame with numeric columns; parameters are
  no longer rounded before back-transformation or before drawing the
  fitted curves. The documentation no longer claims the models are
  ranked.
- [`theme_r4pde()`](https://emdelponte.github.io/r4pde/reference/theme_r4pde.md)
  works without cowplot (in Suggests), falling back to a
  [`theme_classic()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)-based
  equivalent.
- Documentation: `FHBWheat$i` and `$n` are the diseased and assessed
  plant counts used by
  [`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md);
  `BlastWheat$index_mean` is the wheat blast index.
- Tests: the
  [`CompMuCens()`](https://emdelponte.github.io/r4pde/reference/CompMuCens.md)
  test checked a non-existent element.

## r4pde 0.2.2

### New functions

- Added
  [`functional_rate()`](https://emdelponte.github.io/r4pde/reference/functional_rate.md):
  Estimates instantaneous rates of plant disease progress (first
  derivative of trajectory over time) and associated uncertainty from
  epidemic trajectories fitted by
  [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md).
  Supports response and link scale derivatives via GAM linear predictor
  matrix differences, handles boundary conditions without extrapolation,
  accounts for zero-epidemic curves and negative rates, and computes key
  rate phenotypes (`r_max`, `t_r_max`, growth duration, cumulative
  positive growth, and growth windows). Includes S3 methods:
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  [`augment()`](https://emdelponte.github.io/r4pde/reference/augment.md).

### Bug fixes and improvements

- [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md):
  Fixed smooth term label inspection during prediction so that `.unit`
  experimental unit and block random effects are correctly identified
  and supplied to `predict.gam()`.

## r4pde 0.2.1 (2026-06-14)

### New functions

- Added a suite of tools for Disease Suppression Profiles (DSP) to
  evaluate fungicide efficacy over time:
  [`functional_contrast()`](https://emdelponte.github.io/r4pde/reference/functional_contrast.md),
  [`functional_summary()`](https://emdelponte.github.io/r4pde/reference/functional_summary.md),
  [`plot_dsp()`](https://emdelponte.github.io/r4pde/reference/plot_dsp.md),
  [`rank_dsp()`](https://emdelponte.github.io/r4pde/reference/rank_dsp.md),
  [`plot_dsp_rank_heatmap()`](https://emdelponte.github.io/r4pde/reference/plot_dsp_rank_heatmap.md),
  and
  [`simulate_dsp_data()`](https://emdelponte.github.io/r4pde/reference/simulate_dsp_data.md).

### Bug fixes and improvements

- [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md):
  Updated the GAM formula to fit independent smooths per treatment by
  default (`trt + s(time, by=trt)`), avoiding artifactual drops
  (“caimento”) at the tails of flat curves. Added the `global_smooth`
  parameter to allow users to opt-in to the previous Global-Specific
  smoothing behavior.
- [`simulate_dsp_data()`](https://emdelponte.github.io/r4pde/reference/simulate_dsp_data.md):
  Adjusted the amplitudes of the simulated suppression profiles
  (“Early”, “Late”, “Persistent”) to achieve exactly a 70% reduction in
  AUDPC compared to the unsprayed control, representing a realistic
  fungicide efficacy scenario.
- Updated `_pkgdown.yml` to include the newly exported functional PCA
  and DSP functions in the “Functional Analysis” reference index,
  resolving site build errors.

## r4pde 0.2.0 (2026-04-26)

### Major changes

- Introduced a modular functional-analysis workflow based on
  [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md),
  [`functional_distances()`](https://emdelponte.github.io/r4pde/reference/functional_distances.md),
  [`functional_pca()`](https://emdelponte.github.io/r4pde/reference/functional_pca.md),
  [`functional_resistance()`](https://emdelponte.github.io/r4pde/reference/functional_resistance.md),
  and related methods.

### New functions

- [`functional_pca()`](https://emdelponte.github.io/r4pde/reference/functional_pca.md):
  Performs functional principal component analysis on fitted disease
  progress curves to decompose variation among epidemic trajectories
  into orthogonal temporal components. Includes plotting and extractor
  functions.
- [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md):
  Fits genotype-specific epidemic trajectories using GAM, with support
  for genotype-level covariates.
- [`functional_distances()`](https://emdelponte.github.io/r4pde/reference/functional_distances.md):
  Computes pairwise functional distances among fitted curves and
  performs hierarchical clustering and permutation testing.
- [`functional_resistance()`](https://emdelponte.github.io/r4pde/reference/functional_resistance.md):
  Calculates Functional Resistance Index (FRI) and Stability-Adjusted
  Functional Resistance Index (SAFRI) with support for stratified
  rankings and bootstrap-supported classification.
- [`suggest_k()`](https://emdelponte.github.io/r4pde/reference/suggest_k.md):
  Helper function to recommend GAM smoothing parameters (`k_smooth`,
  `k_trt`, `k_env`, `gamma`) based on data structure; supports tidy-eval
  column names and can infer replication from a data frame.

### Enhancements

- Added support for genotype-level auxiliary covariates (e.g.,
  `heading_group`) to adjust functional curves and resistance rankings,
  allowing for better distinction between genetic resistance and
  phenological escape.
- Improved bootstrap methodology for resistance classification with
  support for stratified group comparisons.
- Updated
  [`diagnose_curves()`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md),
  [`plot_curves()`](https://emdelponte.github.io/r4pde/reference/plot_curves.md),
  and
  [`plot_dendrogram()`](https://emdelponte.github.io/r4pde/reference/plot_dendrogram.md)
  to support the new functional analysis object classes.

## r4pde 0.1.2 (2026-04-03)

### New functions

- Added
  [`functional_instability()`](https://emdelponte.github.io/r4pde/reference/functional_instability.md)
  for computing normalized functional instability (NFI) of
  genotype-by-environment epidemic curves, with optional decomposition
  into spatial and temporal components.
- Added
  [`plot_functional_instability()`](https://emdelponte.github.io/r4pde/reference/plot_functional_instability.md)
  for visualizing overall, spatial, and temporal NFI results as bar
  charts.

## r4pde 0.1.1 (2026-02-05)

### Enhancements

- Added
  [`get_era5()`](https://emdelponte.github.io/r4pde/reference/get_era5.md)
  for retrieving ERA5 reanalysis weather data via Open-Meteo API.
- Added
  [`get_brdwgd()`](https://emdelponte.github.io/r4pde/reference/get_brdwgd.md)
  for extracting daily weather data from the Brazilian Daily Weather
  Gridded Data (BR-DWGD) NetCDF files.
- Added plotting and diagnostic helpers for curve models.

## r4pde 0.1.0 (2025-06-21)

### First release

- Initial release of the `r4pde` package.
