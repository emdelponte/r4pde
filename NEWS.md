# r4pde 0.2.3

## New features
- New `sad_*` family for the validation of standard area diagrams (SADs), for studies with a single SAD (without and with the SAD) and for comparisons of two or more SADs or conditions:
  - `sad_data()` organises estimates in long or wide format and checks them.
  - Conventional two-stage workflow: `sad_ccc()` (Lin's concordance and its components per rater), `sad_compare()` (means by condition with Tukey letters, mixed model with rater random when raters are paired), `sad_reliability()` (ICC(A,1) with confidence interval and OCC) and `sad_plots()` (exploratory figures).
  - One-stage hierarchical model: `sad_fit()` (glmmTMB mixed model with rater intercepts and slopes per method, image difficulty and severity-dependent dispersion; gaussian or ordered beta), `sad_curve()` (bias, variance components and expected absolute error along the severity range, and differences between conditions with simultaneous bands; `type = "smooth"` fits a functional GAMM of the error curves), `sad_regions()`, `sad_raters()`, `plot_sad_raters()`, `sad_images()` and `sad_simulate()` (power for planning a study).
- New dataset `SADRustSoybean` (Franceschi et al. 2020) and vignette "Validation of standard area diagrams".
- glmmTMB, emmeans, lmerTest, multcompView, MASS and ggridges are suggested packages, needed only by the `sad_*` functions.
- `functional_curves()` gains `env_smooth` and `ge_term` (with `k_ge`). `env_smooth = TRUE` replaces the environment random intercept with a random environment factor-smooth, `s(time, env, bs = "fs")`, which adjusts for environments that differ in epidemic timing or shape and for treatments that are unbalanced across environments. `ge_term = "intercept"` or `"smooth"` adds a genotype-by-environment random intercept or factor-smooth per treatment x environment curve. Defaults keep the previous model. `k_env` now sets the basis dimension of the environment factor-smooth. Both factor-smooths use a first-order penalty (`m = 1`), so an environment or G x E deviation flattens out where a curve has no assessments instead of bending back up at the start of the time domain.
- Treatment mean curves average environment factor-smooth and G x E deviations over environments instead of dropping them, because factor-smooths are not centred. The same construction is shared by `functional_rate()` and `diagnose_curves()`, so fitted rates match the curves.
- `functional_suppression_profiles()` gains `env_smooth`, `ge_term` and `k_ge`, passed to `functional_curves()` in the joint analysis across environments. `ge_term = "intercept"` or `"smooth"` adds a treatment-by-environment random intercept or factor-smooth per suppression curve, as already available for `functional_instability()`. The treatment mean suppression curves average those deviations over the environments where each treatment was observed. Defaults keep the previous model. With `by_environment = TRUE` the terms are ignored, with a warning.

## Bug fixes
- `functional_instability()` now predicts genotype-by-environment curves with the G x E term, averages out plot-level unit effects, and no longer fails when `functional_curves()` was fitted with `block` or `unit`. It warns when the model has no G x E term, because nFI then mainly reflects each genotype's disease level.
- `functional_distances()` (curve level) and `functional_resistance(group_method = "bootstrap")` supply all model variables when predicting curve-level trajectories; the bootstrap in `functional_resistance()` previously failed when an environment was modelled.
- `functional_distances()`: `perm_unit` is now honoured. Labels are permuted between permutation units (for example genotypes) and mapped back to their curves, instead of being shuffled between curves; this matters when `test_factor` is a genotype-level property such as a resistance class. `perm_strata` is applied at the unit level when strata are constant within units; otherwise whole units are permuted without strata, with a warning.
- `functional_distances()` gains `seed` (default 1) and, like `simulate_dsp_data(seed = 123)`, no longer resets the user's random stream: the previous `.Random.seed` is restored afterwards. The bootstrap seed (`boot_seed`) is handled the same way.
- `functional_resistance(instability = )` accepts the list returned by `functional_instability(return_curves = TRUE)` (it read a non-existent `nfi_table` element) and joins on the treatment column instead of guessing it.
- `functional_pca_clusters(choose_k = "silhouette" | "elbow")` is implemented; it previously used k = 2 silently. The criterion table is attached as attribute `"k_selection"`.
- `reconstruct_curves()` multiplies back by the column scaling when `functional_pca(scale = TRUE)` was used, and returns numeric times.
- `AFSD()`: the summary mean shape index now averages the per-focus `SIF = rows / cols` (it averaged `cols / rows`).
- `functional_contrast()` averages replicated observations explicitly, with a message, instead of letting `approx()` average ties silently, and reports when constant extrapolation is used.
- `functional_summary()`: `persistence` is the total time above `threshold`, summed over sampling intervals with the profile above the threshold at both ends; gaps below the threshold are no longer counted. Results for profiles without gaps are unchanged.
- `fit_gradients()`: `results_table` is a data frame with numeric columns; parameters are no longer rounded before back-transformation or before drawing the fitted curves. The documentation no longer claims the models are ranked.
- `theme_r4pde()` works without cowplot (in Suggests), falling back to a `theme_classic()`-based equivalent.
- Documentation: `FHBWheat$i` and `$n` are the diseased and assessed plant counts used by `BPL()`; `BlastWheat$index_mean` is the wheat blast index.
- Tests: the `CompMuCens()` test checked a non-existent element.
- `functional_suppression_profiles()` works when the treatment column is not named `"treatment"`; the column name was hard-coded in the cluster table and in a join (also in `print()`).
- `functional_suppression_profiles()` with `environment`: the metrics are computed from the fitted suppression curves returned to the scale of `response` (they were on the 0-1 proportion scale used internally by `functional_curves()`). `threshold` and the metrics now have the same units with and without `environment`; previously the default `threshold = 5` gave `persistence = 0` for every treatment with percent data. `protected_area` and `max_suppression` from joint fits on percent data are now 100 times larger than before.

# r4pde 0.2.2

## New functions
- Added `functional_rate()`: Estimates instantaneous rates of plant disease progress (first derivative of trajectory over time) and associated uncertainty from epidemic trajectories fitted by `functional_curves()`. Supports response and link scale derivatives via GAM linear predictor matrix differences, handles boundary conditions without extrapolation, accounts for zero-epidemic curves and negative rates, and computes key rate phenotypes (`r_max`, `t_r_max`, growth duration, cumulative positive growth, and growth windows). Includes S3 methods: `print()`, `summary()`, `plot()`, and `augment()`.

## Bug fixes and improvements
- `functional_curves()`: Fixed smooth term label inspection during prediction so that `.unit` experimental unit and block random effects are correctly identified and supplied to `predict.gam()`.

# r4pde 0.2.1 (2026-06-14)

## New functions
- Added a suite of tools for Disease Suppression Profiles (DSP) to evaluate fungicide efficacy over time: `functional_contrast()`, `functional_summary()`, `plot_dsp()`, `rank_dsp()`, `plot_dsp_rank_heatmap()`, and `simulate_dsp_data()`.

## Bug fixes and improvements
- `functional_curves()`: Updated the GAM formula to fit independent smooths per treatment by default (`trt + s(time, by=trt)`), avoiding artifactual drops ("caimento") at the tails of flat curves. Added the `global_smooth` parameter to allow users to opt-in to the previous Global-Specific smoothing behavior.
- `simulate_dsp_data()`: Adjusted the amplitudes of the simulated suppression profiles ("Early", "Late", "Persistent") to achieve exactly a 70% reduction in AUDPC compared to the unsprayed control, representing a realistic fungicide efficacy scenario.
- Updated `_pkgdown.yml` to include the newly exported functional PCA and DSP functions in the "Functional Analysis" reference index, resolving site build errors.

# r4pde 0.2.0 (2026-04-26)

## Major changes
- Introduced a modular functional-analysis workflow based on `functional_curves()`, `functional_distances()`, `functional_pca()`, `functional_resistance()`, and related methods.

## New functions
- `functional_pca()`: Performs functional principal component analysis on fitted disease progress curves to decompose variation among epidemic trajectories into orthogonal temporal components. Includes plotting and extractor functions.
- `functional_curves()`: Fits genotype-specific epidemic trajectories using GAM, with support for genotype-level covariates.
- `functional_distances()`: Computes pairwise functional distances among fitted curves and performs hierarchical clustering and permutation testing.
- `functional_resistance()`: Calculates Functional Resistance Index (FRI) and Stability-Adjusted Functional Resistance Index (SAFRI) with support for stratified rankings and bootstrap-supported classification.
- `suggest_k()`: Helper function to recommend GAM smoothing parameters (`k_smooth`, `k_trt`, `k_env`, `gamma`) based on data structure; supports tidy-eval column names and can infer replication from a data frame.

## Enhancements
- Added support for genotype-level auxiliary covariates (e.g., `heading_group`) to adjust functional curves and resistance rankings, allowing for better distinction between genetic resistance and phenological escape.
- Improved bootstrap methodology for resistance classification with support for stratified group comparisons.
- Updated `diagnose_curves()`, `plot_curves()`, and `plot_dendrogram()` to support the new functional analysis object classes.

# r4pde 0.1.2 (2026-04-03)

## New functions
- Added `functional_instability()` for computing normalized functional instability (NFI) of genotype-by-environment epidemic curves, with optional decomposition into spatial and temporal components.
- Added `plot_functional_instability()` for visualizing overall, spatial, and temporal NFI results as bar charts.

# r4pde 0.1.1 (2026-02-05)

## Enhancements
- Added `get_era5()` for retrieving ERA5 reanalysis weather data via Open-Meteo API.
- Added `get_brdwgd()` for extracting daily weather data from the Brazilian Daily Weather Gridded Data (BR-DWGD) NetCDF files.
- Added plotting and diagnostic helpers for curve models.

# r4pde 0.1.0 (2025-06-21)

## First release
- Initial release of the `r4pde` package.
