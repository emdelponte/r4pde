# r4pde 0.2.1 (2026-06-14)

## New functions
- Added a suite of tools for Disease Suppression Profiles (DSP) to evaluate fungicide efficacy over time: `functional_contrast()`, `functional_summary()`, `plot_dsp()`, `rank_dsp()`, `plot_dsp_rank_heatmap()`, and `simulate_dsp_data()`.

## Bug fixes and improvements
- `functional_curves()`: Updated the GAM formula to fit independent smooths per treatment by default (`trt + s(time, by=trt)`), avoiding artifactual drops ("caimento") at the tails of flat curves. Added the `global_smooth` parameter to allow users to opt-in to the previous Global-Specific smoothing behavior.
- `simulate_dsp_data()`: Adjusted the amplitudes of the simulated suppression profiles ("Early", "Late", "Persistent") to achieve exactly a 70% reduction in AUDPC compared to the unsprayed control, representing a realistic fungicide efficacy scenario.
- Updated `_pkgdown.yml` to include the newly exported functional PCA and DSP functions in the "Functional Analysis" reference index, resolving site build errors.

# r4pde 0.2.0 (2026-04-26)

## Major changes
- Refactored the monolithic `compare_curves()` workflow into a modular functional analysis API.
- `compare_curves()` is now a soft-deprecated wrapper around the new modular functions.

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
- Added `compare_curves()` for functional comparison and clustering of epidemic curves.
- Added plotting and diagnostic helpers for curve models.

# r4pde 0.1.0 (2025-06-21)

## First release
- Initial release of the `r4pde` package.
