# Changelog

## r4pde 0.2.0 (2026-04-26)

### Major changes

- Refactored the monolithic
  [`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md)
  workflow into a modular functional analysis API.
- [`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md)
  is now a soft-deprecated wrapper around the new modular functions.

### New functions

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
- Added
  [`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md)
  for functional comparison and clustering of epidemic curves.
- Added plotting and diagnostic helpers for curve models.

## r4pde 0.1.0 (2025-06-21)

CRAN release: 2025-07-02

### First release

- Initial release of the `r4pde` package.
