# Package index

## Disease Quantification

Functions for calculating disease indices and survival-based metrics.

- [`DSI()`](https://emdelponte.github.io/r4pde/reference/DSI.md) :
  Calculate the Disease Severity Index (DSI) (class for each unit)
- [`DSI2()`](https://emdelponte.github.io/r4pde/reference/DSI2.md) :
  Calculate the Disease severity Index (DSI) (frequency of each class)
- [`CompMuCens()`](https://emdelponte.github.io/r4pde/reference/CompMuCens.md)
  : Survival analysis for quantitative ordinal scale data.

## Spatial analysis

Functions to detect and quantify spatial structure in disease data.

- [`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md) :
  Analysis of foci structure and dynamics (AFSD)
- [`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md) :
  Binary Power Law Analysis for Spatial Disease Patterns
- [`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md)
  : Count the Number of Ones in Subareas of a Matrix
- [`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md)
  : Random Subgrid Sampling of a Binary Matrix
- [`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md)
  : Test for Spatial Join Count Statistics
- [`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md)
  : Runs Test
- [`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md)
  : Runs Test for Each Row and Column of a Binary Matrix
- [`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md)
  : Boustrophedon Run Test for Binary Matrix
- [`windowpane()`](https://emdelponte.github.io/r4pde/reference/windowpane.md)
  : Window Pane for Epidemiological Analysis
- [`windowpane_tests()`](https://emdelponte.github.io/r4pde/reference/windowpane_tests.md)
  : Windowpane Tests for Correlation Analysis
- [`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)
  : Plot ASFD

## Disease modeling

Functions for fitting gradient models and retrieving climate data.

- [`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md)
  : Fit Gradient Models to Data
- [`get_nasapower()`](https://emdelponte.github.io/r4pde/reference/get_nasapower.md)
  : Fetch NASA POWER Data for Multiple Locations with a Progress Bar

## Datasets

Datasets used in the examples and analyses.

- [`BlastWheat`](https://emdelponte.github.io/r4pde/reference/BlastWheat.md)
  : BlastWheat dataset
- [`BudBlightSoybean`](https://emdelponte.github.io/r4pde/reference/BudBlightSoybean.md)
  : BudBlightSoybean dataset
- [`DidymellaWatermelon`](https://emdelponte.github.io/r4pde/reference/DidymellaWatermelon.md)
  : DidymellaWatermelon dataset
- [`FHBWheat`](https://emdelponte.github.io/r4pde/reference/FHBWheat.md)
  : FHBWheat dataset
- [`FusariumBanana`](https://emdelponte.github.io/r4pde/reference/FusariumBanana.md)
  : FusariumBanana dataset
- [`RustSoybean`](https://emdelponte.github.io/r4pde/reference/RustSoybean.md)
  : RustSoybean dataset
- [`WhiteMoldSoybean`](https://emdelponte.github.io/r4pde/reference/WhiteMoldSoybean.md)
  : WhiteMoldSoybean dataset
- [`SpatialAggregated`](https://emdelponte.github.io/r4pde/reference/SpatialAggregated.md)
  : SpatialAggregated dataset
- [`SpatialRandom`](https://emdelponte.github.io/r4pde/reference/SpatialRandom.md)
  : SpatialRandom dataset

## Miscellaneous

Utility functions for plotting and customization.

- [`theme_r4pde()`](https://emdelponte.github.io/r4pde/reference/theme_r4pde.md)
  : Custom ggplot2 theme based on cowplot::theme_half_open

## Curve comparison

- [`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md)
  : Compare epidemic curves using functional GAM smoothing and
  clustering
- [`plot_curves()`](https://emdelponte.github.io/r4pde/reference/plot_curves.md)
  : Plot environment-adjusted epidemic curves by cluster
- [`plot_dendrogram()`](https://emdelponte.github.io/r4pde/reference/plot_dendrogram.md)
  : Plot functional dendrogram of epidemic curves
- [`plot(`*`<r4pde_compare_curves>`*`)`](https://emdelponte.github.io/r4pde/reference/plot.r4pde_compare_curves.md)
  : Plot method for compare_curves objects

## Diagnostics

- [`diagnose_curves()`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md)
  : Diagnostic tools for functional epidemic curve models
- [`plot_diagnostics()`](https://emdelponte.github.io/r4pde/reference/plot_diagnostics.md)
  : Plot diagnostics for functional epidemic curve models
- [`print(`*`<r4pde_curve_diagnostics>`*`)`](https://emdelponte.github.io/r4pde/reference/print.r4pde_curve_diagnostics.md)
  : Print method for curve diagnostics
