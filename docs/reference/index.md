# Package index

## Disease Quantification

Functions for calculating disease indices and survival-based metrics.

- [`CompMuCens()`](https://emdelponte.github.io/r4pde/reference/CompMuCens.md)
  : Survival analysis for quantitative ordinal scale data
- [`DSI()`](https://emdelponte.github.io/r4pde/reference/DSI.md) :
  Calculate the Disease Severity Index (DSI) (class for each unit)
- [`DSI2()`](https://emdelponte.github.io/r4pde/reference/DSI2.md) :
  Calculate the Disease severity Index (DSI) (frequency of each class)

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
- [`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md)
  : Boustrophedon Run Test for Binary Matrix
- [`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md)
  : Runs Test for Each Row and Column of a Binary Matrix
- [`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)
  : Plot ASFD
- [`windowpane()`](https://emdelponte.github.io/r4pde/reference/windowpane.md)
  : Window Pane for Epidemiological Analysis
- [`windowpane_tests()`](https://emdelponte.github.io/r4pde/reference/windowpane_tests.md)
  : Windowpane Tests for Correlation Analysis

## Disease modeling

Functions for fitting gradient models and retrieving climate data.

- [`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md)
  : Fit Gradient Models to Data
- [`get_brdwgd()`](https://emdelponte.github.io/r4pde/reference/get_brdwgd.md)
  : Fetch BR-DWGD Data for Multiple Locations
- [`get_era5()`](https://emdelponte.github.io/r4pde/reference/get_era5.md)
  : Fetch ERA5 Data from Open-Meteo for Multiple Locations
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
- [`SpatialAggregated`](https://emdelponte.github.io/r4pde/reference/SpatialAggregated.md)
  : SpatialAggregated dataset
- [`SpatialRandom`](https://emdelponte.github.io/r4pde/reference/SpatialRandom.md)
  : SpatialRandom dataset
- [`WhiteMoldSoybean`](https://emdelponte.github.io/r4pde/reference/WhiteMoldSoybean.md)
  : WhiteMoldSoybean dataset

## Miscellaneous

Utility functions for plotting and customization.

- [`theme_r4pde()`](https://emdelponte.github.io/r4pde/reference/theme_r4pde.md)
  : Custom ggplot2 theme based on cowplot::theme_half_open

## Functional Analysis

Functions for functional modeling and comparison of epidemic curves.

- [`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md)
  : Compare epidemic curves using GAM smoothing, functional distances,
  clustering, and optional permutation/bootstrapping
- [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md)
  : Fit genotype-specific epidemic trajectories using GAM
- [`functional_distances()`](https://emdelponte.github.io/r4pde/reference/functional_distances.md)
  : Compute pairwise functional distances and cluster curves
- [`functional_instability()`](https://emdelponte.github.io/r4pde/reference/functional_instability.md)
  : Normalized functional instability from compare_curves output
- [`functional_resistance()`](https://emdelponte.github.io/r4pde/reference/functional_resistance.md)
  : Compute functional resistance and stability-adjusted functional
  resistance
- [`plot(`*`<functional_curves>`*`)`](https://emdelponte.github.io/r4pde/reference/plot.functional_curves.md)
  : Plot functional_curves
- [`plot(`*`<r4pde_compare_curves>`*`)`](https://emdelponte.github.io/r4pde/reference/plot.r4pde_compare_curves.md)
  : Plot method for compare_curves objects
- [`plot_curves()`](https://emdelponte.github.io/r4pde/reference/plot_curves.md)
  : Plot environment-adjusted epidemic curves by cluster
- [`plot_dendrogram()`](https://emdelponte.github.io/r4pde/reference/plot_dendrogram.md)
  : Plot functional dendrogram of epidemic curves
- [`plot_functional_instability()`](https://emdelponte.github.io/r4pde/reference/plot_functional_instability.md)
  : Plot method for functional instability results
- [`print(`*`<functional_curves>`*`)`](https://emdelponte.github.io/r4pde/reference/print.functional_curves.md)
  : Print functional_curves
- [`print(`*`<functional_distances>`*`)`](https://emdelponte.github.io/r4pde/reference/print.functional_distances.md)
  : Print functional_distances
- [`print(`*`<functional_resistance>`*`)`](https://emdelponte.github.io/r4pde/reference/print.functional_resistance.md)
  : Print functional_resistance

## Diagnostics

- [`diagnose_curves()`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md)
  : Diagnostic tools for functional epidemic curve models
- [`plot_diagnostics()`](https://emdelponte.github.io/r4pde/reference/plot_diagnostics.md)
  : Plot diagnostics for functional epidemic curve models
- [`print(`*`<r4pde_curve_diagnostics>`*`)`](https://emdelponte.github.io/r4pde/reference/print.r4pde_curve_diagnostics.md)
  : Print method for curve diagnostics
