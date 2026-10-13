# r4pde

[![CRAN](https://www.r-pkg.org/badges/version/r4pde)](https://CRAN.R-project.org/package=r4pde)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/r4pde)](https://CRAN.R-project.org/package=r4pde)

**r4pde** provides tools for quantitative plant disease epidemiology. It
started as the companion package to the book [R for Plant Disease
Epidemiology](https://r4pde.netlify.app/) and now also implements
methods from recent research on the functional analysis of disease
progress curves. It is aimed at researchers, students and practitioners
who analyse field and greenhouse disease data.

## What the package does

| Area | What you can do | Main functions |
|----|----|----|
| **Disease progress curves (functional analysis)** | Fit smooth, design-adjusted disease progress curves with GAMs; estimate instantaneous rates; compare and cluster whole trajectories; rank resistance and quantify stability across environments | [`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md), [`functional_rate()`](https://emdelponte.github.io/r4pde/reference/functional_rate.md), [`functional_distances()`](https://emdelponte.github.io/r4pde/reference/functional_distances.md), [`functional_pca()`](https://emdelponte.github.io/r4pde/reference/functional_pca.md), [`functional_resistance()`](https://emdelponte.github.io/r4pde/reference/functional_resistance.md), [`functional_instability()`](https://emdelponte.github.io/r4pde/reference/functional_instability.md), [`diagnose_curves()`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md) |
| **Disease management** | Build Disease Suppression Profiles (DSP) for fungicides or other control measures relative to an untreated check, summarise and rank them | [`functional_contrast()`](https://emdelponte.github.io/r4pde/reference/functional_contrast.md), [`functional_summary()`](https://emdelponte.github.io/r4pde/reference/functional_summary.md), [`functional_suppression_profiles()`](https://emdelponte.github.io/r4pde/reference/functional_suppression_profiles.md), [`rank_dsp()`](https://emdelponte.github.io/r4pde/reference/rank_dsp.md), [`plot_dsp()`](https://emdelponte.github.io/r4pde/reference/plot_dsp.md) |
| **Spatial epidemiology** | Test for aggregation (ordinary runs, join counts, binary power law), analyse foci (AFSD), count subareas, fit dispersal gradients | [`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md), [`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md), [`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md), [`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md), [`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md), [`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md) |
| **Disease quantification** | Compute disease severity indices from ordinal scales and compare treatments on censored ordinal data | [`DSI()`](https://emdelponte.github.io/r4pde/reference/DSI.md), [`DSI2()`](https://emdelponte.github.io/r4pde/reference/DSI2.md), [`CompMuCens()`](https://emdelponte.github.io/r4pde/reference/CompMuCens.md) |
| **Weather and environment** | Download daily weather (NASA POWER, ERA5, BR-DWGD) and build window-pane variables for weather–disease models | [`get_nasapower()`](https://emdelponte.github.io/r4pde/reference/get_nasapower.md), [`get_era5()`](https://emdelponte.github.io/r4pde/reference/get_era5.md), [`get_brdwgd()`](https://emdelponte.github.io/r4pde/reference/get_brdwgd.md), [`windowpane()`](https://emdelponte.github.io/r4pde/reference/windowpane.md), [`windowpane_tests()`](https://emdelponte.github.io/r4pde/reference/windowpane_tests.md) |
| **Datasets** | Real epidemiological datasets used in the book and in the examples | `BudBlightSoybean`, `RustSoybean`, `FHBWheat`, `WhiteMoldSoybean`, `SpatialAggregated`, … |

The [reference page](https://emdelponte.github.io/r4pde/reference/)
lists every function by family.

## Installation

Install the stable release from CRAN:

``` r

install.packages("r4pde")
```

The development version is on GitHub. To install it with its
dependencies (including Bioconductor packages), use
[`pak`](https://pak.r-lib.org/):

``` r

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install("emdelponte/r4pde")
```

Or with remotes:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("emdelponte/r4pde")
```

## Learning more

- **Book:** [R for Plant Disease
  Epidemiology](https://r4pde.netlify.app/) explains the epidemiological
  concepts and uses this package throughout.
- **Vignette:** [Functional analysis of disease progress
  curves](https://emdelponte.github.io/r4pde/articles/functional-analysis.html)
  walks through the full workflow, from fitting curves to resistance and
  stability across environments.
- **Website:**
  [emdelponte.github.io/r4pde](https://emdelponte.github.io/r4pde/) has
  the function reference and the changelog.
- **Methods:** many functions follow Madden, Hughes and van den Bosch
  (2007), *The Study of Plant Disease Epidemics*
  ([doi:10.1094/9780890545058](https://doi.org/10.1094/9780890545058)).

## Reporting issues

Please report bugs or suggest features on the [GitHub issue
tracker](https://github.com/emdelponte/r4pde/issues). A small
reproducible example helps a lot.

## License and citation

r4pde is released under the MIT license. To cite it, run:

``` r

citation("r4pde")
```
