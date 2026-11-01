# r4pde <img src="man/figures/logo.png" width="200px" align="right"/>

[![CRAN](https://www.r-pkg.org/badges/version/r4pde)](https://CRAN.R-project.org/package=r4pde) [![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/r4pde)](https://CRAN.R-project.org/package=r4pde)

**r4pde** provides tools for quantitative plant disease epidemiology. It started as the companion package to the book [R for Plant Disease Epidemiology](https://r4pde.netlify.app/) and now also implements methods from recent research on the functional analysis of disease progress curves and on the validation of standard area diagrams. It is aimed at researchers, students and practitioners who analyse field and greenhouse disease data.

## What the package does

| Area | What you can do | Main functions |
|---|---|---|
| **Disease progress curves (functional analysis)** | Fit smooth, design-adjusted disease progress curves with GAMs; estimate instantaneous rates; compare and cluster whole trajectories; rank resistance and quantify stability across environments | `functional_curves()`, `functional_rate()`, `functional_distances()`, `functional_pca()`, `functional_resistance()`, `functional_instability()`, `diagnose_curves()` |
| **Disease management** | Build Disease Suppression Profiles (DSP) for fungicides or other control measures relative to an untreated check, summarise and rank them | `functional_contrast()`, `functional_summary()`, `functional_suppression_profiles()`, `rank_dsp()`, `plot_dsp()` |
| **Spatial epidemiology** | Test for aggregation (ordinary runs, join counts, binary power law), analyse foci (AFSD), count subareas, fit dispersal gradients | `oruns_test()`, `join_count()`, `BPL()`, `AFSD()`, `count_subareas()`, `fit_gradients()` |
| **Disease quantification** | Compute disease severity indices from ordinal scales and compare treatments on censored ordinal data | `DSI()`, `DSI2()`, `CompMuCens()` |
| **Validation of standard area diagrams (SADs)** | Conventional analysis with Lin's concordance per rater, comparison of conditions and inter-rater reliability; one-stage hierarchical model of all estimates, with bias, variation among raters and expected error along the severity range, functional comparison of SADs, rater and image effects, and power simulation for planning studies | `sad_data()`, `sad_ccc()`, `sad_compare()`, `sad_reliability()`, `sad_plots()`, `sad_fit()`, `sad_curve()`, `sad_regions()`, `sad_raters()`, `sad_images()`, `sad_simulate()` |
| **Weather and environment** | Download daily weather (NASA POWER, ERA5, BR-DWGD) and build window-pane variables for weather–disease models | `get_nasapower()`, `get_era5()`, `get_brdwgd()`, `windowpane()`, `windowpane_tests()` |
| **Datasets** | Real epidemiological datasets used in the book and in the examples | `BudBlightSoybean`, `RustSoybean`, `FHBWheat`, `WhiteMoldSoybean`, `SpatialAggregated`, `SADRustSoybean`, … |

The [reference page](https://emdelponte.github.io/r4pde/reference/) lists every function by family.

## Installation

Install the stable release from CRAN:

``` r
install.packages("r4pde")
```

The development version is on GitHub. To install it with its dependencies (including Bioconductor packages), use [`pak`](https://pak.r-lib.org/):

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

- **Book:** [R for Plant Disease Epidemiology](https://r4pde.netlify.app/) explains the epidemiological concepts and uses this package throughout.
- **Vignettes:**
  - [Functional analysis of disease progress curves](https://emdelponte.github.io/r4pde/articles/functional-analysis.html) walks through the full workflow, from fitting curves to resistance and stability across environments.
  - [Validation of standard area diagrams](https://emdelponte.github.io/r4pde/articles/sad-validation.html) walks through the conventional and the one-stage analysis of a SAD validation study, step by step.
- **Website:** [emdelponte.github.io/r4pde](https://emdelponte.github.io/r4pde/) has the function reference and the changelog.
- **Methods:** many functions follow Madden, Hughes and van den Bosch (2007), *The Study of Plant Disease Epidemics* ([doi:10.1094/9780890545058](https://doi.org/10.1094/9780890545058)).

## Reporting issues

Please report bugs or suggest features on the [GitHub issue tracker](https://github.com/emdelponte/r4pde/issues). A small reproducible example helps a lot.

## License and citation

r4pde is released under the MIT license. To cite it, run:

``` r
citation("r4pde")
```
