# Exploratory plots for SAD validation

Produces the exploratory figures of a conventional SAD validation study:
estimates against actual severity, Cb against r, boxplots of Lin's
statistics, errors along the severity range, error densities, the
standard deviation among raters along the severity range and the gain in
accuracy relative to unaided estimates.

## Usage

``` r
sad_plots(x, ccc = sad_ccc(x), baseline = NULL, palette = NULL)
```

## Arguments

- x:

  An object returned by
  [`sad_data()`](https://emdelponte.github.io/r4pde/reference/sad_data.md).

- ccc:

  Statistics per rater, by default `sad_ccc(x)`.

- baseline:

  Optional name of the method used as baseline (for example `"UN"`);
  when given, the plot of gain in \\\rho_c\\ against baseline accuracy
  is added.

- palette:

  Optional named vector of colours, one per condition.

## Value

A named list of ggplot objects: `scatter`, `cb_vs_r`, `box_ccc`,
`error`, `error_density` (requires ggridges), `sd_vs_actual` and, with
`baseline`, `gain_vs_baseline`.

## Examples

``` r
x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
p <- sad_plots(x, baseline = "UN")
p$scatter
```
