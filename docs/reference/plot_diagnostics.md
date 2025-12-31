# Plot diagnostics for functional epidemic curve models

Combines multiple diagnostic plots from
[`diagnose_curves`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md)
into a multi-panel layout.

## Usage

``` r
plot_diagnostics(
  diag,
  layout = c("2x2", "vertical", "horizontal"),
  show_curve = FALSE,
  rel_heights = c(1, 0.7)
)
```

## Arguments

- diag:

  An object of class `"r4pde_curve_diagnostics"`.

- layout:

  Layout of diagnostic panels: `"2x2"`, `"vertical"`, or `"horizontal"`.

- show_curve:

  Logical; whether to include the example adjusted curve.

- rel_heights:

  Relative heights for panels when stacking.

## Value

A composite plot object generated using cowplot.
