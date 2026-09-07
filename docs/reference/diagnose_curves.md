# Diagnostic tools for functional epidemic curve models

Diagnostic tools for functional epidemic curve models fitted with
[`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md).
Produces diagnostic plots without invoking base graphics.

## Usage

``` r
diagnose_curves(x, grid_n = 200)
```

## Arguments

- x:

  An object of class `"functional_curves"`.

- grid_n:

  Number of points used for the diagnostic smooth curve.

## Value

An object of class `"r4pde_curve_diagnostics"` containing:

- diagnostic data,

- k-index checks,

- concurvity estimates,

- ggplot-based diagnostic panels.
