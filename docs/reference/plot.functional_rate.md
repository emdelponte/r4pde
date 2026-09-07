# Plot functional_rate curves and rates

Plot functional_rate curves and rates

## Usage

``` r
# S3 method for class 'functional_rate'
plot(
  x,
  which = NULL,
  type = c("both", "curve", "rate"),
  show_ci = TRUE,
  show_threshold = TRUE,
  ncol = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `"functional_rate"`.

- which:

  Optional character vector specifying which treatments/groups to
  display.

- type:

  Character specifying what to plot: `"both"` (default; fitted curve in
  top panel and derivative in bottom panel), `"curve"` (only
  trajectory), or `"rate"` (only instantaneous rate).

- show_ci:

  Logical; whether to show confidence ribbons. Default is `TRUE`.

- show_threshold:

  Logical; whether to show threshold and zero horizontal lines. Default
  is `TRUE`.

- ncol:

  Optional integer specifying number of facet columns.

- ...:

  Additional arguments.

## Value

A `ggplot` or combined `patchwork` object.
