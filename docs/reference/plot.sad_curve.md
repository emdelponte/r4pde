# Plot curves from sad_curve()

Plot curves from sad_curve()

## Usage

``` r
# S3 method for class 'sad_curve'
plot(x, what = c("bias", "expected_abs_error", "difference"), ...)
```

## Arguments

- x:

  An object returned by
  [`sad_curve()`](https://emdelponte.github.io/r4pde/reference/sad_curve.md).

- what:

  `"bias"` (bias curves with pointwise intervals),
  `"expected_abs_error"` (only for `type = "linear"`) or `"difference"`
  (differences with simultaneous bands; shaded ranges are significant).

- ...:

  Not used.

## Value

A ggplot object.
