# Bias and calibration of each rater

Combines fixed and random effects of a
[`sad_fit()`](https://emdelponte.github.io/r4pde/reference/sad_fit.md)
model into the bias of each rater at the centring severity and the slope
of the rater's estimates on actual severity (1 = perfect calibration),
per condition.

## Usage

``` r
sad_raters(fit)

plot_sad_raters(x, ...)
```

## Arguments

- fit:

  An object returned by
  [`sad_fit()`](https://emdelponte.github.io/r4pde/reference/sad_fit.md)
  with `family = "gaussian"`.

- x:

  An object returned by `sad_raters()`.

- ...:

  Not used.

## Value

`sad_raters()` returns a tibble with `rater`, `condition`, `method`,
`bias_at_center` (percentage points) and `slope`. `plot_sad_raters()`
returns a ggplot object, with lines linking the same rater across
conditions.

## Examples

``` r
# \donttest{
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
  rt <- sad_raters(sad_fit(x))
  plot_sad_raters(rt)
}
#> Structure 'full' did not converge; trying the next one.
#> Structure 'noim' did not converge; trying the next one.

# }
```
