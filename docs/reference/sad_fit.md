# Fit the one-stage hierarchical model for SAD validation

Fits a single mixed model to all severity estimates of a SAD validation
study, instead of first summarising each rater by Lin's statistics. The
model describes how the mean estimate, the variation among raters and
the residual variation change along the severity range in each
condition.

## Usage

``` r
sad_fit(
  x,
  structure = c("auto", "full", "noim", "sep", "base"),
  dispersion = c("quad", "lin"),
  family = c("gaussian", "ordbeta"),
  center = 10,
  scale = 10
)
```

## Arguments

- x:

  An object returned by
  [`sad_data()`](https://emdelponte.github.io/r4pde/reference/sad_data.md).

- structure:

  Random-effects structure for raters. `"auto"` (default) tries
  `"full"`, `"noim"`, `"sep"` and `"base"` in this order and keeps the
  first that converges with a positive definite Hessian and finite
  standard errors. See Details.

- dispersion:

  Model for the residual variance: `"quad"` (default), a quadratic
  function of severity in each condition, which allows larger errors at
  intermediate severity; or `"lin"`, a log-linear function.

- family:

  `"gaussian"` (estimates in %) or `"ordbeta"` (ordered beta regression
  on the 0-1 scale, which handles the bounds at 0 and 100% and accepts
  exact zeros and ones). `"ordbeta"` is useful as a sensitivity
  analysis.

- center, scale:

  Severity is entered as `z = (actual - center) / scale`. The defaults
  (10% and 10 percentage points) give intercepts that refer to 10%
  severity and help convergence.

## Value

An object of class `sad_fit`: a list with the `model`, the `structure`
used, the `data` (with `z`), `center`, `scale`, `family`, `dispersion`,
`conditions` and `methods`.

## Details

The mean model is `estimate ~ z * condition`, so each condition has its
own intercept (bias at `center`) and slope (calibration). Random-effects
structures:

- full:

  `(0 + method + method:z | rater) + (1 | image) + (1 | image:condition)`:
  a random intercept and slope per rater and method, all correlated, an
  image effect (intrinsic difficulty) and an image x condition effect.

- noim:

  as `full` without the image x condition effect.

- sep:

  `(0 + method | rater) + (0 + method:z | rater) + (1 | image)`:
  intercepts and slopes uncorrelated.

- base:

  `(1 + z | rater) + (1 | image)`: rater effects common to all methods.

The dispersion formula is `~ condition * (z + I(z^2))` or
`~ condition * z`. Models are fitted with
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html).

## See also

[`sad_curve()`](https://emdelponte.github.io/r4pde/reference/sad_curve.md),
[`sad_raters()`](https://emdelponte.github.io/r4pde/reference/sad_raters.md),
[`sad_images()`](https://emdelponte.github.io/r4pde/reference/sad_images.md),
[`sad_compare()`](https://emdelponte.github.io/r4pde/reference/sad_compare.md)

## Examples

``` r
# \donttest{
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
  fit <- sad_fit(x)
  fit
}
#> Structure 'full' did not converge; trying the next one.
#> Structure 'noim' did not converge; trying the next one.
#> One-stage model for SAD validation
#>   family: gaussian | dispersion: quad | rater structure: sep
#>   severity centred at 10%, in units of 10 percentage points
#>   40 raters, 50 images, 4000 estimates; conditions: UN_Old, UN_New, SAD_Old, SAD_New
#>   AIC = 28896.6
# }
```
