# Simulate power for planning a SAD validation study

Simulates paired studies (each rater assesses the same images without
and with a SAD) and compares the power of the conventional two-stage
approach, a paired t test on \\\rho_c\\, with the likelihood ratio test
of the one-stage model for any effect of the SAD on mean or dispersion
(4 df).

## Usage

``` r
sad_simulate(
  n_rater = c(10, 20, 40),
  n_image = c(20, 40),
  nsim = 200,
  shift = -3,
  slope_change = 0,
  precision = 1,
  bias0 = 3,
  slope0 = 1.1,
  sd_rater = 5,
  sd_slope = 0.18,
  sd_image = 3.5,
  sigma0 = 8,
  rho_rater = 0.5,
  severity = function(n) stats::runif(n, 0.5, 80),
  cores = 1,
  seed = 1
)
```

## Arguments

- n_rater, n_image:

  Vectors with the numbers of raters and images; all combinations are
  simulated.

- nsim:

  Number of simulated studies per scenario.

- shift:

  Change in bias at 10% severity caused by the SAD (percentage points).

- slope_change:

  Change in slope caused by the SAD.

- precision:

  Multiplier of the rater and residual standard deviations with the SAD
  (\< 1 means the SAD makes estimates more precise).

- bias0, slope0:

  Bias at 10% severity and slope without the SAD.

- sd_rater, sd_slope:

  Standard deviations of rater intercepts and slopes.

- sd_image:

  Standard deviation of image difficulty.

- sigma0:

  Residual standard deviation at 10% severity without the SAD.

- rho_rater:

  Correlation of a rater's effects with and without the SAD.

- severity:

  Function of `n` that draws the actual severities.

- cores:

  Number of cores
  ([`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html);
  ignored on Windows).

- seed:

  Random seed; the user's random number stream is restored.

## Value

A tibble with one row per scenario and the rejection rates at 5%
(`power_two_stage`, `power_one_stage`) and the number of converged fits.

## Details

Default parameters are typical values estimated from ten published SAD
validation studies. Estimates are truncated to 0-100%. The one-stage
model uses the `"sep"` structure. With fewer than about 20 raters the
likelihood ratio test can be liberal; a parametric bootstrap is then
safer.

## Examples

``` r
# \donttest{
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  sad_simulate(n_rater = 10, n_image = 20, nsim = 5)
}
#> # A tibble: 1 × 5
#>   n_rater n_image power_two_stage power_one_stage  n_ok
#>     <dbl>   <dbl>           <dbl>           <dbl> <int>
#> 1      10      20               0               0     5
# }
```
