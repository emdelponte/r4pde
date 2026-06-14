# \#' Simulate Disease Suppression Profile Data

Simulates a dataset for demonstrating functional analysis with Disease
Suppression Profiles. Returns a data frame with varying treatments and
replicates over multiple time points.

## Usage

``` r
simulate_dsp_data()
```

## Value

A `tibble` with columns `treatment`, `rep`, `time`, `y_control`,
`dsp_true`, `severity_true`, and `severity`.
