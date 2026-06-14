# Summarize Disease Suppression Profiles

Extracts functional descriptors from Disease Suppression Profiles
(DSPs).

## Usage

``` r
functional_summary(
  object,
  time = "time",
  response = "DSP",
  treatment = "treatment",
  group = NULL,
  threshold = 0,
  positive_only = TRUE
)
```

## Arguments

- object:

  An object of class `"functional_dsp"` or a data frame with a DSP
  column.

- time:

  Character string naming the time column.

- response:

  Character string naming the DSP column.

- treatment:

  Character string naming the treatment column.

- group:

  Optional character vector of grouping variables.

- threshold:

  Numeric; minimum suppression threshold for calculating persistence.
  Default is 0.

- positive_only:

  Logical; if `TRUE`, negative DSP values are set to 0 before
  area/centroid calculations. Default is `TRUE`.

## Value

A `tibble` with one row per group/treatment containing:

- `protected_area`: Trapezoidal integral of DSP(t).

- `mean_suppression`: Temporal mean of DSP(t).

- `max_suppression`: Maximum value of DSP(t).

- `time_max_suppression`: Time at which maximum suppression occurs.

- `persistence`: Total duration where DSP(t) \> threshold.

- `centroid`: Temporal centroid defined as \\\int t \cdot DSP(t) dt /
  \int DSP(t) dt\\.

- `energy`: Integral of \\DSP(t)^2\\.

- `early_area, mid_area, late_area`: Protected area in the first,
  second, and third tertiles of the time period.

- `late_decline_slope`: Slope of DSP(t) in the late period.

- `n_time_points`: Number of time points evaluated.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using a generic data frame with custom column names
dat <- data.frame(
  DAA = c(0, 10, 20),
  fungicide = c("A", "A", "A"),
  contrast = c(0, 10, 0)
)
functional_summary(
  dat, 
  time = "DAA", 
  response = "contrast", 
  treatment = "fungicide"
)
} # }
```
