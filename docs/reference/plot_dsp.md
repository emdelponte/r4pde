# Plot Disease Suppression Profiles

Plots functional contrast curves (DSP) over time.

## Usage

``` r
plot_dsp(
  object,
  time = NULL,
  response = NULL,
  treatment = NULL,
  group = NULL,
  facet = NULL,
  show_zero = TRUE,
  linewidth = 1,
  alpha = 1,
  ylab = "Disease suppression profile",
  xlab = "Time"
)
```

## Arguments

- object:

  An object of class `"functional_dsp"`.

- time:

  Character string for the time column (auto-detected if object is
  `"functional_dsp"`).

- response:

  Character string for the DSP column.

- treatment:

  Character string for the treatment column.

- group:

  Optional character vector of grouping variables.

- facet:

  Optional character string specifying a variable for faceting.

- show_zero:

  Logical; whether to draw a horizontal reference line at y = 0. Default
  is `TRUE`.

- linewidth:

  Numeric; width of the lines.

- alpha:

  Numeric; transparency of the lines.

- ylab:

  Character string; label for the y-axis.

- xlab:

  Character string; label for the x-axis.

## Value

A `ggplot` object.
