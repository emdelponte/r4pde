# Functional Contrast (Disease Suppression Profiles)

Calculates functional contrast curves (Disease Suppression Profiles -
DSP) between a reference treatment (e.g., untreated control) and other
fungicide treatments over time. The DSP is defined by default as:
\$\$DSP(t) = y\_{control}(t) - y\_{treatment}(t)\$\$ where \\y(t)\\ is
the disease intensity at time \\t\\. The goal is to evaluate the
intensity, timing, and persistence of the fungicide protection.

## Usage

``` r
functional_contrast(
  data,
  time = "time",
  response = "response",
  treatment = "treatment",
  reference = "Control",
  group = NULL,
  smooth = FALSE,
  grid = NULL,
  contrast = c("difference", "relative"),
  keep_reference = FALSE
)
```

## Arguments

- data:

  A data frame containing disease progress observations, or an object of
  class `"functional_curves"`.

- time:

  Character string naming the time variable.

- response:

  Character string naming the response variable.

- treatment:

  Character string naming the treatment variable.

- reference:

  Character string specifying the reference treatment level. Default is
  `"Control"`.

- group:

  Optional character vector of grouping variables (e.g.,
  `c("experiment", "location")`).

- smooth:

  Logical; if `TRUE`, the function could apply smoothing (currently
  placeholder).

- grid:

  Optional numeric vector of common time points for interpolation. If
  provided, curves are interpolated.

- contrast:

  Character specifying the contrast type. Currently `"difference"` is
  implemented.

- keep_reference:

  Logical; if `TRUE`, includes the reference treatment in the output
  (with DSP = 0). Default is `FALSE`.

## Value

An object of class `"functional_dsp"` (inheriting from `"tbl_df"`),
containing:

- grouping variables (if any)

- treatment variable

- time variable

- `y_reference`: disease intensity of the reference

- `y_treatment`: disease intensity of the treatment

- `DSP`: the computed contrast

- `contrast_type`: type of contrast used

- `reference`: name of the reference treatment

## See also

[`functional_summary`](https://emdelponte.github.io/r4pde/reference/functional_summary.md),
[`plot_dsp`](https://emdelponte.github.io/r4pde/reference/plot_dsp.md)
