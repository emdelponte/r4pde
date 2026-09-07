# Plot functional suppression profiles

Plot functional suppression profiles

## Usage

``` r
# S3 method for class 'functional_suppression_profiles'
plot(
  x,
  type = c("dendrogram", "profiles", "heatmap", "rank", "all"),
  show_cut = TRUE,
  show_points = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class `"functional_suppression_profiles"`.

- type:

  Character string specifying the plot type. One of `"dendrogram"`,
  `"profiles"`, `"heatmap"`, `"rank"`, or `"all"`.

- show_cut:

  Logical; whether to display the cluster cut height in dendrogram.
  Default is `TRUE`.

- show_points:

  Logical; whether to overlay observed DSP values as points in the
  profiles plot. Default is `FALSE`.

- ...:

  Additional arguments passed to specific plot functions.

## Value

A `ggplot` object or a list of `ggplot` objects.

## Details

The plot method visualizes the components of a Functional Suppression
Profile. The dendrogram defines the Functional Suppression Profiles
based on temporal similarity. The heatmap is ordered by the dendrogram
to help interpret these functional groups using summary metrics.
