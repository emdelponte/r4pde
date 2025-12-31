# Plot environment-adjusted epidemic curves by cluster

Plots environment-adjusted mean epidemic curves for each treatment,
colored according to functional cluster membership.

The returned object is a `ggplot` and can be further modified using
standard ggplot2 layers.

## Usage

``` r
plot_curves(x, label_fun = NULL, palette = NULL, alpha = 0.9, linewidth = 1.1)
```

## Arguments

- x:

  An object of class `"r4pde_compare_curves"`.

- label_fun:

  Optional function to modify treatment labels.

- palette:

  Optional named vector of colors for clusters.

- alpha:

  Line transparency.

- linewidth:

  Line width.

## Value

A `ggplot` object.
