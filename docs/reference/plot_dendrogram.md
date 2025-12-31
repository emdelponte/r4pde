# Plot functional dendrogram of epidemic curves

Visualizes the hierarchical clustering of treatments based on functional
distances among epidemic curves.

## Usage

``` r
plot_dendrogram(x, label_fun = NULL, palette = NULL, show_cut = TRUE)
```

## Arguments

- x:

  An object of class `"r4pde_compare_curves"`.

- label_fun:

  Optional function to modify treatment labels.

- palette:

  Optional named vector of colors for clusters.

- show_cut:

  Logical; whether to display the cluster cut height.

## Value

A `ggplot` object.
