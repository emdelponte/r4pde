# Plot method for compare_curves objects

Produces a paired visualization of environment-adjusted mean curves and
the corresponding functional dendrogram.

## Usage

``` r
# S3 method for class 'r4pde_compare_curves'
plot(x, label_fun = NULL, palette = NULL, ...)
```

## Arguments

- x:

  An object of class `"r4pde_compare_curves"`.

- label_fun:

  Optional function to modify treatment labels.

- palette:

  Optional named vector of colors for clusters.

- ...:

  Additional arguments (currently ignored).

## Value

A patchwork object combining curves and dendrogram.
