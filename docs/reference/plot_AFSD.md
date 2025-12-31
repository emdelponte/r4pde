# Plot ASFD

This function creates a tile plot of the foci (cluster) identified by
the AFSD function. It colors each cell in a foci and labels the centroid
of each cluster with the foci ID. The 'ggplot2' package is used for the
plot, and will be automatically installed if not already present.

## Usage

``` r
plot_AFSD(df)
```

## Arguments

- df:

  A dataframe containing at least three columns: 'x', 'y', and
  'cluster_id'. 'x' and 'y' are spatial coordinates and 'cluster_id' is
  the cluster identifier to which each cell belongs.

## Value

A ggplot object with the scatter plot of foci (clusters).

## See also

Other Spatial analysis:
[`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md),
[`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md),
[`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md),
[`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md),
[`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md)

## Examples

``` r
df <- data.frame(x = sample(1:100, 500, replace = TRUE),
                 y = sample(1:100, 500, replace = TRUE),
                 i = sample(0:1, 500, replace = TRUE, prob = c(0.7, 0.3)))

# Perform the AFSD
result <- AFSD(df)
# Plot the foci
plot_AFSD(result[[3]])

```
