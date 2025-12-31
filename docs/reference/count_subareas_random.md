# Random Subgrid Sampling of a Binary Matrix

Randomly samples submatrices (quadrats) of specified size from a binary
matrix, and returns the positions, submatrices, and count of 1s in each
sampled quadrat.

## Usage

``` r
count_subareas_random(matrix_data, sub_rows = 3, sub_cols = 3, n_samples = 100)
```

## Arguments

- matrix_data:

  A binary matrix of 0s and 1s.

- sub_rows:

  Number of rows in each subgrid sample.

- sub_cols:

  Number of columns in each subgrid sample.

- n_samples:

  Number of subgrid samples to draw.

## Value

A list of sampled subgrids. Each element is a list with:

- position:

  Row and column start position of the sample.

- submatrix:

  The sampled subgrid matrix.

- count:

  Number of 1s in the sampled submatrix.

## See also

Other Spatial analysis:
[`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md),
[`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md),
[`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md),
[`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)
