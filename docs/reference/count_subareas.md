# Count the Number of Ones in Subareas of a Matrix

This function takes a binary matrix (0s and 1s) and divides it into
rectangular subareas, counting the number of ones in each. Subareas are
defined by the number of rows and columns specified by the user. If the
matrix dimensions are not perfectly divisible by the subarea size, edge
subareas may be smaller.

## Usage

``` r
count_subareas(matrix_data, sub_rows, sub_cols)
```

## Arguments

- matrix_data:

  A matrix of 0s and 1s to analyze.

- sub_rows:

  Number of rows in each subarea.

- sub_cols:

  Number of columns in each subarea.

## Value

A matrix where each cell corresponds to a subarea and contains the count
of ones.

## See also

Other Spatial analysis:
[`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md),
[`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md),
[`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md),
[`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)

## Examples

``` r
set.seed(123)
mat <- matrix(sample(c(0, 1), 12 * 16, replace = TRUE), nrow = 16, ncol = 12)
count_matrix <- count_subareas(mat, sub_rows = 3, sub_cols = 3)
print(count_matrix)
#>      [,1] [,2] [,3] [,4]
#> [1,]    3    4    3    7
#> [2,]    3    3    6    6
#> [3,]    3    4    6    2
#> [4,]    6    4    4    5
#> [5,]    3    6    7    4
#> [6,]    2    1    1    0
```
