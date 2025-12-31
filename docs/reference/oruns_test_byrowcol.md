# Runs Test for Each Row and Column of a Binary Matrix

Applies the ordinary runs test to each row and column of a binary matrix
individually.

## Usage

``` r
oruns_test_byrowcol(mat)
```

## Arguments

- mat:

  A binary matrix (containing 0s and 1s, and possibly NAs).

## Value

A list with four elements:

- row_results:

  Data frame with test results for each row.

- col_results:

  Data frame with test results for each column.

- row_summary:

  Percentage summary of interpretation for rows.

- col_summary:

  Percentage summary of interpretation for columns.

## See also

[`oruns_test`](https://emdelponte.github.io/r4pde/reference/oruns_test.md)

Other Spatial analysis:
[`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md),
[`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md),
[`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md),
[`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md),
[`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)
