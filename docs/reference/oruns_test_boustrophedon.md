# Boustrophedon Run Test for Binary Matrix

Applies the ordinary runs test to a binary matrix using
boustrophedon-style traversal. The function supports two modes: row-wise
and column-wise boustrophedon. Each traversal flattens the matrix into a
1D sequence which is then tested using `oruns_test`.

## Usage

``` r
oruns_test_boustrophedon(mat)
```

## Arguments

- mat:

  A binary matrix (containing 0s and 1s, and possibly NAs).

## Value

A list with two elements:

- rowwise_boustrophedon:

  List containing the sequence and result of `oruns_test` for row-wise
  traversal.

- colwise_boustrophedon:

  List containing the sequence and result of `oruns_test` for
  column-wise traversal.

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
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)
