# Runs Test

Perform a runs test on the input data to test for clustering or
randomness.

## Usage

``` r
oruns_test(x)
```

## Arguments

- x:

  A numeric vector representing the input data

## Value

an `r4pde.oruns` object.

An `r4pde.oruns` object is a `list` containing:

- U, number of runs,

- EU, expected number of runs,

- sU, standard deviation of the expected number of runs

- Z, Z-score of the observed number of runs,

- pvalue, the p-value of the Z-score, and

- result, the test result of either "aggregation or clustering" or
  "randomness"

## See also

Other Spatial analysis:
[`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md),
[`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md),
[`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md),
[`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)

## Examples

``` r
oruns_test(c(1, 0, 1, 1, 0, 1, 0, 0, 1, 1))
#> Ordinary Runs Test of Data Sequence:
#> -------------------------------------
#> Total Number of Runs (U): 7
#> Expected Number of Runs (EU): 5.80
#> Standard Deviation of Runs (sU): 1.42
#> Z-score: 0.84
#> P-value: 0.1996
#> 
#> Interpretation:
#> Based on the Z-score, the sequence exhibits 'randomness'.
```
