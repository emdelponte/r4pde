# Prepare severity estimates for SAD validation

Organises visual severity estimates made with and without a standard
area diagram (SAD), or with two or more SADs, into the long format used
by the `sad_*` family, and checks the data for common problems.

## Usage

``` r
sad_data(
  data,
  rater = "rater",
  image = "leaf",
  actual = "actual",
  estimate = "estimate",
  method = "method",
  wide = NULL,
  group = NULL,
  sep = "_"
)
```

## Arguments

- data:

  A data frame in long or wide format.

- rater, image, actual:

  Names of the columns with the rater identifier, the image (leaf,
  fruit, ...) identifier and the actual severity (%).

- estimate:

  Name of the column with the estimates (long format).

- method:

  Name of the column with the assessment method (long format), for
  example `"unaided"` and `"aided"`.

- wide:

  Character vector with the estimate columns in wide format. Each column
  becomes one level of `method`. If the vector is named, the names are
  used as the method labels, e.g. `c(UN = "unaided", SAD = "aided")`.

- group:

  Optional name of a column that splits raters into groups, for example
  which SAD each rater used. When given, `condition` is `method` and
  `group` pasted together.

- sep:

  Separator used to build `condition` from `method` and `group`.

## Value

A tibble of class `sad_data` with columns `rater`, `image`, `actual`,
`estimate`, `method`, `group`, `condition` and `error` (estimate minus
actual, in percentage points), and attribute `design`.

## Details

Each combination of method and group is a *condition*. The function
warns when an image has more than one actual severity value, when
estimates fall outside 0 to 100%, when a rater x image x condition
combination is repeated, and when a rater lacks one of the methods. The
design is labelled `"paired"` when raters assessed the images under more
than one condition (the usual unaided then aided study), and
`"independent"` otherwise.

## See also

[`sad_ccc()`](https://emdelponte.github.io/r4pde/reference/sad_ccc.md),
[`sad_compare()`](https://emdelponte.github.io/r4pde/reference/sad_compare.md),
[`sad_fit()`](https://emdelponte.github.io/r4pde/reference/sad_fit.md)

## Examples

``` r
x <- sad_data(SADRustSoybean, rater = "rater", image = "leaf",
              wide = c(UN = "unaided", SAD = "aided"), group = "sad")
x
#> SAD validation data (paired design): 40 raters, 50 images, 4000 estimates
#> Conditions: UN_Old, UN_New, SAD_Old, SAD_New 
#> # A tibble: 4,000 × 8
#>    rater image actual estimate method group condition  error
#>    <fct> <fct>  <dbl>    <dbl> <fct>  <fct> <fct>      <dbl>
#>  1 1     1       0.25      5   UN     Old   UN_Old     4.75 
#>  2 1     1       0.25      0.6 SAD    Old   SAD_Old    0.35 
#>  3 1     2       2.5      15   UN     Old   UN_Old    12.5  
#>  4 1     2       2.5       5   SAD    Old   SAD_Old    2.5  
#>  5 1     3       7.24     10   UN     Old   UN_Old     2.76 
#>  6 1     3       7.24      7   SAD    Old   SAD_Old   -0.240
#>  7 1     4       7.31     20   UN     Old   UN_Old    12.7  
#>  8 1     4       7.31     10   SAD    Old   SAD_Old    2.69 
#>  9 1     5       9.07     15   UN     Old   UN_Old     5.93 
#> 10 1     5       9.07     18   SAD    Old   SAD_Old    8.93 
#> # ℹ 3,990 more rows
```
