# Lin's concordance statistics per rater and condition

First stage of the conventional SAD validation: computes Lin's
concordance correlation coefficient and its components for each rater
under each condition.

## Usage

``` r
sad_ccc(x)
```

## Arguments

- x:

  An object returned by
  [`sad_data()`](https://emdelponte.github.io/r4pde/reference/sad_data.md).

## Value

A tibble with one row per rater and condition and columns `n`, `rho_c`,
`r`, `Cb`, `v`, `u`, `intercept` and `slope`.

## Details

\\\rho_c = r \times C_b\\, where \\r\\ is Pearson's correlation
(precision) and \\C_b = 2 / (v + 1/v + u^2)\\ is the bias correction
factor (accuracy). \\v = \sigma\_{est}/\sigma\_{act}\\ is the scale
shift and \\u = (\mu\_{est} - \mu\_{act}) /
\sqrt{\sigma\_{est}\sigma\_{act}}\\ is the location shift, negative when
the rater underestimates. Note that
[`epiR::epi.ccc()`](https://rdrr.io/pkg/epiR/man/epi.ccc.html) reports
the location shift with the opposite sign. The intercept and slope of
the regression of estimates on actual values are also returned.

## References

Lin, L. I. (1989). A concordance correlation coefficient to evaluate
reproducibility. Biometrics 45:255-268.

## See also

[`sad_compare()`](https://emdelponte.github.io/r4pde/reference/sad_compare.md)

## Examples

``` r
x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
head(sad_ccc(x))
#> # A tibble: 6 × 12
#>   rater group method condition     n rho_c     r    Cb     v       u intercept
#>   <fct> <fct> <fct>  <fct>     <int> <dbl> <dbl> <dbl> <dbl>   <dbl>     <dbl>
#> 1 1     Old   UN     UN_Old       50 0.901 0.942 0.956  1.05  0.297      7.22 
#> 2 2     Old   UN     UN_Old       50 0.906 0.923 0.982  1.18 -0.0998    -5.85 
#> 3 3     Old   UN     UN_Old       50 0.829 0.915 0.907  1.34  0.346      0.259
#> 4 4     Old   UN     UN_Old       50 0.923 0.927 0.996  1.09  0.0301     0.211
#> 5 5     Old   UN     UN_Old       50 0.771 0.916 0.841  1.36  0.531      4.35 
#> 6 6     Old   UN     UN_Old       50 0.869 0.897 0.969  1.18  0.191      2.31 
#> # ℹ 1 more variable: slope <dbl>
```
