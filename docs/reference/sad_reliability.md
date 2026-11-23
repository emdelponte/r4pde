# Inter-rater reliability by condition

Computes the intraclass correlation coefficient ICC(A,1) with its
confidence interval and the overall concordance correlation coefficient
(OCC) among raters for each condition, as in Table 3 of Franceschi et
al. (2020).

## Usage

``` r
sad_reliability(x, conf = 0.95)
```

## Arguments

- x:

  An object returned by
  [`sad_data()`](https://emdelponte.github.io/r4pde/reference/sad_data.md).

- conf:

  Confidence level for the ICC interval.

## Value

A tibble with one row per condition and columns `ICC`, `ICC_lower`,
`ICC_upper`, `OCC`, `n_raters` and `n_images`.

## Details

ICC(A,1) is the two-way random effects, absolute agreement, single rater
coefficient of McGraw and Wong (1996). OCC is the overall concordance
correlation coefficient of Barnhart et al. (2002), computed as in
[`epiR::epi.occc()`](https://rdrr.io/pkg/epiR/man/epi.occc.html). Only
raters that assessed all images of a condition are used.

## References

McGraw, K. O., and Wong, S. P. (1996). Forming inferences about some
intraclass correlation coefficients. Psychological Methods 1:30-46.

Barnhart, H. X., Haber, M., and Song, J. (2002). Overall concordance
correlation coefficient for evaluating agreement among multiple
observers. Biometrics 58:1020-1027.

## Examples

``` r
x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
sad_reliability(x)
#> # A tibble: 4 × 7
#>   condition   ICC ICC_lower ICC_upper   OCC n_raters n_images
#>   <fct>     <dbl>     <dbl>     <dbl> <dbl>    <int>    <int>
#> 1 UN_Old    0.749     0.650     0.833 0.747       20       50
#> 2 UN_New    0.739     0.644     0.823 0.736       20       50
#> 3 SAD_Old   0.762     0.673     0.840 0.760       20       50
#> 4 SAD_New   0.942     0.918     0.962 0.941       20       50
```
