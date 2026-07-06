# Binary Power Law Analysis for Spatial Disease Patterns

This function calculates the Binary Power Law (BPL) parameters for
spatial disease patterns, fits a linear model, and performs a hypothesis
test for the slope.

## Usage

``` r
BPL(data)
```

## Arguments

- data:

  A data frame containing the following columns:

  - `field`: The field identifier.

  - `n`: The number of observations in each quadrat.

  - `i`: The incidence count in each quadrat.

## Value

A list containing the following elements:

- `summary`: A data frame summarizing the input data by field, including
  total observations (`n_total`), mean incidence (`incidence_mean`),
  observed variance (`V`), and binomial variance (`Vbin`).

- `model_summary`: A summary of the linear model fitted to the
  log-transformed variances.

- `hypothesis_test`: The result of the hypothesis test for the slope
  being equal to 1.

- `ln_Ap`: The intercept of the linear model, representing the natural
  logarithm of the parameter \\ A_p \\.

- `slope`: The slope of the linear model.

## Details

The function performs the following steps:

1.  Summarizes the data by field to calculate the total number of
    observations (`n_total`), mean incidence (`incidence_mean`),
    observed variance (`V`), and binomial variance (`Vbin`).

2.  Log-transforms the variances.

3.  Fits a linear model to the log-transformed variances.

4.  Tests the hypothesis that the slope of the linear model is equal to
    1.

## See also

Other Spatial analysis:
[`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md),
[`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md),
[`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md),
[`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)

## Examples

``` r
# \donttest{
# Example usage with a sample data frame
result <- BPL(FHBWheat)
print(result$summary)
#> # A tibble: 153 × 7
#>    field n_total incidence_mean      V   Vbin log_V log_Vbin
#>    <int>   <int>          <dbl>  <dbl>  <dbl> <dbl>    <dbl>
#>  1     1     200          0.52  0.0333 0.0250 -3.40    -3.69
#>  2     2     200          0.535 0.0529 0.0249 -2.94    -3.69
#>  3     3     200          0.3   0.0274 0.021  -3.60    -3.86
#>  4     4     200          0.42  0.0206 0.0244 -3.88    -3.71
#>  5     5     200          0.35  0.01   0.0228 -4.61    -3.78
#>  6     6     200          0.365 0.0466 0.0232 -3.07    -3.76
#>  7     7     200          0.315 0.0245 0.0216 -3.71    -3.84
#>  8     8     200          0.505 0.0352 0.0250 -3.35    -3.69
#>  9     9     200          0.56  0.0204 0.0246 -3.89    -3.70
#> 10    10     200          0.43  0.0243 0.0245 -3.72    -3.71
#> # ℹ 143 more rows
print(result$model_summary)
#> 
#> Call:
#> lm(formula = log_V ~ log_Vbin, data = df_summary)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.06015 -0.23599 -0.00801  0.25043  1.80832 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   0.3158     0.1860   1.698   0.0916 .  
#> log_Vbin      1.0205     0.0420  24.295   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 0.3935 on 151 degrees of freedom
#> Multiple R-squared:  0.7963, Adjusted R-squared:  0.7949 
#> F-statistic: 590.3 on 1 and 151 DF,  p-value: < 2.2e-16
#> 
print(result$hypothesis_test)
#> 
#> Linear hypothesis test:
#> log_Vbin = 1
#> 
#> Model 1: restricted model
#> Model 2: log_V ~ log_Vbin
#> 
#>   Res.Df    RSS Df Sum of Sq      F Pr(>F)
#> 1    152 23.422                           
#> 2    151 23.385  1  0.036968 0.2387 0.6259
print(paste("ln(Ap):", result$ln_Ap))
#> [1] "ln(Ap): 0.315806312326611"
print(paste("Slope (b):", result$slope))
#> [1] "Slope (b): 1.02052263172857"
# }
```
