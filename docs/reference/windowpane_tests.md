# Windowpane Tests for Correlation Analysis

This function performs bootstrapped correlation analysis for multiple
predictors against a response variable. It applies the Simes method for
global significance testing and calculates individual correlations,
p-values, and bootstrap statistics.

## Usage

``` r
windowpane_tests(
  data,
  response_var,
  corr_type = "spearman",
  R = 1000,
  global_alpha = 0.05,
  individual_alpha = 0.005
)
```

## Arguments

- data:

  A data frame containing the predictors and the response variable.

- response_var:

  A string representing the name of the response variable in the data
  frame.

- corr_type:

  A string specifying the correlation method to use; options are
  "spearman" (default), "pearson", or "kendall".

- R:

  An integer indicating the number of bootstrap replications. Default is
  1000.

- global_alpha:

  A numeric value representing the global alpha level for the Simes
  correction. Default is 0.05.

- individual_alpha:

  A numeric value for the individual alpha threshold for testing
  individual predictors. Default is 0.005.

## Value

A list containing the following elements:

- results:

  A data frame with columns: `variable`, `correlation`, `p_value`,
  `mean_corr`, `sd_corr`, `median_corr`, `rank`, `simes_threshold`,
  `significant_simes`, and `individual_significant`.

- summary_table:

  A data frame summarizing the global p-value (*Pg*) and maximum
  correlation.

- global_significant:

  A logical value indicating whether the global test is significant.

## Details

The function calculates correlations between the response variable and
each predictor in the data frame, using bootstrapping to generate mean,
standard deviation, and median estimates of the correlation. The Simes
correction is applied to control for multiple testing, providing a
global p-value (*Pg*). The function also returns the maximum observed
correlation.
