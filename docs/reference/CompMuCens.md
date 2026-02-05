# Survival analysis for quantitative ordinal scale data

Survival analysis for quantitative ordinal scale data

## Usage

``` r
CompMuCens(dat, scale, grade = TRUE, ckData = FALSE)
```

## Arguments

- dat:

  Data frame containing the data to be processed.

- scale:

  A numeric vector indicating the scale or order of classes.

- grade:

  Logical. If TRUE, uses the class value. If FALSE, uses the NPE
  (Non-Parametric Estimate).

- ckData:

  Logical. If TRUE, returns the input data along with the results. If
  FALSE, returns only the results.

## Value

A list containing the score statistic, hypothesis tests, adjusted
significance level, and conclusion based on pairwise comparisons.

## Details

This function supports analysis of quantitative ordinal scale data via
interval-censored methods. The approach follows the workflow described
by Chiang et al. (2023) and the reference implementation provided in the
CompMuCens repository.

## References

Chiang, K.S., Chang, Y.M., Liu, H.I., Lee, J.Y., El Jarroudi, M. and
Bock, C. (2023). Survival Analysis as a Basis to Test Hypotheses When
Using Quantitative Ordinal Scale Disease Severity Data. Phytopathology.
<https://apsjournals.apsnet.org/doi/abs/10.1094/PHYTO-02-23-0055-R>

## See also

Other Disease quantification:
[`DSI()`](https://emdelponte.github.io/r4pde/reference/DSI.md),
[`DSI2()`](https://emdelponte.github.io/r4pde/reference/DSI2.md)

## Examples

``` r
if (requireNamespace("interval", quietly = TRUE)) {
  trAs <- c(5,4,2,5,5,4,4,2,5,2,2,3,4,3,2,2,6,2,2,4,2,4,2,4,5,3,4,2,2,3)
  trBs <- c(5,3,2,4,4,5,4,5,4,4,6,4,5,5,5,2,6,2,3,5,2,6,4,3,2,5,3,5,4,5)
  trCs <- c(2,3,1,4,1,1,4,1,1,3,2,1,4,1,1,2,5,2,1,3,1,4,2,2,2,4,2,3,2,2)
  trDs <- c(5,5,4,5,5,6,6,4,6,4,3,5,5,6,4,6,5,6,5,4,5,5,5,3,5,6,5,5,5,6)
  inputData <- data.frame(
    treatment = c(rep("A",30), rep("B",30), rep("C",30), rep("D",30)),
    x = c(trAs, trBs, trCs, trDs)
  )
  CompMuCens(dat = inputData,
             scale = c(0,3,6,12,25,50,75,88,94,97,100,100),
             ckData = TRUE)
}
```
