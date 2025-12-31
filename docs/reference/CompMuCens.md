# Survival analysis for quantitative ordinal scale data.

Survival analysis for quantitative ordinal scale data.

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

Returns a list containing the score statistic, hypothesis tests,
adjusted significance level, and conclusion based on pairwise
comparisons.

## Details

To assist plant pathologists in analyzing quantitative ordinal scale
data and encourage the uptake of the interval-censored analysis method,
Chiang and collaborators have developed this function and provided
comprehensive explanation of the program code used to implement class
ratings analyzed through this method in this repository:
https://github.com/StatisticalMethodsinPlantProtection/CompMuCens
According to results in the paper, the method can be applied to reduce
the risk of type II errors when considering quantitative ordinal data,
which are widely used in plant pathology and related disciplines.The
function starts by converting the data into a censored data format and
performs multiple pairwise comparisons to determine significance using
the score statistic method.

## References

Chiang, K.S., Chang, Y.M., Liu, H.I., Lee, J.Y., El Jarroudi, M. and
Bock, C., 2023. Survival Analysis as a Basis to Test Hypotheses When
Using Quantitative Ordinal Scale Disease Severity Data. Phytopathology,
in press. Available at:
https://apsjournals.apsnet.org/doi/abs/10.1094/PHYTO-02-23-0055-R

## See also

Other Disease quantification:
[`DSI()`](https://emdelponte.github.io/r4pde/reference/DSI.md),
[`DSI2()`](https://emdelponte.github.io/r4pde/reference/DSI2.md)

## Examples

``` r
# Entering your data as ordinal rating scores
trAs=c(5,4,2,5,5,4,4,2,5,2,2,3,4,3,2,2,6,2,2,4,2,4,2,4,5,3,4,2,2,3)
trBs=c(5,3,2,4,4,5,4,5,4,4,6,4,5,5,5,2,6,2,3,5,2,6,4,3,2,5,3,5,4,5)
trCs=c(2,3,1,4,1,1,4,1,1,3,2,1,4,1,1,2,5,2,1,3,1,4,2,2,2,4,2,3,2,2)
trDs=c(5,5,4,5,5,6,6,4,6,4,3,5,5,6,4,6,5,6,5,4,5,5,5,3,5,6,5,5,5,6)
# Data shaping into input format
inputData = data.frame(treatment=c(rep("A",30),rep("B",30),rep("C",30),
rep("D",30)), x=c(trAs, trBs, trCs, trDs))
# Perform analysis using CompMuCens() function
CompMuCens(dat=inputData, scale=c(0,3,6,12,25,50,75,88,94,97,100,100),ckData=TRUE)
#> $inputData
#>     treatment ClassValue intervals
#> 1           A          5  [12, 25]
#> 2           A          4  [ 6, 12]
#> 3           A          2  [ 0,  3]
#> 4           A          5  [12, 25]
#> 5           A          5  [12, 25]
#> 6           A          4  [ 6, 12]
#> 7           A          4  [ 6, 12]
#> 8           A          2  [ 0,  3]
#> 9           A          5  [12, 25]
#> 10          A          2  [ 0,  3]
#> 11          A          2  [ 0,  3]
#> 12          A          3  [ 3,  6]
#> 13          A          4  [ 6, 12]
#> 14          A          3  [ 3,  6]
#> 15          A          2  [ 0,  3]
#> 16          A          2  [ 0,  3]
#> 17          A          6  [25, 50]
#> 18          A          2  [ 0,  3]
#> 19          A          2  [ 0,  3]
#> 20          A          4  [ 6, 12]
#> 21          A          2  [ 0,  3]
#> 22          A          4  [ 6, 12]
#> 23          A          2  [ 0,  3]
#> 24          A          4  [ 6, 12]
#> 25          A          5  [12, 25]
#> 26          A          3  [ 3,  6]
#> 27          A          4  [ 6, 12]
#> 28          A          2  [ 0,  3]
#> 29          A          2  [ 0,  3]
#> 30          A          3  [ 3,  6]
#> 31          B          5  [12, 25]
#> 32          B          3  [ 3,  6]
#> 33          B          2  [ 0,  3]
#> 34          B          4  [ 6, 12]
#> 35          B          4  [ 6, 12]
#> 36          B          5  [12, 25]
#> 37          B          4  [ 6, 12]
#> 38          B          5  [12, 25]
#> 39          B          4  [ 6, 12]
#> 40          B          4  [ 6, 12]
#> 41          B          6  [25, 50]
#> 42          B          4  [ 6, 12]
#> 43          B          5  [12, 25]
#> 44          B          5  [12, 25]
#> 45          B          5  [12, 25]
#> 46          B          2  [ 0,  3]
#> 47          B          6  [25, 50]
#> 48          B          2  [ 0,  3]
#> 49          B          3  [ 3,  6]
#> 50          B          5  [12, 25]
#> 51          B          2  [ 0,  3]
#> 52          B          6  [25, 50]
#> 53          B          4  [ 6, 12]
#> 54          B          3  [ 3,  6]
#> 55          B          2  [ 0,  3]
#> 56          B          5  [12, 25]
#> 57          B          3  [ 3,  6]
#> 58          B          5  [12, 25]
#> 59          B          4  [ 6, 12]
#> 60          B          5  [12, 25]
#> 61          C          2  [ 0,  3]
#> 62          C          3  [ 3,  6]
#> 63          C          1         0
#> 64          C          4  [ 6, 12]
#> 65          C          1         0
#> 66          C          1         0
#> 67          C          4  [ 6, 12]
#> 68          C          1         0
#> 69          C          1         0
#> 70          C          3  [ 3,  6]
#> 71          C          2  [ 0,  3]
#> 72          C          1         0
#> 73          C          4  [ 6, 12]
#> 74          C          1         0
#> 75          C          1         0
#> 76          C          2  [ 0,  3]
#> 77          C          5  [12, 25]
#> 78          C          2  [ 0,  3]
#> 79          C          1         0
#> 80          C          3  [ 3,  6]
#> 81          C          1         0
#> 82          C          4  [ 6, 12]
#> 83          C          2  [ 0,  3]
#> 84          C          2  [ 0,  3]
#> 85          C          2  [ 0,  3]
#> 86          C          4  [ 6, 12]
#> 87          C          2  [ 0,  3]
#> 88          C          3  [ 3,  6]
#> 89          C          2  [ 0,  3]
#> 90          C          2  [ 0,  3]
#> 91          D          5  [12, 25]
#> 92          D          5  [12, 25]
#> 93          D          4  [ 6, 12]
#> 94          D          5  [12, 25]
#> 95          D          5  [12, 25]
#> 96          D          6  [25, 50]
#> 97          D          6  [25, 50]
#> 98          D          4  [ 6, 12]
#> 99          D          6  [25, 50]
#> 100         D          4  [ 6, 12]
#> 101         D          3  [ 3,  6]
#> 102         D          5  [12, 25]
#> 103         D          5  [12, 25]
#> 104         D          6  [25, 50]
#> 105         D          4  [ 6, 12]
#> 106         D          6  [25, 50]
#> 107         D          5  [12, 25]
#> 108         D          6  [25, 50]
#> 109         D          5  [12, 25]
#> 110         D          4  [ 6, 12]
#> 111         D          5  [12, 25]
#> 112         D          5  [12, 25]
#> 113         D          5  [12, 25]
#> 114         D          3  [ 3,  6]
#> 115         D          5  [12, 25]
#> 116         D          6  [25, 50]
#> 117         D          5  [12, 25]
#> 118         D          5  [12, 25]
#> 119         D          5  [12, 25]
#> 120         D          6  [25, 50]
#> 
#> $U.Score
#> # A tibble: 4 × 2
#>   treatment  score
#>   <chr>      <dbl>
#> 1 D         -15.1 
#> 2 B          -4.54
#> 3 A           4.23
#> 4 C          15.4 
#> 
#> $Hypothesis.test
#> # A tibble: 3 × 4
#>   treat1 treat2 `p-value for H0: treat1 <= treat2` p-value for H0: treat1 = tr…¹
#>   <chr>  <chr>                               <dbl> <lgl>                        
#> 1 D      B                                0.00188  NA                           
#> 2 B      A                                0.0113   NA                           
#> 3 A      C                                0.000746 NA                           
#> # ℹ abbreviated name: ¹​`p-value for H0: treat1 = treat2`
#> 
#> $adj.Signif
#> [1] 0.01666667
#> 
#> $Conclusion
#> [1] "c(\"D>B\", \">A\", \">C\")"
#> 
```
