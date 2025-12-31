# Calculate the Disease severity Index (DSI) (frequency of each class)

This function calculates the Disease Severity Index (DSI) given a vector
of classes, a vector of frequencies, and a maximum possible class value.
The DSI is calculated as a weighted sum of class values, where each
class is multiplied by its corresponding frequency, then divided by the
product of the total frequency and maximum class value, and finally
multiplied by 100 to get a percentage.

## Usage

``` r
DSI2(class, freq, max)
```

## Arguments

- class:

  A numeric vector representing the classes.

- freq:

  A numeric vector representing the frequency of each class. Must be the
  same length as 'class'.

- max:

  A numeric value representing the maximum possible class value.

## Value

Returns a single numeric value representing the DSI.

## See also

Other Disease quantification:
[`CompMuCens()`](https://emdelponte.github.io/r4pde/reference/CompMuCens.md),
[`DSI()`](https://emdelponte.github.io/r4pde/reference/DSI.md)

## Examples

``` r
DSI2(c(0, 1, 2, 3, 4), c(2, 0, 5, 0, 5), 4)
#> [1] 62.5
```
