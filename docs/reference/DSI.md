# Calculate the Disease Severity Index (DSI) (class for each unit)

This function calculates the Disease Severity Index (DSI) based on the
provided unit, class, and maximum class value. The DSI is computed by
aggregating the classes, calculating weights by multiplying the
frequency of each class by the class itself, and then dividing the sum
of these weights by the product of the total number of entries and the
maximum class value, then multiplying by 100.

## Usage

``` r
DSI(unit, class, max)
```

## Arguments

- unit:

  A vector representing the units.

- class:

  A vector representing the classes corresponding to the units.

- max:

  A numeric value representing the maximum possible class value.

## Value

Returns a single numeric value representing the DSI.

## See also

Other Disease quantification:
[`CompMuCens()`](https://emdelponte.github.io/r4pde/reference/CompMuCens.md),
[`DSI2()`](https://emdelponte.github.io/r4pde/reference/DSI2.md)

## Examples

``` r
# Example usage:
unit <- c(1, 2, 3, 4, 5, 6)
class <- c(1, 2, 1, 2, 3, 1)
max <- 3
DSI(unit, class, max)
#> [1] 55.55556
```
