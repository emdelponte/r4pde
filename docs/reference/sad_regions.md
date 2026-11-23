# Severity ranges where conditions differ

Summarises the differences of a
[`sad_curve()`](https://emdelponte.github.io/r4pde/reference/sad_curve.md)
object into the severity ranges where the simultaneous band excludes
zero.

## Usage

``` r
sad_regions(curve)
```

## Arguments

- curve:

  An object returned by
  [`sad_curve()`](https://emdelponte.github.io/r4pde/reference/sad_curve.md).

## Value

A tibble with `contrast`, `from`, `to` (actual severity, %) and
`direction` of the difference in bias (`"a > b"` or `"a < b"`).
