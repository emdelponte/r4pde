# Intrinsic difficulty of each image

Extracts the image random effects of a
[`sad_fit()`](https://emdelponte.github.io/r4pde/reference/sad_fit.md)
model: the shift, common to all raters and conditions, of the estimates
of each image relative to what its actual severity predicts. Large
absolute values flag images that are hard to assess, or whose actual
severity may be wrong.

## Usage

``` r
sad_images(fit)
```

## Arguments

- fit:

  An object returned by
  [`sad_fit()`](https://emdelponte.github.io/r4pde/reference/sad_fit.md).

## Value

A tibble with `image`, `actual`, `difficulty` (percentage points for the
gaussian family) and its conditional standard error `se`, ordered by
absolute difficulty.
