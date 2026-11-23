# SADRustSoybean dataset

Severity of soybean rust estimated by 40 raters on 50 leaves, first
without and then with a standard area diagram (SAD). Raters 1 to 20 used
the old SAD and raters 21 to 40 the new SAD proposed by Franceschi et
al. (2020). Use it with
[`sad_data()`](https://emdelponte.github.io/r4pde/reference/sad_data.md)
and the other `sad_*` functions.

## Usage

``` r
SADRustSoybean
```

## Format

A data frame with 2000 rows and 6 columns:

- sad:

  SAD used by the rater in the aided assessment: `"Old"` or `"New"`

- rater:

  Rater identifier (1 to 40)

- leaf:

  Leaf identifier (1 to 50)

- actual:

  Actual severity (%) measured by image analysis

- unaided:

  Severity (%) estimated without the SAD

- aided:

  Severity (%) estimated with the SAD

## Source

Franceschi, V. T., Alves, K. S., Mazaro, S. M., Godoy, C. V., Duarte, H.
S. S., and Del Ponte, E. M. (2020). A new standard area diagram set for
assessment of severity of soybean rust improves accuracy of estimates
and optimizes resource use. Plant Pathology 69:495-505.
