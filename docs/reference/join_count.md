# Test for Spatial Join Count Statistics

The function `join_count` calculates spatial join count statistics for a
binary matrix, identifying patterns of aggregation or randomness.

## Usage

``` r
join_count(matrix_data, verbose = TRUE)
```

## Arguments

- matrix_data:

  A binary matrix (with elements 0 and 1) representing the spatial
  distribution of two types of points: 0 for healthy plants (H) and 1
  for diseased plants (D). This matrix reflects the geographical
  distribution or layout of plants in the studied area.

- verbose:

  Logical. If TRUE (default), prints a formatted message to the console.

## Value

A comprehensive, rich-text formatted string of results that includes:

- Statistical counts of specific binary sequences (e.g., "01 or 10",
  "11")

- Expected counts under the assumption of Complete Spatial Randomness
  (CSR)

- Standard deviations and Z-scores (ZHD for "01 or 10" sequences, ZDD
  for "11" sequences)

- Interpretation of whether the spatial distribution for each sequence
  type is "Aggregated" or "Not Aggregated" based on Z-scores

- A summary explaining the implications of these statistics and patterns

The return value aims to provide a clear understanding of the spatial
arrangement's characteristics, aiding in further spatial analysis or
research.

## Details

The function conducts an analysis by first counting the occurrence of
specific sequences ("01 or 10" and "11" - equivalent to HD and DD) in
the binary matrix. It then calculates expected values, standard
deviations, and Z-scores to determine the spatial randomness or
aggregation. The analysis considers both horizontal and vertical
adjacency (rook case) in the matrix.

## References

Madden, L. V., Hughes, G., & van den Bosch, F. (2007). The Study of
Plant Disease Epidemics. The American Phytopathological Society.

## See also

Other Spatial analysis:
[`AFSD()`](https://emdelponte.github.io/r4pde/reference/AFSD.md),
[`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md),
[`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md),
[`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)
