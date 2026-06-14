# Rank Treatments based on Functional DSP Descriptors

Generates rankings of treatments based on extracted DSP functional
descriptors.

## Usage

``` r
rank_dsp(
  summary_data,
  metrics = c("protected_area", "max_suppression", "persistence", "centroid"),
  higher_is_better = TRUE,
  group = NULL,
  ties.method = "average"
)
```

## Arguments

- summary_data:

  The tibble output from
  [`functional_summary`](https://emdelponte.github.io/r4pde/reference/functional_summary.md).

- metrics:

  Character vector of metrics to rank.

- higher_is_better:

  Logical vector or named logical vector indicating if higher values are
  better for each metric. If a single logical, applies to all.

- group:

  Optional character vector of grouping variables.

- ties.method:

  Character string specifying how to handle ties. Default is
  `"average"`.

## Value

A tibble with ranks for each specified metric, and an `average_rank`.
