# Window Pane for Epidemiological Analysis

This function calculates summary statistics within specified windows
around a given end date in a dataset, facilitating epidemiological
analysis. It allows backward, forward, or both directions of window
calculations based on a user-defined variable and window lengths.

## Usage

``` r
windowpane(
  data,
  end_date_col,
  date_col,
  variable,
  summary_type,
  threshold = NULL,
  window_lengths,
  direction = "backward",
  group_by_cols = NULL,
  date_format = "%Y-%m-%d"
)
```

## Arguments

- data:

  A data frame containing the input data.

- end_date_col:

  A string specifying the name of the column representing the end date.

- date_col:

  A string specifying the name of the column representing the date
  variable.

- variable:

  A string specifying the name of the column for which summary
  statistics are calculated.

- summary_type:

  A string specifying the type of summary to calculate. Options are
  "mean", "sum", "above_threshold", or "below_threshold".

- threshold:

  Optional numeric value used when `summary_type` is "above_threshold"
  or "below_threshold".

- window_lengths:

  A numeric vector specifying the window lengths (in days) for the
  calculations.

- direction:

  A string specifying the direction of the window. Options are
  "backward" (default), "forward", or "both".

- group_by_cols:

  Optional vector of strings specifying column names for grouping the
  data.

- date_format:

  A string specifying the format of the date columns. Default is
  "%Y-%m-%d".

## Value

A data frame with the calculated summary values for each window.

## See also

Other Disease modeling:
[`get_nasapower()`](https://emdelponte.github.io/r4pde/reference/get_nasapower.md)
