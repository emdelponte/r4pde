# Plot DSP Rankings Heatmap

Creates a heatmap visualization of treatment rankings based on DSP
descriptors.

## Usage

``` r
plot_dsp_rank_heatmap(ranks_data, treatment = "treatment", group = NULL)
```

## Arguments

- ranks_data:

  Tibble output from
  [`rank_dsp`](https://emdelponte.github.io/r4pde/reference/rank_dsp.md).

- treatment:

  Character string naming the treatment column (for the y-axis).

- group:

  Optional character string for faceting by group.

## Value

A `ggplot` object.
