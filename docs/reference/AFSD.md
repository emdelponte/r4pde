# Analysis of foci structure and dynamics (AFSD)

This function performs the analysis of a simple method introduced by
Nelson (1996) and expanded by Laranjeira et al. (1998). The function
assumes the dataframe supplied as input has columns 'x', 'y', and 'i',
where 'x' and 'y' are spatial coordinates and 'i' is a disease indicator
variable (1 if diseased, otherwise 0). The function performs several
steps including filtering rows where 'i' is 1, converting to an
adjacency matrix, and creating foci using igraph. It then calculates
various statistics about the foci and returns these in a list.

## Usage

``` r
AFSD(df)
```

## Arguments

- df:

  A dataframe containing at least three columns: 'x', 'y', and 'i'. 'x'
  and 'y' represent spatial coordinates and 'i' is a disease indicator
  (1 if diseased, otherwise 0).

## Value

A list containing: cluster_summary2: a dataframe summarizing the number
and size of foci, and proportions of diseased plants. cluster_df: a
dataframe containing foci information, including size and number of rows
and columns in each foci. df_clustered: the original dataframe with an
added 'focus_id' column, showing which foci each row belongs to.

## See also

Other Spatial analysis:
[`BPL()`](https://emdelponte.github.io/r4pde/reference/BPL.md),
[`count_subareas()`](https://emdelponte.github.io/r4pde/reference/count_subareas.md),
[`count_subareas_random()`](https://emdelponte.github.io/r4pde/reference/count_subareas_random.md),
[`fit_gradients()`](https://emdelponte.github.io/r4pde/reference/fit_gradients.md),
[`join_count()`](https://emdelponte.github.io/r4pde/reference/join_count.md),
[`oruns_test()`](https://emdelponte.github.io/r4pde/reference/oruns_test.md),
[`oruns_test_boustrophedon()`](https://emdelponte.github.io/r4pde/reference/oruns_test_boustrophedon.md),
[`oruns_test_byrowcol()`](https://emdelponte.github.io/r4pde/reference/oruns_test_byrowcol.md),
[`plot_AFSD()`](https://emdelponte.github.io/r4pde/reference/plot_AFSD.md)

## Examples

``` r
# Generate a sample dataframe
set.seed(123)
df <- data.frame(x = sample(1:100, 500, replace = TRUE),
                 y = sample(1:100, 500, replace = TRUE),
                 i = sample(0:1, 500, replace = TRUE, prob = c(0.7, 0.3)))

# Perform the AFSD
result <- AFSD(df)
```
