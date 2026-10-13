# Compute pairwise functional distances and cluster curves

Computes L2 functional distances among treatment mean trajectories
evaluated on a common time grid from a `functional_curves` object. Also
computes curve-level distances, performs hierarchical clustering, and
optionally performs permutation testing.

## Usage

``` r
functional_distances(
  object,
  cluster_k = 4,
  hc_method = "ward.D2",
  test_factor = NULL,
  n_perm = 999,
  seed = 1,
  test_mode = c("auto", "none", "global", "global_pairwise"),
  min_strata_for_pairwise = 6,
  perm_unit = NULL,
  perm_strata = NULL,
  bootstrap = FALSE,
  boot_B = 399,
  boot_seed = 1,
  boot_ci = c(0.025, 0.975),
  show_progress = TRUE,
  curve_level = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `functional_curves`.

- cluster_k:

  Number of clusters used to cut the hierarchical tree.

- hc_method:

  Clustering linkage method passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

- test_factor:

  Optional a priori grouping used for a distance-based permutation test.

- n_perm:

  Number of permutations for the distance-based test.

- seed:

  Integer seed for the permutation test, or `NULL` to use the current
  random stream. The caller's random stream is restored afterwards.

- test_mode:

  Character; permutation test mode.

- min_strata_for_pairwise:

  Minimum number of strata levels required for pairwise tests.

- perm_unit:

  Optional character naming the unit at which group labels are permuted
  (default: the curve unit). Use it when `test_factor` is a property of
  a coarser unit, such as a genotype-level resistance class: all curves
  of a unit keep the same label in every permutation, so curves of the
  same genotype are not treated as exchangeable.

- perm_strata:

  Optional character naming a factor defining restricted permutation
  strata. When `perm_unit` is coarser than the curve, strata are
  honoured only if they are constant within each permutation unit;
  otherwise whole units are permuted without strata, with a warning.

- bootstrap:

  Logical; if `TRUE`, compute bootstrap envelopes and distance
  summaries.

- boot_B:

  Integer number of bootstrap replicates.

- boot_seed:

  Integer seed for bootstrap resampling. The caller's random stream is
  restored afterwards.

- boot_ci:

  Length-2 numeric vector of quantiles for bootstrap intervals.

- show_progress:

  Logical; whether to show progress.

- curve_level:

  Logical; whether to compute curve-level distance matrix.

- ...:

  Additional arguments.

## Value

An object of class `"functional_distances"` containing:

- `distance`: treatment-by-treatment functional distance matrix;

- `distance_table`: pairwise distances as a data frame;

- `hc`: `hclust` object;

- `clusters`: treatment cluster assignments;

- `curve_distance`: curve-by-curve distance matrix;

- `test`: permutation test results;

- `bootstrap`: bootstrap results if requested.

## Examples

``` r
if (FALSE) { # \dontrun{
fd <- functional_distances(fc, cluster_k = 4)
print(fd)
plot_dendrogram(fd)
plot_curves(fd)
} # }
```
