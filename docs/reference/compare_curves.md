# Compare epidemic curves using functional GAM smoothing and clustering

Fits a generalized additive model (GAM) to disease progress curves using
a beta regression (with automatic fallback to quasibinomial), derives
environment-adjusted mean curves for each treatment, computes functional
distances among treatments, and performs hierarchical clustering of
epidemic trajectories.

The function is designed for replicated plant disease epidemics measured
repeatedly over time, optionally across environments and blocks.

## Usage

``` r
compare_curves(
  data,
  time,
  response,
  treatment,
  environment = NULL,
  block = NULL,
  unit = NULL,
  response_scale = c("proportion", "percent"),
  eps = 1e-04,
  min_points = 5,
  grid_n = 140,
  env_ref = NULL,
  k_smooth = 10,
  k_env = 4,
  k_trt = 4,
  gamma = 1.4,
  discrete = TRUE,
  cluster_k = 4,
  hc_method = "ward.D2",
  ...
)
```

## Arguments

- data:

  A data.frame containing epidemic curve data.

- time:

  Character string giving the name of the time variable.

- response:

  Character string giving the name of the response variable (disease
  severity or incidence).

- treatment:

  Character string giving the name of the treatment or cultivar factor.

- environment:

  Optional character string giving the environment factor.

- block:

  Optional character string giving the block or plot factor.

- unit:

  Optional character string giving a unique experimental unit
  identifier. If NULL, a unit is constructed from the interaction of
  environment, treatment, and block (when available).

- response_scale:

  Character string specifying the response scale: either `"proportion"`
  (0–1) or `"percent"` (0–100).

- eps:

  Numeric small constant retained for API compatibility.

- min_points:

  Minimum number of observations required per epidemic curve.

- grid_n:

  Number of time points used to evaluate predicted mean curves.

- env_ref:

  Reference environment level for adjusted predictions. If NULL, the
  first level is used.

- k_smooth:

  Basis dimension for the global smooth over time.

- k_env:

  Basis dimension for the environment-specific smooth.

- k_trt:

  Basis dimension for the treatment-specific smooth.

- gamma:

  Penalization parameter passed to
  [`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html).

- discrete:

  Logical; whether to use discrete fitting in `bam()`.

- cluster_k:

  Number of clusters used to cut the hierarchical tree.

- hc_method:

  Clustering method passed to
  [`hclust()`](https://rdrr.io/r/stats/hclust.html).

- ...:

  Reserved for future extensions.

## Value

An object of class `"r4pde_compare_curves"` containing:

- fitted GAM object,

- environment-adjusted predicted curves,

- functional distance matrix,

- hierarchical clustering,

- treatment cluster assignments.

## Details

The response is internally transformed using the Smithson–Verkuilen
adjustment to avoid exact 0 and 1 values. Functional distances are
computed using an L2 norm over predicted mean curves evaluated on a
common time grid.

## See also

[`plot_curves`](https://emdelponte.github.io/r4pde/reference/plot_curves.md),
[`plot_dendrogram`](https://emdelponte.github.io/r4pde/reference/plot_dendrogram.md),
[`diagnose_curves`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md)
