# Compare epidemic curves using GAM smoothing, functional distances, clustering, and optional permutation/bootstrapping

Fits a generalized additive model (GAM) to replicated disease progress
curves. By default, the response is modeled with beta regression (logit
link) after a Smithson–Verkuilen adjustment to handle 0 and 1 values; if
beta precision (theta) estimation fails, the model falls back to a
quasibinomial GAM fit on the unclipped proportion response.

From the fitted model, the function derives (i) environment-adjusted
treatment mean curves on a common time grid, (ii) L2 functional
distances among treatment mean trajectories followed by hierarchical
clustering, and (iii) curve-level fitted trajectories (excluding unit
random effects) to build a curve-by-curve distance matrix. Optionally,
it performs a restricted permutation test on the curve-level distance
matrix for a priori group labels, and computes bootstrap envelopes and
distance uncertainty summaries from resampled predicted curves.

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
  k_trt = 6,
  gamma = 1.4,
  discrete = TRUE,
  cluster_k = 4,
  hc_method = "ward.D2",
  test_factor = NULL,
  n_perm = 999,
  test_mode = c("auto", "none", "global", "global_pairwise"),
  min_strata_for_pairwise = 6,
  perm_unit = NULL,
  perm_strata = NULL,
  bootstrap = FALSE,
  boot_B = 399,
  boot_seed = 1,
  boot_ci = c(0.025, 0.975),
  bootstrap_mode = c("predicted", "refit"),
  ...
)
```

## Arguments

- data:

  A data.frame containing disease progress observations in long format.

- time:

  Character string naming the time variable (numeric or coercible to
  numeric).

- response:

  Character string naming the response variable (disease
  incidence/severity).

- treatment:

  Character string naming the treatment/cultivar factor.

- environment:

  Optional character string naming an environment factor (e.g.,
  site-year).

- block:

  Optional character string naming a blocking factor.

- unit:

  Optional character string naming a unique experimental unit (curve)
  identifier. If `NULL`, a unit is constructed from the interaction of
  available factors (environment, treatment, and block).

- response_scale:

  Character string specifying the response scale: `"proportion"` (0–1)
  or `"percent"` (0–100).

- eps:

  Numeric small constant retained for API compatibility (not used for
  beta handling).

- min_points:

  Minimum number of observations required per curve (unit) to be
  retained.

- grid_n:

  Number of time points for evaluating predicted curves on a common
  grid.

- env_ref:

  Reference environment level used for adjusted treatment mean
  predictions. If `NULL`, the first level of `environment` is used.

- k_smooth:

  Basis dimension for the global smooth of time.

- k_env:

  Basis dimension for the environment-specific smooth (if `environment`
  is provided).

- k_trt:

  Basis dimension for the treatment-specific smooth.

- gamma:

  Penalization parameter passed to
  [`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html).

- discrete:

  Logical; whether to use discrete (approximate) fitting in
  [`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html).

- cluster_k:

  Number of clusters used to cut the hierarchical tree for treatments.

- hc_method:

  Clustering linkage method passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

- test_factor:

  Optional a priori grouping used for a distance-based permutation test.
  Either (i) a single character naming a column in `data` that maps
  uniquely to the permutation unit (see `perm_unit`), or (ii) a named
  vector whose names are permutation-unit identifiers and whose values
  are group labels.

- n_perm:

  Number of permutations for the distance-based test.

- test_mode:

  Character; permutation test mode. One of `"auto"`, `"none"`,
  `"global"`, or `"global_pairwise"`. In `"auto"`, pairwise tests may be
  disabled when restricted strata support is limited.

- min_strata_for_pairwise:

  Minimum number of strata levels (e.g., blocks) required to enable
  pairwise tests under `test_mode = "auto"`.

- perm_unit:

  Optional character naming the unit at which group labels are permuted
  (e.g., genotype/cultivar). If `NULL`, permutation is performed at the
  curve (unit) level.

- perm_strata:

  Optional character naming a factor defining restricted permutation
  strata (e.g., block or environment). If `NULL`, defaults to `block`
  when provided.

- bootstrap:

  Logical; if `TRUE`, compute bootstrap envelopes and distance summaries
  from resampled predicted curves.

- boot_B:

  Integer number of bootstrap replicates.

- boot_seed:

  Integer seed for bootstrap resampling.

- boot_ci:

  Length-2 numeric vector of quantiles for bootstrap intervals (e.g.,
  `c(0.025, 0.975)`).

- bootstrap_mode:

  Character; bootstrap strategy. Currently supports `"predicted"`
  (resampling predicted curves). `"refit"` is not implemented in this
  version.

- ...:

  Reserved for future extensions.

## Value

An object of class `"r4pde_compare_curves"` containing:

- `gam`: fitted GAM object and `family_used`;

- `pred`: environment-adjusted treatment mean curves on the common grid;

- `distance`: treatment-by-treatment functional distance matrix and
  `hc`;

- `clusters`: treatment cluster assignments and summary scores;

- `curve_distance`: curve-by-curve distance matrix;

- `test`: optional distance-based permutation test results;

- `bootstrap`: optional bootstrap envelopes and distance summaries;

- diagnostic fields including `settings` and captured beta-fit warnings.

## Details

Functional distances among treatments are computed as an L2 norm over
the predicted mean curves evaluated on a common time grid. Curve-level
distances are computed analogously using fitted trajectories (excluding
unit random effects), and are scaled by the median non-zero distance for
numerical stability; interpret these distances as relative rather than
absolute units.

The permutation test uses a PERMANOVA-like pseudo-\\F\\ statistic
computed from the curve-level distance matrix, with optional restricted
permutations within strata. When enabled, pairwise comparisons are
adjusted using Holm's method.

## See also

[`plot_curves`](https://emdelponte.github.io/r4pde/reference/plot_curves.md),
[`plot_dendrogram`](https://emdelponte.github.io/r4pde/reference/plot_dendrogram.md),
[`diagnose_curves`](https://emdelponte.github.io/r4pde/reference/diagnose_curves.md)
