# Fit genotype-specific epidemic trajectories using GAM

Fits a generalized additive model (GAM) to replicated disease progress
curves, returning adjusted treatment mean curves evaluated on a common
time grid.

## Usage

``` r
functional_curves(
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
  family_try = c("betar", "quasibinomial"),
  show_progress = TRUE,
  covariates = NULL,
  include_covariates = FALSE,
  covariate_smooths = FALSE,
  ...
)
```

## Arguments

- data:

  A data.frame containing disease progress observations.

- time:

  Character string naming the time variable.

- response:

  Character string naming the response variable.

- treatment:

  Character string naming the treatment/cultivar factor.

- environment:

  Optional character string naming an environment factor.

- block:

  Optional character string naming a blocking factor.

- unit:

  Optional character string naming a unique experimental unit
  identifier.

- response_scale:

  Character string specifying the response scale: `"proportion"` or
  `"percent"`.

- eps:

  Numeric small constant retained for API compatibility.

- min_points:

  Minimum number of observations required per curve.

- grid_n:

  Number of time points for evaluating predicted curves.

- env_ref:

  Reference environment level used for adjusted predictions.

- k_smooth:

  Basis dimension for the global smooth of time.

- k_env:

  Basis dimension for the environment-specific smooth.

- k_trt:

  Basis dimension for the treatment-specific smooth.

- gamma:

  Penalization parameter passed to
  [`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html).

- discrete:

  Logical; whether to use discrete (approximate) fitting.

- family_try:

  Character string specifying the GAM family to try.

- show_progress:

  Logical; whether to show progress.

- covariates:

  Optional character vector of genotype-level covariates (e.g.,
  `c("heading_group")`).

- include_covariates:

  Logical; whether to include covariates as fixed effects in the model.

- covariate_smooths:

  Logical; if TRUE and `include_covariates = TRUE`, adds smooth
  interactions `s(time, by=covariate)` for factor covariates.

- ...:

  Additional arguments.

## Value

An object of class `"functional_curves"` containing:

- `gam`: fitted GAM object and `family_used`;

- `curves`: environment-adjusted treatment mean curves on the common
  grid;

- `grid`: the time grid used;

- `observed_data`: the processed input data;

- `genotype_info`: a tibble with one row per genotype containing
  covariates;

- `vars`: variable names used;

- `settings`: model settings;

- `warnings_betar`: any warnings caught during beta fitting;

- `plot_mean`: ggplot object of the mean curves.

## Details

Genotype-level covariates are descriptors that do not vary within a
genotype, such as phenological groups. For example, in wheat blast
studies, a cultivar's `heading_group` (early, intermediate, late) can be
supplied. This helps distinguish between true genetic resistance and
phenological escape, as cultivars with different heading dates may
encounter different infection-risk windows in the same environment. When
`include_covariates = TRUE`, the model accounts for these covariates,
yielding heading-adjusted functional resistance.

## Examples

``` r
if (FALSE) { # \dontrun{
fc <- functional_curves(
  data = my_data,
  time = "time_var",
  response = "severity",
  treatment = "cultivar",
  environment = "env",
  covariates = c("heading_group"),
  include_covariates = TRUE,
  covariate_smooths = TRUE
)
plot(fc)
} # }
```
