# Suggest GAM smoothing parameters for epidemic curve models

A helper function to guide selection of the `k_smooth`, `k_trt`,
`k_env`, and `gamma` parameters used by
[`functional_curves`](https://emdelponte.github.io/r4pde/reference/functional_curves.md).
Can infer the effective replication from a data frame or accept scalar
values directly when the data are not yet available.

For sparse epidemic data it is strongly recommended to use
`rule = "minimum"` so that k values are based on the least-replicated
treatment-by-environment combination, guarding against over-fitting.

## Usage

``` r
suggest_k(
  data = NULL,
  time = NULL,
  treatment = NULL,
  environment = NULL,
  n_time = NULL,
  n_env = NULL,
  rule = c("minimum", "median"),
  smoothness = c("conservative", "moderate", "flexible")
)
```

## Arguments

- data:

  Optional data frame. If supplied, `time`, `treatment`, and optionally
  `environment` are used to compute summaries of the number of unique
  time points and environments per treatment.

- time:

  Unquoted column name for the time variable, or a character string
  naming the column.

- treatment:

  Unquoted column name for the treatment / cultivar variable, or a
  character string naming the column.

- environment:

  Optional unquoted column name for the environment variable, or a
  character string naming the column.

- n_time:

  Integer. Number of unique time points to use directly (ignored when
  `data` is supplied).

- n_env:

  Integer. Number of unique environments to use directly (ignored when
  `data` is supplied, or when there is no environment variable).

- rule:

  Character string. How to summarise the distribution of replication
  counts across treatment-by-environment combinations: `"minimum"`
  (default, conservative) or `"median"`.

- smoothness:

  Character string controlling how liberal the recommendations are:
  `"conservative"` (default), `"moderate"`, or `"flexible"`.

## Value

A named list with the following elements:

- `time_summary`:

  Named numeric vector with minimum, median, and maximum unique time
  points per treatment-by-environment combination.

- `environment_summary`:

  Named numeric vector with minimum, median, and maximum unique
  environments per treatment (or `NULL` when no environment variable is
  given).

- `effective_n_time`:

  The effective number of time points chosen by `rule`.

- `effective_n_env`:

  The effective number of environments chosen by `rule`, or `NULL`.

- `k_smooth`:

  Recommended basis dimension for the global time smooth.

- `k_trt`:

  Recommended basis dimension for the treatment-specific smooth.

- `k_env`:

  Recommended basis dimension for the environment random effect.

- `gamma`:

  Recommended penalisation multiplier.

- `message`:

  A short interpretation message.

## Details

The function computes, for every treatment-by-environment combination,
the number of unique time points at which observations are available. It
then summarises these counts using the chosen `rule` to obtain
`effective_n_time`. Similarly it computes, per treatment, the number of
unique environments, and summarises using the same `rule` to obtain
`effective_n_env`.

Recommended k values follow these heuristics:

- `k_smooth`: global smooth over time; capped below
  `effective_n_time - 1`.

- `k_trt`: treatment-specific smooth; capped below `k_smooth`.

- `k_env`: environment random effect; capped below `effective_n_env`.

- `gamma`: penalty multiplier; higher values encourage smoother fits.

## Examples

``` r
# Using explicit values
suggest_k(n_time = 5, n_env = 8)
#> Error in suggest_k(n_time = 5, n_env = 8): could not find function "suggest_k"

# Inferring from a data frame (unquoted column names)
df <- data.frame(
  time      = rep(1:6, times = 6),
  cultivar  = rep(c("A", "B", "C"), each = 12),
  env       = rep(rep(c("E1", "E2"), each = 6), times = 3),
  severity  = runif(36, 0, 0.5)
)
suggest_k(
  data        = df,
  time        = time,
  treatment   = cultivar,
  environment = env,
  rule        = "minimum",
  smoothness  = "conservative"
)
#> Error in suggest_k(data = df, time = time, treatment = cultivar, environment = env,     rule = "minimum", smoothness = "conservative"): could not find function "suggest_k"
```
