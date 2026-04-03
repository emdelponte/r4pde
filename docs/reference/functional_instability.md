# Normalized functional instability from compare_curves output

Computes normalized functional instability (NFI) for each treatment
based on genotype-by-environment predicted curves extracted from a
[`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md)
object. Optionally, instability can be decomposed into spatial and
temporal components if the environment identifier can be split into
location and year.

## Usage

``` r
functional_instability(
  x,
  n_time = 200,
  env_sep = NULL,
  env_names = c("location", "year"),
  return_curves = FALSE
)
```

## Arguments

- x:

  An object returned by
  [`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md).

- n_time:

  Number of points in the prediction grid over the time domain.

- env_sep:

  Optional separator used to split `env` into spatial and temporal
  components, for example `"_"` or `"-"`. If `NULL`, only the overall
  NFI is returned.

- env_names:

  A character vector of length 2 giving the names of the components
  obtained after splitting `env`. Default is `c("location", "year")`.

- return_curves:

  Logical. If `TRUE`, the predicted genotype-by-environment curves are
  also returned.

## Value

If `return_curves = FALSE`, a tibble with one row per genotype and
columns:

- geno:

  Treatment or genotype identifier.

- n_env:

  Number of environments used in the calculation.

- FI:

  Absolute functional instability.

- mean_energy:

  Integrated squared mean epidemic curve.

- nFI:

  Normalized functional instability. Lower values indicate greater
  stability.

If `env_sep` is provided, the returned tibble also includes:

- FI_space, mean_energy_space, nFI_space:

  Spatial component of instability.

- FI_time, mean_energy_time, nFI_time:

  Temporal component of instability.

If `return_curves = TRUE`, a list is returned with two elements:

- metrics:

  The tibble described above.

- curves:

  A tibble of predicted genotype-by-environment curves with columns
  `geno`, `env`, `time`, `mu`, and `eta`.

## Details

The overall instability metric is defined as the mean integrated squared
deviation of each genotype-by-environment curve from the
genotype-specific mean curve across environments, normalized by the
integrated squared mean curve.

Let \\f\_{ge}(t)\\ denote the predicted epidemic trajectory of genotype
\\g\\ in environment \\e\\, and let \\\bar f_g(t)\\ denote the mean
trajectory of genotype \\g\\ across environments. Functional instability
is computed as:

\$\$ FI_g = \frac{1}{E_g}\sum\_{e=1}^{E_g} \int_T \left(f\_{ge}(t)-\bar
f_g(t)\right)^2 dt \$\$

and normalized as:

\$\$ nFI_g = \frac{FI_g}{\int_T \bar f_g(t)^2 dt} \$\$

Numerical integration is performed with the trapezoidal rule on a
regular prediction grid over the observed time domain.

## See also

[`compare_curves()`](https://emdelponte.github.io/r4pde/reference/compare_curves.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m1 <- r4pde::compare_curves(
  data = dat_ready,
  time = "time",
  response = "y",
  treatment = "geno",
  environment = "env",
  cluster_k = 4
)

# Overall instability
res <- functional_instability(m1)
res

# Overall + spatial and temporal components
# Assuming env has values such as "PF_2021"
res2 <- functional_instability(m1, env_sep = "_")
res2

# Return curves used in the calculations
out <- functional_instability(m1, env_sep = "_", return_curves = TRUE)
out$metrics
head(out$curves)
} # }
```
