# Functional principal component analysis of disease progress curves

Performs functional principal component analysis on fitted disease
progress curves returned by
[`functional_curves`](https://emdelponte.github.io/r4pde/reference/functional_curves.md).
The function decomposes variation among epidemic trajectories into
orthogonal temporal components and returns curve-level scores,
eigenfunctions, variance explained, and reconstructed curves.

## Usage

``` r
functional_pca(object, ...)

# S3 method for class 'functional_curves'
functional_pca(
  object,
  n_components = NULL,
  var_explained = 0.95,
  center = TRUE,
  scale = FALSE,
  method = c("pca_on_grid"),
  ...
)

# S3 method for class 'functional_dsp'
functional_pca(
  object,
  n_components = NULL,
  var_explained = 0.95,
  center = TRUE,
  scale = FALSE,
  method = c("pca_on_grid"),
  ...
)
```

## Arguments

- object:

  An object returned by
  [`functional_curves`](https://emdelponte.github.io/r4pde/reference/functional_curves.md).

- ...:

  Additional arguments for future extensions.

- n_components:

  Optional integer number of functional principal components to retain.

- var_explained:

  Cumulative variance threshold used when `n_components = NULL`.

- center:

  Logical; whether to center curves before PCA. Default TRUE.

- scale:

  Logical; whether to scale grid columns before PCA. Default FALSE.

- method:

  Character; method for FPCA, currently only `"pca_on_grid"` is
  supported.

## Value

An object of class `"r4pde_functional_pca"` containing:

- `scores`: Tibble of curve-level FPC scores.

- `eigenfunctions`: Tibble of eigenfunction values across the time grid.

- `variance`: Tibble of eigenvalues, variance explained, and cumulative
  variance.

- `mean_curve`: Tibble of the mean curve.

- `reconstructed`: Tibble of reconstructed curves using retained
  components.

- `input_curves`: Tibble of original fitted curves.

- `pca`: The underlying `prcomp` object.

- `settings`: List of settings used.

- `call`: The matched call.

## Details

The function uses the fitted curves from
[`functional_curves()`](https://emdelponte.github.io/r4pde/reference/functional_curves.md)
and does not refit the disease progress model. The first implementation
uses PCA on a common prediction grid. FPC scores can be analyzed as
functional epidemiological traits in downstream models.

Interpretation:

- FPC1 often captures the largest mode of variation, commonly overall
  epidemic intensity or speed.

- Later FPCs may capture timing of disease onset, curve crossing, late
  acceleration, or other shape-related deviations.

- Interpretation must be data-driven and should be based on
  eigenfunction plots and mean +/- perturbation plots.

## Examples

``` r
if (FALSE) { # \dontrun{
curves <- functional_curves(...)
fpca <- functional_pca(curves, var_explained = 0.95)

print(fpca)
plot(fpca, type = "scree")
plot(fpca, type = "components")
plot(fpca, type = "scores", components = c(1, 2))

scores <- get_fpca_scores(fpca)
} # }
```
