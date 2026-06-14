# Plot functional PCA results

Plot functional PCA results

## Usage

``` r
# S3 method for class 'r4pde_functional_pca'
plot(
  x,
  type = c("scree", "components", "scores", "reconstruction", "mean"),
  components = c(1, 2),
  curve_id = NULL,
  ...
)

# S3 method for class 'r4pde_functional_pca'
autoplot(object, ...)
```

## Arguments

- x:

  An object of class `"r4pde_functional_pca"`.

- type:

  Type of plot: "scree", "components", "scores", "reconstruction", or
  "mean".

- components:

  Integer vector of length 2 indicating which components to plot for
  scores and mean perturbation.

- curve_id:

  Optional vector of curve IDs to include in reconstruction plot.

- ...:

  Additional arguments.

- object:

  An object of class `"r4pde_functional_pca"`.
