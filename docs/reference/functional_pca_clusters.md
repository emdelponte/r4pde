# Cluster curves based on FPCA scores

Cluster curves based on FPCA scores

## Usage

``` r
functional_pca_clusters(
  x,
  k = NULL,
  components = NULL,
  method = c("kmeans", "hclust"),
  choose_k = c("none", "silhouette", "elbow"),
  ...
)
```

## Arguments

- x:

  An object of class `"r4pde_functional_pca"`.

- k:

  Number of clusters.

- components:

  Integer vector of components to use.

- method:

  Clustering method: "kmeans" or "hclust".

- choose_k:

  Method used to choose `k` when `k = NULL`: `"silhouette"` (largest
  mean silhouette width over k = 2 to 10) or `"elbow"` (largest drop
  below the line joining the within-cluster sums of squares of k = 1 and
  k = 10). `"none"` requires `k`. The criterion table is attached to the
  result as attribute `"k_selection"`.

- ...:

  Additional arguments passed to clustering functions.
