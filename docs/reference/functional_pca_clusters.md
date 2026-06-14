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

  Method for suggesting k: "none", "silhouette", or "elbow".

- ...:

  Additional arguments passed to clustering functions.
