# Functional suppression profiles of disease control treatments

Computes disease suppression trajectories relative to a reference
treatment, estimates functional distances among suppression curves,
clusters treatments according to their temporal suppression profiles,
and summarizes each profile using functional suppression metrics.

## Usage

``` r
functional_suppression_profiles(
  data,
  reference,
  time = "time",
  response = "severity",
  treatment = "treatment",
  environment = NULL,
  by_environment = FALSE,
  env_ref = NULL,
  threshold = 5,
  metrics = c("protected_area", "max_suppression", "persistence", "centroid"),
  dist_method = "euclidean",
  hclust_method = "average",
  k = "auto",
  cut_height = NULL,
  ...
)
```

## Arguments

- data:

  A data frame or tibble containing the disease progress data.

- reference:

  Character string specifying the reference treatment name.

- time:

  Character string specifying the time column. Default is `"time"`.

- response:

  Character string specifying the response column. Default is
  `"severity"`.

- treatment:

  Character string specifying the treatment column. Default is
  `"treatment"`.

- environment:

  Character string specifying the environment column (optional).

- by_environment:

  Logical; if `TRUE`, the analysis is run separately for each
  environment. Default is `FALSE`.

- env_ref:

  Character string specifying the reference environment for joint
  analysis (optional).

- threshold:

  Numeric threshold for the persistence metric. Default is `5`.

- metrics:

  Character vector of metrics to include in the ranking. Default is
  `c("protected_area", "max_suppression", "persistence", "centroid")`.

- dist_method:

  Character string specifying the distance method. Default is
  `"euclidean"`.

- hclust_method:

  Character string specifying the hierarchical clustering method.
  Default is `"average"`.

- k:

  Integer specifying the number of functional profiles to create, or
  `"auto"` (default). When `"auto"`, the optimal number of profiles is
  estimated by maximizing the average silhouette width.

- cut_height:

  Numeric specifying the cut height for clustering. If provided,
  overrides `k`.

- ...:

  Additional arguments passed to
  [`functional_contrast`](https://emdelponte.github.io/r4pde/reference/functional_contrast.md).

## Value

A list of class `"functional_suppression_profiles"` containing:

- `call`: The matched call.

- `reference`: The reference treatment name.

- `contrast`: The contrast data from `functional_contrast`.

- `dsp_curves`: The fitted curves from `functional_curves`.

- `distances`: The distances object from `functional_distances`.

- `hclust`: The hierarchical clustering object.

- `clusters`: A tibble with treatment cluster assignments (internal).

- `profiles`: A lightweight tibble mapping treatments to profiles.

- `classification`: A tibble with the final classification combining
  profiles, mean rank, and functional metrics, ordered by mean rank.

- `summary`: A tibble with functional summary metrics.

- `ranking`: A tibble with metric rankings.

- `cluster_summary`: A tibble combining metrics, ranks, and cluster
  assignments.

- `parameters`: A list of input parameters.

## Details

A Functional Suppression Profile (FSP) is the temporal trajectory of
disease suppression produced by a treatment relative to an untreated or
reference control. Functional distances and hierarchical clustering are
computed from smoothed disease suppression trajectories rather than raw
observed suppression values. Functional summary metrics are subsequently
used to interpret the resulting suppression profiles.

The workflow includes contrasting treatments against a reference,
fitting functional curves, calculating distances, clustering, and
summarizing the suppression using metrics like protected area, maximum
suppression, and persistence.

Treatments belonging to the same functional suppression profile are
displayed using the same colour throughout all graphical summaries. This
visual consistency facilitates interpretation of the relationship
between suppression trajectories, functional clustering, and suppression
metrics.

The `classification` table combines functional profile membership,
suppression metrics, and mean functional rank. Functional profiles are
defined from distances among smoothed suppression trajectories, whereas
mean rank summarizes the overall performance across functional
suppression metrics.

## See also

[`functional_contrast`](https://emdelponte.github.io/r4pde/reference/functional_contrast.md),
[`functional_curves`](https://emdelponte.github.io/r4pde/reference/functional_curves.md),
[`functional_distances`](https://emdelponte.github.io/r4pde/reference/functional_distances.md),
[`functional_summary`](https://emdelponte.github.io/r4pde/reference/functional_summary.md),
[`rank_dsp`](https://emdelponte.github.io/r4pde/reference/rank_dsp.md),
[`plot_dendrogram`](https://emdelponte.github.io/r4pde/reference/plot_dendrogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim_dat <- tibble::tibble(
  treatment = rep(c("Control", "A", "B", "C"), each = 6),
  time = rep(seq(0, 25, by = 5), times = 4),
  severity = c(
    c(5, 10, 20, 35, 50, 65),
    c(3, 5, 10, 18, 30, 40),
    c(4, 6, 12, 22, 35, 45),
    c(2, 4, 8, 14, 22, 32)
  )
)
fsp <- functional_suppression_profiles(
  data = sim_dat,
  reference = "Control",
  time = "time",
  response = "severity",
  treatment = "treatment",
  threshold = 5,
  k = 2
)

fsp
summary(fsp)
plot(fsp, type = "dendrogram")
plot(fsp, type = "profiles")
plot(fsp, type = "heatmap")
plot(fsp, type = "rank")
} # }
```
