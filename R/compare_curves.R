#' Compare epidemic curves using GAM smoothing, functional distances, clustering, and optional permutation/bootstrapping
#'
#' @description
#' Fits a generalized additive model (GAM) to replicated disease progress curves.
#' By default, the response is modeled with beta regression (logit link) after a
#' Smithson–Verkuilen adjustment to handle 0 and 1 values; if beta precision
#' (theta) estimation fails, the model falls back to a quasibinomial GAM fit on
#' the unclipped proportion response.
#'
#' From the fitted model, the function derives (i) environment-adjusted treatment
#' mean curves on a common time grid, (ii) L2 functional distances among
#' treatment mean trajectories followed by hierarchical clustering, and (iii)
#' curve-level fitted trajectories (excluding unit random effects) to build a
#' curve-by-curve distance matrix. Optionally, it performs a restricted
#' permutation test on the curve-level distance matrix for a priori group labels,
#' and computes bootstrap envelopes and distance uncertainty summaries from
#' resampled predicted curves.
#'
#' @param data A data.frame containing disease progress observations in long format.
#' @param time Character string naming the time variable (numeric or coercible to numeric).
#' @param response Character string naming the response variable (disease incidence/severity).
#' @param treatment Character string naming the treatment/cultivar factor.
#' @param environment Optional character string naming an environment factor (e.g., site-year).
#' @param block Optional character string naming a blocking factor.
#' @param unit Optional character string naming a unique experimental unit (curve) identifier.
#'   If \code{NULL}, a unit is constructed from the interaction of available factors
#'   (environment, treatment, and block).
#' @param response_scale Character string specifying the response scale:
#'   \code{"proportion"} (0--1) or \code{"percent"} (0--100).
#' @param eps Numeric small constant retained for API compatibility (not used for beta handling).
#' @param min_points Minimum number of observations required per curve (unit) to be retained.
#' @param grid_n Number of time points for evaluating predicted curves on a common grid.
#' @param env_ref Reference environment level used for adjusted treatment mean predictions.
#'   If \code{NULL}, the first level of \code{environment} is used.
#' @param k_smooth Basis dimension for the global smooth of time.
#' @param k_env Basis dimension for the environment-specific smooth (if \code{environment} is provided).
#' @param k_trt Basis dimension for the treatment-specific smooth.
#' @param gamma Penalization parameter passed to \code{mgcv::bam()}.
#' @param discrete Logical; whether to use discrete (approximate) fitting in \code{mgcv::bam()}.
#' @param family_try Character string specifying the GAM family to try:
#'   \code{"betar"} (beta regression) or \code{"quasibinomial"}. If \code{"betar"}
#'   is selected but precision estimation fails, the model falls back to
#'   \code{"quasibinomial"}.
#' @param cluster_k Number of clusters used to cut the hierarchical tree for treatments.
#' @param hc_method Clustering linkage method passed to \code{stats::hclust()}.
#'
#' @param test_factor Optional a priori grouping used for a distance-based permutation test.
#'   Either (i) a single character naming a column in \code{data} that maps uniquely to the
#'   permutation unit (see \code{perm_unit}), or (ii) a named vector whose names are
#'   permutation-unit identifiers and whose values are group labels.
#' @param n_perm Number of permutations for the distance-based test.
#' @param test_mode Character; permutation test mode. One of \code{"auto"}, \code{"none"},
#'   \code{"global"}, or \code{"global_pairwise"}. In \code{"auto"}, pairwise tests may be
#'   disabled when restricted strata support is limited.
#' @param min_strata_for_pairwise Minimum number of strata levels (e.g., blocks) required to
#'   enable pairwise tests under \code{test_mode = "auto"}.
#' @param perm_unit Optional character naming the unit at which group labels are permuted
#'   (e.g., genotype/cultivar). If \code{NULL}, permutation is performed at the curve (unit) level.
#' @param perm_strata Optional character naming a factor defining restricted permutation strata
#'   (e.g., block or environment). If \code{NULL}, defaults to \code{block} when provided.
#'
#' @param bootstrap Logical; if \code{TRUE}, compute bootstrap envelopes and distance summaries
#'   from resampled predicted curves.
#' @param boot_B Integer number of bootstrap replicates.
#' @param boot_seed Integer seed for bootstrap resampling.
#' @param boot_ci Length-2 numeric vector of quantiles for bootstrap intervals (e.g., \code{c(0.025, 0.975)}).
#' @param bootstrap_mode Character; bootstrap strategy. Currently supports \code{"predicted"}
#'   (resampling predicted curves). \code{"refit"} is not implemented in this version.
#' @param show_progress Logical; whether to show progress bars for long-running tasks.
#' @param curve_level Logical; whether to compute curve-level distance matrix. If \code{NULL} (default),
#'   it is automatically set to \code{TRUE} if \code{test_factor} is provided.
#' @param ... Reserved for future extensions.
#'
#' @details
#' Functional distances among treatments are computed as an L2 norm over the
#' predicted mean curves evaluated on a common time grid. Curve-level distances
#' are computed analogously using fitted trajectories (excluding unit random effects),
#' and are scaled by the median non-zero distance for numerical stability; interpret
#' these distances as relative rather than absolute units.
#'
#' The permutation test uses a PERMANOVA-like pseudo-\eqn{F} statistic computed from
#' the curve-level distance matrix, with optional restricted permutations within
#' strata. When enabled, pairwise comparisons are adjusted using Holm's method.
#'
#' @return An object of class \code{"r4pde_compare_curves"} containing:
#' \itemize{
#'   \item \code{gam}: fitted GAM object and \code{family_used};
#'   \item \code{pred}: environment-adjusted treatment mean curves on the common grid;
#'   \item \code{distance}: treatment-by-treatment functional distance matrix and \code{hc};
#'   \item \code{clusters}: treatment cluster assignments and summary scores;
#'   \item \code{curve_distance}: curve-by-curve distance matrix;
#'   \item \code{test}: optional distance-based permutation test results;
#'   \item \code{bootstrap}: optional bootstrap envelopes and distance summaries;
#'   \item diagnostic fields including \code{settings} and captured beta-fit warnings.
#' }
#'
#' @seealso
#' \code{\link{plot_curves}},
#' \code{\link{plot_dendrogram}},
#' \code{\link{diagnose_curves}}
#'
#' @export
compare_curves <- function(
    data,
    time,
    response,
    treatment,
    environment = NULL,
    block = NULL,
    unit = NULL,
    response_scale = c("proportion", "percent"),
    eps = 1e-4,
    min_points = 5,
    grid_n = 140,
    env_ref = NULL,
    k_smooth = 10,
    k_env = 4,
    k_trt = 6,
    gamma = 1.4,
    discrete = TRUE,
    family_try = c("betar","quasibinomial"),
    cluster_k = 4,
    hc_method = "ward.D2",

    # distance-based permutation test (optional)
    test_factor = NULL,
    n_perm = 999,
    test_mode = c("auto","none","global","global_pairwise"),
    min_strata_for_pairwise = 6,

    # crucial for phenotyping
    perm_unit = NULL,

    # restricted permutation strata
    perm_strata = NULL,

    # bootstrap + plots
    bootstrap = FALSE,
    boot_B = 399,
    boot_seed = 1,
    boot_ci = c(0.025, 0.975),
    bootstrap_mode = c("predicted","refit"),
    show_progress = TRUE,
    curve_level = NULL, 
    ...

){
  warning("`compare_curves()` is soft-deprecated. Please use the new workflow with `functional_curves()` and `functional_distances()`.", call. = FALSE)
  
  response_scale <- match.arg(response_scale)
  test_mode <- match.arg(test_mode)
  bootstrap_mode <- match.arg(bootstrap_mode)
  family_try <- match.arg(family_try)
  
  fc <- functional_curves(
    data = data, time = time, response = response, treatment = treatment,
    environment = environment, block = block, unit = unit,
    response_scale = response_scale, eps = eps, min_points = min_points,
    grid_n = grid_n, env_ref = env_ref, k_smooth = k_smooth, k_env = k_env,
    k_trt = k_trt, gamma = gamma, discrete = discrete, family_try = family_try,
    show_progress = show_progress, ...
  )
  
  fd <- functional_distances(
    object = fc, cluster_k = cluster_k, hc_method = hc_method,
    test_factor = test_factor, n_perm = n_perm, test_mode = test_mode,
    min_strata_for_pairwise = min_strata_for_pairwise,
    perm_unit = perm_unit, perm_strata = perm_strata,
    bootstrap = bootstrap, boot_B = boot_B, boot_seed = boot_seed, boot_ci = boot_ci,
    show_progress = show_progress, curve_level = curve_level, ...
  )
  
  if (bootstrap_mode == "refit" && bootstrap) {
    stop("bootstrap_mode='refit' not included in this streamlined version. Use bootstrap_mode='predicted'.")
  }
  
  # Return original structure
  out <- list(
    call = match.call(),
    data = fc$observed_data,
    vars = fc$vars,
    settings = list(
      response_scale = response_scale, eps = eps, min_points = min_points,
      grid_n = grid_n, env_ref = env_ref,
      k_smooth = k_smooth, k_env = k_env, k_trt = k_trt,
      k_smooth_eff = fc$settings$k_smooth_eff,
      k_env_eff = fc$settings$k_env_eff,
      k_trt_eff = fc$settings$k_trt_eff,
      gamma = gamma, discrete = discrete,
      family_try = family_try,
      cluster_k = cluster_k, hc_method = hc_method,
      test_factor = test_factor, n_perm = n_perm,
      test_mode = test_mode,
      perm_unit = perm_unit,
      perm_strata = perm_strata,
      bootstrap = bootstrap, boot_B = boot_B, boot_seed = boot_seed,
      boot_ci = boot_ci, bootstrap_mode = bootstrap_mode,
      min_strata_for_pairwise = min_strata_for_pairwise,
      show_progress = show_progress,
      curve_level = curve_level
    ),
    family_used = fc$family_used,
    gam = fc$gam,
    pred = fc$curves,
    distance = fd$distance_matrix,
    curve_distance = fd$curve_distance,
    test = fd$test,
    hc = fd$hc,
    clusters = fd$clusters,
    scores = fc$scores,
    warnings_betar = fc$warnings_betar,
    bootstrap = fd$bootstrap,
    plot_mean = fc$plot_mean,
    plot_envelope = if(bootstrap && !is.null(fd$bootstrap)) {
      env_tbl <- fd$bootstrap$envelopes
      .time <- fc$vars$time
      .trt <- fc$vars$treatment
      ggplot2::ggplot(
        env_tbl,
        ggplot2::aes(
          x = .data[[.time]],
          y = .data[["mu_boot_mean"]],
          color = .data[[.trt]],
          fill = .data[[.trt]]
        )
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data[["mu_lo"]], ymax = .data[["mu_hi"]]),
          alpha = 0.20,
          colour = NA
        ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::theme_classic()
    } else NULL
  )
  
  class(out) <- "r4pde_compare_curves"
  out
}

