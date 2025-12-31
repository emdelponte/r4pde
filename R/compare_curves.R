#' Compare epidemic curves using functional GAM smoothing and clustering
#'
#' @description
#' Fits a generalized additive model (GAM) to disease progress curves using
#' a beta regression (with automatic fallback to quasibinomial),
#' derives environment-adjusted mean curves for each treatment,
#' computes functional distances among treatments,
#' and performs hierarchical clustering of epidemic trajectories.
#'
#' The function is designed for replicated plant disease epidemics
#' measured repeatedly over time, optionally across environments and blocks.
#'
#' @param data A data.frame containing epidemic curve data.
#' @param time Character string giving the name of the time variable.
#' @param response Character string giving the name of the response variable
#'   (disease severity or incidence).
#' @param treatment Character string giving the name of the treatment or cultivar factor.
#' @param environment Optional character string giving the environment factor.
#' @param block Optional character string giving the block or plot factor.
#' @param unit Optional character string giving a unique experimental unit identifier.
#'   If NULL, a unit is constructed from the interaction of environment, treatment,
#'   and block (when available).
#' @param response_scale Character string specifying the response scale:
#'   either \code{"proportion"} (0–1) or \code{"percent"} (0–100).
#' @param eps Numeric small constant retained for API compatibility.
#' @param min_points Minimum number of observations required per epidemic curve.
#' @param grid_n Number of time points used to evaluate predicted mean curves.
#' @param env_ref Reference environment level for adjusted predictions.
#'   If NULL, the first level is used.
#' @param k_smooth Basis dimension for the global smooth over time.
#' @param k_env Basis dimension for the environment-specific smooth.
#' @param k_trt Basis dimension for the treatment-specific smooth.
#' @param gamma Penalization parameter passed to \code{mgcv::bam()}.
#' @param discrete Logical; whether to use discrete fitting in \code{bam()}.
#' @param cluster_k Number of clusters used to cut the hierarchical tree.
#' @param hc_method Clustering method passed to \code{hclust()}.
#' @param ... Reserved for future extensions.
#'
#' @details
#' The response is internally transformed using the Smithson–Verkuilen adjustment
#' to avoid exact 0 and 1 values. Functional distances are computed using an
#' L2 norm over predicted mean curves evaluated on a common time grid.
#'
#' @return An object of class \code{"r4pde_compare_curves"} containing:
#' \itemize{
#'   \item fitted GAM object,
#'   \item environment-adjusted predicted curves,
#'   \item functional distance matrix,
#'   \item hierarchical clustering,
#'   \item treatment cluster assignments.
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
    eps = 1e-4,            # kept for API compatibility; not used for beta handling
    min_points = 5,
    grid_n = 140,
    env_ref = NULL,
    k_smooth = 10,
    k_env = 4,
    k_trt = 4,
    gamma = 1.4,
    discrete = TRUE,
    cluster_k = 4,
    hc_method = "ward.D2",
    ...
){
  response_scale <- match.arg(response_scale)

  # deps
  req <- c("dplyr","tidyr","tibble","purrr","ggplot2","mgcv")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if(length(miss)) stop("Missing required packages: ", paste(miss, collapse = ", "))

  # ---- validate columns ----
  dat <- data
  needed <- c(time, response, treatment)
  if(!is.null(environment)) needed <- c(needed, environment)
  if(!is.null(block))       needed <- c(needed, block)
  if(!is.null(unit))        needed <- c(needed, unit)
  needed <- unique(needed)

  bad <- setdiff(needed, names(dat))
  if(length(bad)) stop("These columns are missing from `data`: ", paste(bad, collapse = ", "))

  .time <- time; .y <- response; .trt <- treatment
  .env  <- environment; .blk <- block; .unit <- unit

  df <- tibble::as_tibble(dat)

  # coerce
  df[[.trt]] <- as.factor(df[[.trt]])
  if(!is.null(.env)) df[[.env]] <- as.factor(df[[.env]])
  if(!is.null(.blk)) df[[.blk]] <- as.factor(df[[.blk]])
  df[[.time]] <- as.numeric(df[[.time]])

  # ---- response handling (robust) ----
  y_raw <- as.numeric(df[[.y]])

  if (response_scale == "percent") {
    y_raw <- y_raw / 100
  } else {
    mx <- suppressWarnings(max(y_raw, na.rm = TRUE))
    if (is.finite(mx) && mx > 1) {
      warning("`response_scale='proportion'` but max(response) > 1; assuming percent and dividing by 100.")
      y_raw <- y_raw / 100
    }
  }

  # constrain to [0,1]
  y_raw <- pmin(pmax(y_raw, 0), 1)

  # Smithson–Verkuilen adjustment (avoid 0/1 for beta regression)
  n_eff <- sum(is.finite(y_raw))
  df$.y <- (y_raw * (n_eff - 1) + 0.5) / n_eff

  # ---- unit id ----
  unit_used <- .unit
  if(is.null(.unit)){
    parts <- list()
    if(!is.null(.env)) parts <- c(parts, list(df[[.env]]))
    parts <- c(parts, list(df[[.trt]]))
    if(!is.null(.blk)) parts <- c(parts, list(df[[.blk]]))
    unit_id <- if(length(parts) == 1) parts[[1]] else interaction(parts, drop = TRUE, sep = "|")
    df$.unit <- as.factor(unit_id)
    .unit <- ".unit"
    unit_used <- ".unit"
  } else {
    df[[.unit]] <- as.factor(df[[.unit]])
  }

  # filter missing essentials
  keep <- is.finite(df[[.time]]) & is.finite(df$.y) & !is.na(df[[.trt]]) & !is.na(df[[.unit]])
  if(!is.null(.env)) keep <- keep & !is.na(df[[.env]])
  df <- df[keep, , drop = FALSE]

  # keep curves with >= min_points
  n_by_unit <- df |>
    dplyr::group_by(.data[[.unit]]) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  df <- df |>
    dplyr::left_join(n_by_unit, by = stats::setNames(.unit, .unit)) |>
    dplyr::filter(.data[["n"]] >= min_points) |>
    dplyr::select(-.data[["n"]])

  if(nrow(df) == 0) stop("No data left after filtering. Check `min_points` and missingness.")

  # ---- model formula (stabilized) ----
  f_terms <- c(
    sprintf("s(%s, k=%s)", .time, k_smooth),
    sprintf("s(%s, %s, bs='fs', k=%s, m=1)", .time, .trt, k_trt),
    sprintf("s(%s, bs='re')", .unit)
  )
  if(!is.null(.env)){
    f_terms <- c(f_terms, sprintf("s(%s, %s, bs='fs', k=%s, m=1)", .time, .env, k_env))
  }
  f_gam <- stats::as.formula(paste0(".y ~ ", paste(f_terms, collapse = " + ")))

  # ---- fit: betar with automatic fallback to quasibinomial ----
  fit_betar <- function() {
    mgcv::bam(
      formula  = f_gam,
      family   = mgcv::betar(link = "logit"),
      data     = df,
      method   = "fREML",
      discrete = discrete,
      gamma    = gamma,
      select   = TRUE
    )
  }

  fit_qb <- function() {
    mgcv::bam(
      formula  = f_gam,
      family   = stats::quasibinomial(link = "logit"),
      data     = df,
      method   = "fREML",
      discrete = discrete,
      gamma    = gamma,
      select   = TRUE
    )
  }

  w <- character()
  m_try <- withCallingHandlers(
    try(fit_betar(), silent = TRUE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )

  theta_failed <- inherits(m_try, "try-error") ||
    any(grepl("theta estimation", w, fixed = TRUE))

  if(theta_failed){
    warning("beta-GAM precision (theta) estimation failed; refitting with quasibinomial().")
    m_gam <- fit_qb()
    fam_used <- "quasibinomial"
  } else {
    m_gam <- m_try
    fam_used <- "betar"
  }

  # ---- env-adjusted mean curves by treatment ----
  tmin <- min(df[[.time]], na.rm = TRUE)
  tmax <- max(df[[.time]], na.rm = TRUE)
  t_grid <- seq(tmin, tmax, length.out = grid_n)

  if(!is.null(.env)){
    env_levels <- levels(df[[.env]])
    if(is.null(env_ref)) env_ref <- env_levels[1]
    if(!(env_ref %in% env_levels)) stop("`env_ref` not found in environment levels.")
  }

  trt_levels <- levels(df[[.trt]])

  pred_trt <- purrr::map_dfr(trt_levels, function(trt_i){
    newd <- tibble::tibble(tmp_time = t_grid)
    names(newd)[1] <- .time
    newd[[.trt]] <- factor(trt_i, levels = trt_levels)

    if(!is.null(.env)){
      newd[[.env]] <- factor(env_ref, levels = levels(df[[.env]]))
    }
    newd[[.unit]] <- df[[.unit]][1]

    # exclude environment smooth and unit RE
    excl <- c(sprintf("s(%s)", .unit))
    if(!is.null(.env)) excl <- c(excl, sprintf("s(%s,%s)", .time, .env))

    mu <- mgcv::predict.gam(m_gam, newdata = newd, type = "response", exclude = excl)

    tibble::tibble(
      treatment = trt_i,
      time = t_grid,
      mu = as.numeric(mu)
    )
  })
  names(pred_trt) <- c(.trt, .time, "mu")

  trt_score <- pred_trt |>
    dplyr::group_by(.data[[.trt]]) |>
    dplyr::summarise(mean_mu = mean(.data[["mu"]]), .groups = "drop") |>
    dplyr::arrange(.data[["mean_mu"]])

  # ---- functional distance + clustering ----
  mat <- pred_trt |>
    dplyr::select(.data[[.time]], .data[[.trt]], .data[["mu"]]) |>
    tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["mu"]]) |>
    dplyr::arrange(.data[[.time]])

  dt_grid <- mean(diff(mat[[.time]]))
  X <- as.matrix(mat[, setdiff(names(mat), .time), drop = FALSE])

  D <- as.matrix(stats::dist(t(X), method = "euclidean")) * sqrt(dt_grid)
  hc <- stats::hclust(stats::as.dist(D), method = hc_method)
  cl_raw <- stats::cutree(hc, k = cluster_k)

  cluster_tbl <- tibble::tibble(
    !!.trt := names(cl_raw),
    cluster_raw = as.integer(cl_raw)
  ) |>
    dplyr::left_join(trt_score, by = stats::setNames(.trt, .trt)) |>
    dplyr::group_by(.data[["cluster_raw"]]) |>
    dplyr::summarise(cluster_mean = mean(.data[["mean_mu"]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data[["cluster_mean"]]) |>
    dplyr::mutate(cluster = dplyr::row_number()) |>
    dplyr::select(.data[["cluster_raw"]], .data[["cluster"]])

  trt_cluster <- tibble::tibble(
    !!.trt := names(cl_raw),
    cluster_raw = as.integer(cl_raw)
  ) |>
    dplyr::left_join(cluster_tbl, by = "cluster_raw") |>
    dplyr::left_join(trt_score, by = stats::setNames(.trt, .trt))

  out <- list(
    call = match.call(),
    data = df,
    vars = list(time = .time, response = .y, treatment = .trt,
                environment = .env, block = .blk, unit = unit_used),
    settings = list(
      response_scale = response_scale, eps = eps, min_points = min_points,
      grid_n = grid_n, env_ref = env_ref,
      k_smooth = k_smooth, k_env = k_env, k_trt = k_trt,
      gamma = gamma, discrete = discrete,
      cluster_k = cluster_k, hc_method = hc_method
    ),
    family_used = fam_used,
    gam = m_gam,
    pred = pred_trt,
    distance = D,
    hc = hc,
    clusters = trt_cluster,
    scores = trt_score,
    warnings_betar = w
  )
  class(out) <- "r4pde_compare_curves"
  out
}
