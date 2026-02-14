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
    cluster_k = 4,
    hc_method = "ward.D2",

    # distance-based permutation test (optional)
    test_factor = NULL,              # column name in data OR named vector by perm_unit id
    n_perm = 999,
    test_mode = c("auto","none","global","global_pairwise"),
    min_strata_for_pairwise = 6,      # e.g., blocks >= 6 before enabling pairwise in auto

    # crucial for phenotyping: what is the "unit of label"?
    # - default NULL => curve-level unit (.unit)
    # - set perm_unit = "cultivar" (or similar) to permute labels at genotype level
    perm_unit = NULL,

    # restricted permutation strata (e.g. blocks, or environments). default: block if provided
    perm_strata = NULL,

    # bootstrap + plots
    bootstrap = FALSE,
    boot_B = 399,
    boot_seed = 1,
    boot_ci = c(0.025, 0.975),
    bootstrap_mode = c("predicted","refit"),  # default predicted (fast)
    show_progress = TRUE,
    curve_level = NULL, 
    ...

){
  response_scale <- match.arg(response_scale)
  test_mode <- match.arg(test_mode)
  bootstrap_mode <- match.arg(bootstrap_mode)

  # deps
  req <- c("dplyr", "tidyr", "tibble", "purrr", "ggplot2", "mgcv", "rlang", "progress")
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

  # ---- response handling ----
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

  y_raw <- pmin(pmax(y_raw, 0), 1)
  df$.y_raw <- y_raw

  n_eff <- sum(is.finite(y_raw))
  if (!is.finite(n_eff) || n_eff < 2) stop("Not enough finite response values to fit the model.")
  df$.y <- (y_raw * (n_eff - 1) + 0.5) / n_eff

  # ---- unit id ----
  unit_used <- .unit
  if(is.null(.unit)){
    parts <- list()
    if(!is.null(.env)) parts <- c(parts, list(df[[.env]]))
    parts <- c(parts, list(df[[.trt]]))
    if(!is.null(.blk)) parts <- c(parts, list(df[[.blk]]))

    unit_id <- if(length(parts) == 1) {
      parts[[1]]
    } else {
      do.call(interaction, c(parts, list(drop = TRUE, sep = "|")))
    }

    df$.unit <- as.factor(unit_id)
    .unit <- ".unit"
    unit_used <- ".unit"
  } else {
    df[[.unit]] <- as.factor(df[[.unit]])
  }

  # filter missing essentials
  keep <- is.finite(df[[.time]]) &
    is.finite(df$.y) &
    is.finite(df$.y_raw) &
    !is.na(df[[.trt]]) &
    !is.na(df[[.unit]])
  if(!is.null(.env)) keep <- keep & !is.na(df[[.env]])
  df <- df[keep, , drop = FALSE]

  # keep curves with >= min_points
  n_by_unit <- df |>
    dplyr::group_by(.data[[.unit]]) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  df <- df |>
    dplyr::left_join(n_by_unit, by = .unit) |>
    dplyr::filter(.data[["n"]] >= min_points) |>
    dplyr::select(-.data[["n"]])

  if(nrow(df) == 0) stop("No data left after filtering. Check `min_points` and missingness.")

  # store original levels (fixed!)
  trt_levels0  <- levels(df[[.trt]])
  env_levels0  <- if(!is.null(.env)) levels(df[[.env]]) else NULL
  blk_levels0  <- if(!is.null(.blk)) levels(df[[.blk]]) else NULL
  unit_levels0 <- levels(df[[.unit]])

  # ---- choose safe k ----
  k_smooth_eff <- min(k_smooth, max(3, length(unique(df[[.time]])) - 1))
  k_trt_eff    <- min(k_trt,    max(2, length(unique(df[[.trt]]))  - 1))
  k_env_eff <- NULL
  if(!is.null(.env)) k_env_eff <- min(k_env, max(2, length(unique(df[[.env]])) - 1))

  # ---- model formula ----
  # Using 'by' smooths + trt main effect is generally more stable and faster than 'fs' 
  f_terms <- c(
    sprintf("%s", .trt),                                 # genotype main effect
    sprintf("s(%s, k=%s)", .time, k_smooth_eff),          # global smooth
    sprintf("s(%s, by=%s, k=%s)", .time, .trt, k_trt_eff) # genotype deviations
  )

  # environment: random intercept (adjustment), not a full time-by-env smooth
  if(!is.null(.env)){
    f_terms <- c(f_terms, sprintf("s(%s, bs='re')", .env))
  }

  # unit RE: only keep if you truly have replication within env×genotype (blocks/plots)
  if(!is.null(.blk) || !is.null(unit)){
    f_terms <- c(f_terms, sprintf("s(%s, bs='re')", .unit))
  }

  f_gam_beta <- stats::as.formula(paste0(".y ~ ", paste(f_terms, collapse = " + ")))
  f_gam_qb   <- stats::as.formula(paste0(".y_raw ~ ", paste(f_terms, collapse = " + ")))

  fit_betar <- function(dat_fit) {
    mgcv::bam(
      formula  = f_gam_beta,
      family   = mgcv::betar(link = "logit"),
      data     = dat_fit,
      method   = "fREML",
      discrete = discrete,
      gamma    = gamma,
      select   = TRUE
    )
  }
  fit_qb <- function(dat_fit) {
    mgcv::bam(
      formula  = f_gam_qb,
      family   = stats::quasibinomial(link = "logit"),
      data     = dat_fit,
      method   = "fREML",
      discrete = discrete,
      gamma    = gamma,
      select   = TRUE
    )
  }

  if (show_progress) message("Fitting GAM model (this may take a moment)...")
  
  w <- character()
  m_try <- withCallingHandlers(
    try(fit_betar(df), silent = TRUE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  theta_failed <- inherits(m_try, "try-error") || any(grepl("theta estimation", w, fixed = TRUE))

  if(theta_failed){
    warning("beta-GAM precision (theta) estimation failed; refitting with quasibinomial().")
    m_gam <- fit_qb(df)
    fam_used <- "quasibinomial"
  } else {
    m_gam <- m_try
    fam_used <- "betar"
  }

  # ---- common time grid ----
  tmin <- min(df[[.time]], na.rm = TRUE)
  tmax <- max(df[[.time]], na.rm = TRUE)
  t_grid <- seq(tmin, tmax, length.out = grid_n)

  if(!is.null(.env)){
    if(is.null(env_ref)) env_ref <- env_levels0[1]
    if(!(env_ref %in% env_levels0)) stop("`env_ref` not found in environment levels.")
  }

  # ---- treatment mean curves (env-adjusted) ----
  # Global adjusted means exclude environment and unit random effects
  pred_trt <- purrr::map_dfr(trt_levels0, function(trt_i){
    newd <- tibble::tibble(tmp_time = t_grid)
    names(newd)[1] <- .time
    newd[[.trt]] <- factor(trt_i, levels = trt_levels0)

    if(!is.null(.env)) newd[[.env]] <- factor(env_levels0[1], levels = env_levels0)

    # exclude random effects for "global adjusted" mean
    excl <- character()
    if(!is.null(.env)) excl <- c(excl, sprintf("s(%s)", .env))

    # if unit RE is in the model, exclude it too
    if(any(grepl(sprintf("s\\(%s\\)", .unit), names(m_gam$smooth), fixed = TRUE))) {
      newd[[.unit]] <- factor(unit_levels0[1], levels = unit_levels0)
      excl <- c(excl, sprintf("s(%s)", .unit))
    }

    mu <- mgcv::predict.gam(m_gam, newdata = newd, type = "response", exclude = excl)

    tibble::tibble(
      !!.trt := trt_i,
      !!.time := t_grid,
      mu = as.numeric(mu)
    )
  })

  trt_score <- pred_trt |>
    dplyr::group_by(.data[[.trt]]) |>
    dplyr::summarise(mean_mu = mean(.data[["mu"]]), .groups = "drop") |>
    dplyr::arrange(.data[["mean_mu"]])

  trt_order <- trt_score[[.trt]]

  plot_mean <- pred_trt |>
    dplyr::mutate(!!.trt := factor(.data[[.trt]], levels = trt_order)) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[.time]], y = .data[["mu"]], color = .data[[.trt]])) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_classic()

  # ---- functional distance + clustering (treatment means) ----
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
    dplyr::left_join(trt_score, by = .trt) |>
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
    dplyr::left_join(trt_score, by = .trt)

  # ------------------------------------------------------------
  # Curve-level fitted trajectories (exclude unit RE) + distance matrix
  # ------------------------------------------------------------
  D_curve <- NULL
  
  # Auto-resolve curve_level
  if (is.null(curve_level)) {
    curve_level <- !is.null(test_factor)
  }

  if (isTRUE(curve_level)) {
    col_name1 <- function(x){
      if (is.null(x)) return(NULL)
      if (is.character(x) && length(x) == 1L) return(x)
      if (rlang::is_quosure(x)) return(rlang::as_name(rlang::quo_get_expr(x)))
      if (rlang::is_symbol(x))  return(rlang::as_name(x))
      if (is.list(x) && length(x) == 1L) return(col_name1(x[[1]]))
      stop("Could not coerce column spec to a single name. Got: ", paste(class(x), collapse = "/"))
    }

    env_col <- col_name1(.env)
    blk_col <- col_name1(.blk)
    ab_cols <- c(env_col, blk_col)
    ab_cols <- ab_cols[!is.na(ab_cols) & nzchar(ab_cols)]

    curve_meta <- df |>
      dplyr::distinct(.data[[unit_used]], .data[[.trt]], dplyr::across(dplyr::all_of(ab_cols))) |>
      dplyr::mutate(
        curve_id_chr = as.character(.data[[unit_used]]),
        trt_chr      = as.character(.data[[.trt]])
      )

    dup <- curve_meta |>
      dplyr::count(curve_id_chr) |>
      dplyr::filter(n > 1)
    if (nrow(dup) > 0) {
      stop("Each curve (unit) must map to one treatment/block/env. Duplicates for: ",
           paste(dup$curve_id_chr, collapse = ", "))
    }

    new_curve_grid <- tidyr::crossing(
      curve_id_chr = curve_meta$curve_id_chr,
      !!.time := t_grid
    ) |>
      dplyr::left_join(curve_meta, by = "curve_id_chr") |>
      dplyr::mutate(
        !!unit_used := factor(curve_id_chr, levels = unit_levels0),
        !!.trt      := factor(trt_chr, levels = trt_levels0)
      )

    if (!is.null(.env) && .env %in% names(new_curve_grid)) {
      new_curve_grid[[.env]] <- factor(new_curve_grid[[.env]], levels = env_levels0)
    }
    if (!is.null(.blk) && .blk %in% names(new_curve_grid)) {
      new_curve_grid[[.blk]] <- factor(new_curve_grid[[.blk]], levels = blk_levels0)
    }

    # Exclude unit RE to get smooth curve-level fitted trajectories
    excl_c <- character()
    if(any(grepl(sprintf("s\\(%s\\)", unit_used), names(m_gam$smooth), fixed = TRUE))) {
       excl_c <- c(excl_c, sprintf("s(%s)", unit_used))
    }

    new_curve_grid$y_hat <- as.numeric(
      mgcv::predict.gam(
        m_gam,
        newdata = new_curve_grid,
        type = "response",
        exclude = excl_c
      )
    )

    Yhat <- new_curve_grid |>
      dplyr::select(curve_id_chr, !!.time, y_hat) |>
      tidyr::pivot_wider(names_from = !!.time, values_from = y_hat) |>
      dplyr::arrange(curve_id_chr)

    curve_ids <- Yhat$curve_id_chr
    Ymat <- as.matrix(Yhat[, -1, drop = FALSE])

    dt <- mean(diff(t_grid))
    D_curve <- as.matrix(stats::dist(Ymat, method = "euclidean")) * sqrt(dt)
    rownames(D_curve) <- curve_ids
    colnames(D_curve) <- curve_ids

    # scaling for numerical stability (keep, but be explicit in interpretation)
    D_curve <- D_curve / stats::median(D_curve[D_curve > 0])
  }

  # ------------------------------------------------------------
  # Permutation test: GLOBAL by default; pairwise only if supported
  # ------------------------------------------------------------
  test_res <- NULL

  if (!is.null(test_factor) && test_mode != "none") {

    # default strata: block if provided
    if (is.null(perm_strata)) perm_strata <- block

    # determine permutation unit (label lives on what?)
    perm_unit_eff <- perm_unit
    if (is.null(perm_unit_eff)) perm_unit_eff <- unit_used

    if (!perm_unit_eff %in% names(df)) {
      stop("`perm_unit` (or default unit) not found in data: ", perm_unit_eff)
    }

    # helper: restricted permutation within strata
    permute_within_strata <- function(g, strata) {
      gperm <- g
      for (s in unique(strata)) {
        idx <- which(strata == s)
        gperm[idx] <- sample(g[idx], length(idx), replace = FALSE)
      }
      gperm
    }

    permanova_F <- function(Dmat, grp) {
      n <- nrow(Dmat)
      grp <- droplevels(factor(grp))
      k <- nlevels(grp)
      if (n < 3 || k < 2) return(list(F = NA_real_, df1 = NA_integer_, df2 = NA_integer_))
      df1 <- k - 1
      df2 <- n - k
      if (df2 <= 0) return(list(F = NA_real_, df1 = df1, df2 = df2))

      A <- -0.5 * (Dmat^2)
      H <- diag(n) - matrix(1, n, n)/n
      G <- H %*% A %*% H
      Xg <- stats::model.matrix(~ grp)
      P  <- Xg %*% solve(t(Xg) %*% Xg) %*% t(Xg)

      SS_between <- sum(diag(P %*% G))
      SS_total   <- sum(diag(G))
      SS_within  <- SS_total - SS_between
      Fstat <- (SS_between/df1) / (SS_within/df2)
      list(F = as.numeric(Fstat), df1 = df1, df2 = df2)
    }

    # Build mapping from perm_unit to group label + (optional) strata label
    if (is.character(test_factor) && length(test_factor) == 1L) {
      tf_col <- test_factor
      if (!tf_col %in% names(df)) stop("`test_factor` not found in data: ", tf_col)

      map_u <- df |>
        dplyr::distinct(.data[[perm_unit_eff]], .data[[tf_col]],
                        .keep_all = FALSE) |>
        dplyr::filter(!is.na(.data[[tf_col]])) |>
        dplyr::mutate(u_chr = as.character(.data[[perm_unit_eff]]))

      dup <- map_u |>
        dplyr::count(u_chr) |>
        dplyr::filter(n > 1)
      if (nrow(dup) > 0) {
        stop("`test_factor` must be unique per `perm_unit`. Duplicates for: ",
             paste(dup$u_chr, collapse = ", "))
      }

      g_u <- map_u[[tf_col]]
      names(g_u) <- map_u$u_chr

    } else {
      # assume named vector by perm_unit
      g_u <- test_factor
      if (is.null(names(g_u))) stop("If `test_factor` is a vector, it must be named by `perm_unit` ids.")
    }

    # Now map each curve_id to perm_unit id, then to group
    map_curve_to_u <- df |>
      dplyr::distinct(.data[[unit_used]], .data[[perm_unit_eff]], .keep_all = FALSE) |>
      dplyr::mutate(curve_id_chr = as.character(.data[[unit_used]]),
                    u_chr = as.character(.data[[perm_unit_eff]]))

    dup2 <- map_curve_to_u |>
      dplyr::count(curve_id_chr) |>
      dplyr::filter(n > 1)
    if (nrow(dup2) > 0) stop("Each curve must map to one perm_unit. Duplicates for: ", paste(dup2$curve_id_chr, collapse = ", "))

    # aligned to D_curve
    curve_ids2 <- rownames(D_curve)
    u_curve <- map_curve_to_u$u_chr
    names(u_curve) <- map_curve_to_u$curve_id_chr
    u_curve <- u_curve[curve_ids2]
    if (any(is.na(u_curve))) stop("Missing perm_unit ids for some curves (units).")

    g_curve <- g_u[u_curve]
    if (any(is.na(g_curve))) {
      miss_u <- unique(u_curve[is.na(g_curve)])
      stop("Missing group labels for perm_unit id(s): ", paste(miss_u, collapse = ", "))
    }
    g_curve <- droplevels(factor(g_curve))
    if (nlevels(g_curve) < 2) stop("`test_factor` must have at least 2 groups.")

    # strata vector aligned to curves (optional)
    strata_curve <- NULL
    if (!is.null(perm_strata)) {
      if (!perm_strata %in% names(df)) stop("`perm_strata` not found in data: ", perm_strata)
      map_curve_to_s <- df |>
        dplyr::distinct(.data[[unit_used]], .data[[perm_strata]], .keep_all = FALSE) |>
        dplyr::mutate(curve_id_chr = as.character(.data[[unit_used]]))
      dup3 <- map_curve_to_s |>
        dplyr::count(curve_id_chr) |>
        dplyr::filter(n > 1)
      if (nrow(dup3) > 0) stop("Each curve must map to one perm_strata. Duplicates for: ", paste(dup3$curve_id_chr, collapse = ", "))

      strata_curve <- map_curve_to_s[[perm_strata]]
      names(strata_curve) <- map_curve_to_s$curve_id_chr
      strata_curve <- droplevels(factor(strata_curve[curve_ids2]))
      if (any(is.na(strata_curve))) stop("Missing strata labels for some curves.")
    }

    # auto mode: disable pairwise when strata levels are too few (e.g., 4 blocks)
    pairwise_allowed <- TRUE
    if (test_mode == "auto") {
      if (!is.null(strata_curve)) {
        n_strata <- nlevels(strata_curve)
        if (n_strata < min_strata_for_pairwise) pairwise_allowed <- FALSE
      }
    }
    if (test_mode == "global") pairwise_allowed <- FALSE
    if (test_mode == "global_pairwise") pairwise_allowed <- TRUE

    if (test_mode == "auto" && !pairwise_allowed) {
      warning("Permutation pairwise tests disabled (auto): limited strata support (e.g., few blocks). Returning global test only.")
    }

    # global
    obs <- permanova_F(D_curve, g_curve)
    set.seed(1)
    
    if (show_progress) {
      pb_glob <- progress::progress_bar$new(
        format = "  Global permutation test [:bar] :percent eta: :eta",
        total = n_perm, clear = FALSE, width = 60, force = TRUE, show_after = 0)
    }
    
    F_perm <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      if (show_progress) pb_glob$tick()
      gperm <- if (is.null(strata_curve)) {
        sample(g_curve, replace = FALSE)
      } else {
        permute_within_strata(g_curve, strata_curve)
      }
      F_perm[i] <- permanova_F(D_curve, gperm)$F
    }
    
    p_global <- (1 + sum(F_perm >= obs$F, na.rm = TRUE)) / (n_perm + 1)

    # pairwise (optional)
    pairwise <- NULL
    lev <- levels(g_curve)
    if (pairwise_allowed && length(lev) > 2) {
      pairs <- utils::combn(lev, 2, simplify = FALSE)
      n_pairs <- length(pairs)
      
      if (show_progress) {
        pb_pair <- progress::progress_bar$new(
          format = "  Pairwise permutation tests [:bar] :percent eta: :eta",
          total = n_pairs * n_perm, clear = FALSE, width = 60, force = TRUE, show_after = 0)
      }

      pairwise <- purrr::map_dfr(pairs, function(pp) {
        keep2 <- g_curve %in% pp
        D2 <- D_curve[keep2, keep2, drop = FALSE]
        g2 <- droplevels(g_curve[keep2])
        s2 <- if (is.null(strata_curve)) NULL else droplevels(strata_curve[keep2])

        ob2 <- permanova_F(D2, g2)
        
        Fp2 <- numeric(n_perm)
        for (j in seq_len(n_perm)) {
          if (show_progress) pb_pair$tick()
          gp <- if (is.null(s2)) {
            sample(g2, replace = FALSE)
          } else {
            permute_within_strata(g2, s2)
          }
          Fp2[j] <- permanova_F(D2, gp)$F
        }
        
        p2 <- (1 + sum(Fp2 >= ob2$F, na.rm = TRUE)) / (n_perm + 1)
        tibble::tibble(group1 = pp[1], group2 = pp[2], F = ob2$F, p = p2)
      }) |>
        dplyr::mutate(p_adj = stats::p.adjust(p, method = "holm"))
    }

    test_res <- list(
      type = "permutation_distance_test",
      mode = test_mode,
      perm_unit = perm_unit_eff,
      perm_strata = perm_strata,
      n_perm = n_perm,
      groups = g_curve,
      strata = strata_curve,
      global = c(F = obs$F, p = p_global, df1 = obs$df1, df2 = obs$df2),
      pairwise = pairwise,
      perm_F = F_perm
    )
  }

  # ------------------------------------------------------------
  # Bootstrap on predicted curves (FAST + correct w/ replacement)
  # ------------------------------------------------------------
  boot_res <- NULL
  plot_envelope <- NULL

  if (isTRUE(bootstrap)) {

    boot_ci <- sort(boot_ci)
    if (length(boot_ci) != 2 || any(!is.finite(boot_ci)) || boot_ci[1] <= 0 || boot_ci[2] >= 1) {
      stop("`boot_ci` must be a length-2 numeric vector strictly between 0 and 1, e.g. c(0.025, 0.975).")
    }
    if (!is.numeric(boot_B) || boot_B < 50) stop("`boot_B` should be >= 50 (e.g., 199, 399, 999).")

    if (bootstrap_mode == "refit") {
      warning("bootstrap_mode='refit' is slow (fits bam() boot_B times). Prefer bootstrap_mode='predicted'.")
    }

    if (bootstrap_mode == "predicted") {

      # resampling unit: blocks if available, else curves (units)
      set.seed(boot_seed)
      if(!is.null(.blk)) {
        ids0 <- blk_levels0
        id_col <- .blk
        boot_unit_label <- "block"
      } else {
        ids0 <- unit_levels0
        id_col <- unit_used
        boot_unit_label <- "unit"
      }

      if (show_progress) {
        pb_boot <- progress::progress_bar$new(
          format = "  Bootstrapping predictions [:bar] :percent eta: :eta",
          total = boot_B, clear = FALSE, width = 60, force = TRUE, show_after = 0)
      }

      for (b in seq_len(boot_B)) {
        if (show_progress) pb_boot$tick()
        ids_b <- sample(ids0, replace = TRUE)

        w_draw <- tibble::tibble(id = factor(ids_b, levels = ids0)) |>
          dplyr::count(id, name = "w")

        grid_b <- new_curve_grid |>
          dplyr::mutate(id = factor(.data[[id_col]], levels = ids0)) |>
          dplyr::left_join(w_draw, by = "id") |>
          dplyr::filter(!is.na(.data[["w"]]))

        boot_pred_list[[b]] <- grid_b |>
          dplyr::group_by(.data[[.trt]], .data[[.time]]) |>
          dplyr::summarise(
            mu = stats::weighted.mean(.data[["y_hat"]], w = .data[["w"]], na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::mutate(.iter = b)
      }

      boot_pred <- dplyr::bind_rows(boot_pred_list)

      env_tbl <- boot_pred |>
        dplyr::group_by(.data[[.trt]], .data[[.time]]) |>
        dplyr::summarise(
          mu_boot_mean = mean(.data[["mu"]], na.rm = TRUE),
          mu_lo = stats::quantile(.data[["mu"]], probs = boot_ci[1], na.rm = TRUE),
          mu_hi = stats::quantile(.data[["mu"]], probs = boot_ci[2], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(!!.trt := factor(.data[[.trt]], levels = trt_order))

      dt_bt <- mean(diff(t_grid))

      boot_dist_one_iter <- function(df_iter) {
        mat_i <- df_iter |>
          dplyr::select(.data[[.time]], .data[[.trt]], .data[["mu"]]) |>
          tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["mu"]]) |>
          dplyr::arrange(.data[[.time]])

        X_i <- as.matrix(mat_i[, setdiff(names(mat_i), .time), drop = FALSE])
        coln <- colnames(X_i)
        if (ncol(X_i) < 2) return(NULL)

        pairs <- utils::combn(coln, 2, simplify = FALSE)
        purrr::map_dfr(pairs, function(pp){
          d <- sqrt(sum((X_i[, pp[1]] - X_i[, pp[2]])^2) * dt_bt)
          tibble::tibble(group1 = pp[1], group2 = pp[2], dist = as.numeric(d))
        })
      }

      if (show_progress) {
        pb_dist <- progress::progress_bar$new(
          format = "  Calculating bootstrap distances [:bar] :percent eta: :eta",
          total = boot_B, clear = FALSE, width = 60, force = TRUE, show_after = 0)
      }

      dist_boot <- boot_pred |>
        dplyr::group_by(.iter) |>
        dplyr::group_modify(~ {
          if (show_progress) pb_dist$tick()
          boot_dist_one_iter(.x)
        }) |>
        dplyr::ungroup()

      dist_summary <- dist_boot |>
        dplyr::group_by(.data[["group1"]], .data[["group2"]]) |>
        dplyr::summarise(
          dist_mean = mean(.data[["dist"]], na.rm = TRUE),
          dist_med  = stats::median(.data[["dist"]], na.rm = TRUE),
          dist_lo   = stats::quantile(.data[["dist"]], probs = boot_ci[1], na.rm = TRUE),
          dist_hi   = stats::quantile(.data[["dist"]], probs = boot_ci[2], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(.data[["dist_mean"]]))

      plot_envelope <- ggplot2::ggplot(
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

      boot_res <- list(
        method = "bootstrap_on_predicted_curves",
        unit = boot_unit_label,
        id_col = id_col,
        B = boot_B,
        seed = boot_seed,
        ci = boot_ci,
        envelopes = env_tbl,
        distances = dist_summary,
        distances_boot = dist_boot
      )
    }

    # If you still want the old slow refit bootstrap sometimes:
    if (bootstrap_mode == "refit") {
      # keep minimal placeholder: you can paste your old refit bootstrap here if needed
      stop("bootstrap_mode='refit' not included in this streamlined version. Use bootstrap_mode='predicted'.")
    }
  }

  out <- list(
    call = match.call(),
    data = df,
    vars = list(time = .time, response = .y, treatment = .trt,
                environment = .env, block = .blk, unit = unit_used),
    settings = list(
      response_scale = response_scale, eps = eps, min_points = min_points,
      grid_n = grid_n, env_ref = env_ref,
      k_smooth = k_smooth, k_env = k_env, k_trt = k_trt,
      k_smooth_eff = k_smooth_eff,
      k_env_eff = k_env_eff,
      k_trt_eff = k_trt_eff,
      gamma = gamma, discrete = discrete,
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
    family_used = fam_used,
    gam = m_gam,
    pred = pred_trt,
    distance = D,
    curve_distance = D_curve,
    test = test_res,
    hc = hc,
    clusters = trt_cluster,
    scores = trt_score,
    warnings_betar = w,
    bootstrap = boot_res,
    plot_mean = plot_mean,
    plot_envelope = plot_envelope
  )
  class(out) <- "r4pde_compare_curves"
  out
}

