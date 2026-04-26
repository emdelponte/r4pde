#' Fit genotype-specific epidemic trajectories using GAM
#'
#' @description
#' Fits a generalized additive model (GAM) to replicated disease progress curves,
#' returning adjusted treatment mean curves evaluated on a common time grid.
#'
#' @param data A data.frame containing disease progress observations.
#' @param time Character string naming the time variable.
#' @param response Character string naming the response variable.
#' @param treatment Character string naming the treatment/cultivar factor.
#' @param environment Optional character string naming an environment factor.
#' @param block Optional character string naming a blocking factor.
#' @param unit Optional character string naming a unique experimental unit identifier.
#' @param response_scale Character string specifying the response scale: \code{"proportion"} or \code{"percent"}.
#' @param eps Numeric small constant retained for API compatibility.
#' @param min_points Minimum number of observations required per curve.
#' @param grid_n Number of time points for evaluating predicted curves.
#' @param env_ref Reference environment level used for adjusted predictions.
#' @param k_smooth Basis dimension for the global smooth of time.
#' @param k_env Basis dimension for the environment-specific smooth.
#' @param k_trt Basis dimension for the treatment-specific smooth.
#' @param gamma Penalization parameter passed to \code{mgcv::bam()}.
#' @param discrete Logical; whether to use discrete (approximate) fitting.
#' @param family_try Character string specifying the GAM family to try.
#' @param show_progress Logical; whether to show progress.
#' @param covariates Optional character vector of genotype-level covariates (e.g., \code{c("heading_group")}).
#' @param include_covariates Logical; whether to include covariates as fixed effects in the model.
#' @param covariate_smooths Logical; if TRUE and \code{include_covariates = TRUE}, adds smooth interactions \code{s(time, by=covariate)} for factor covariates.
#' @param ... Additional arguments.
#'
#' @details
#' Genotype-level covariates are descriptors that do not vary within a genotype, such as phenological groups.
#' For example, in wheat blast studies, a cultivar's \code{heading_group} (early, intermediate, late)
#' can be supplied. This helps distinguish between true genetic resistance and phenological escape, as
#' cultivars with different heading dates may encounter different infection-risk windows in the same
#' environment. When \code{include_covariates = TRUE}, the model accounts for these covariates,
#' yielding heading-adjusted functional resistance.
#'
#' @return An object of class \code{"functional_curves"} containing:
#' \itemize{
#'   \item \code{gam}: fitted GAM object and \code{family_used};
#'   \item \code{curves}: environment-adjusted treatment mean curves on the common grid;
#'   \item \code{grid}: the time grid used;
#'   \item \code{observed_data}: the processed input data;
#'   \item \code{genotype_info}: a tibble with one row per genotype containing covariates;
#'   \item \code{vars}: variable names used;
#'   \item \code{settings}: model settings;
#'   \item \code{warnings_betar}: any warnings caught during beta fitting;
#'   \item \code{plot_mean}: ggplot object of the mean curves.
#' }
#'
#' @examples
#' \dontrun{
#' fc <- functional_curves(
#'   data = my_data,
#'   time = "time_var",
#'   response = "severity",
#'   treatment = "cultivar",
#'   environment = "env",
#'   covariates = c("heading_group"),
#'   include_covariates = TRUE,
#'   covariate_smooths = TRUE
#' )
#' plot(fc)
#' }
#' @export
functional_curves <- function(
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
    show_progress = TRUE,
    covariates = NULL,
    include_covariates = FALSE,
    covariate_smooths = FALSE,
    ...
) {
  response_scale <- match.arg(response_scale)
  family_try <- match.arg(family_try)

  # Validate columns
  dat <- data
  needed <- c(time, response, treatment)
  if(!is.null(environment)) needed <- c(needed, environment)
  if(!is.null(block))       needed <- c(needed, block)
  if(!is.null(unit))        needed <- c(needed, unit)
  if(!is.null(covariates))  needed <- c(needed, covariates)
  needed <- unique(needed)

  bad <- setdiff(needed, names(dat))
  if(length(bad)) stop("These columns are missing from `data`: ", paste(bad, collapse = ", "))

  .time <- time; .y <- response; .trt <- treatment
  .env  <- environment; .blk <- block; .unit <- unit

  df <- tibble::as_tibble(dat)

  df[[.trt]] <- as.factor(df[[.trt]])
  if(!is.null(.env)) df[[.env]] <- as.factor(df[[.env]])
  if(!is.null(.blk)) df[[.blk]] <- as.factor(df[[.blk]])
  if(!is.null(covariates)) {
    for (cov in covariates) {
      if (is.character(df[[cov]])) df[[cov]] <- as.factor(df[[cov]])
    }
  }
  df[[.time]] <- as.numeric(df[[.time]])

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

  keep <- is.finite(df[[.time]]) &
    is.finite(df$.y) &
    is.finite(df$.y_raw) &
    !is.na(df[[.trt]]) &
    !is.na(df[[.unit]])
  if(!is.null(.env)) keep <- keep & !is.na(df[[.env]])
  if(!is.null(covariates)) {
    for (cov in covariates) {
      keep <- keep & !is.na(df[[cov]])
    }
  }
  df <- df[keep, , drop = FALSE]

  if(!is.null(covariates)) {
    genotype_info <- df |> dplyr::select(dplyr::all_of(c(.trt, covariates))) |> dplyr::distinct()
    if (nrow(genotype_info) > length(unique(df[[.trt]]))) {
      for (cov in covariates) {
        check <- df |> dplyr::group_by(.data[[.trt]]) |> dplyr::summarise(n = dplyr::n_distinct(.data[[cov]]))
        if (any(check$n > 1)) {
          stop(sprintf("Covariate '%s' has multiple values within at least one genotype. Covariates must be genotype-level.", cov))
        }
      }
    }
  } else {
    genotype_info <- tibble::tibble(!!.trt := unique(df[[.trt]]))
  }

  n_by_unit <- df |>
    dplyr::group_by(.data[[.unit]]) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  df <- df |>
    dplyr::left_join(n_by_unit, by = .unit) |>
    dplyr::filter(.data[["n"]] >= min_points) |>
    dplyr::select(-.data[["n"]])

  if(nrow(df) == 0) stop("No data left after filtering. Check `min_points` and missingness.")

  trt_levels0  <- levels(df[[.trt]])
  env_levels0  <- if(!is.null(.env)) levels(df[[.env]]) else NULL
  blk_levels0  <- if(!is.null(.blk)) levels(df[[.blk]]) else NULL
  unit_levels0 <- levels(df[[.unit]])

  k_smooth_eff <- min(k_smooth, max(3, length(unique(df[[.time]])) - 1))
  k_trt_eff    <- min(k_trt,    max(2, length(unique(df[[.trt]]))  - 1))
  k_env_eff <- NULL
  if(!is.null(.env)) k_env_eff <- min(k_env, max(2, length(unique(df[[.env]])) - 1))

  f_terms <- c(
    sprintf("%s", .trt),
    sprintf("s(%s, k=%s)", .time, k_smooth_eff),
    sprintf("s(%s, by=%s, k=%s)", .time, .trt, k_trt_eff)
  )

  if (include_covariates && !is.null(covariates)) {
    for (cov in covariates) {
      f_terms <- c(f_terms, sprintf("%s", cov))
      if (covariate_smooths && is.factor(df[[cov]])) {
        f_terms <- c(f_terms, sprintf("s(%s, by=%s, k=%s)", .time, cov, k_trt_eff))
      }
    }
  }

  if(!is.null(.env)){
    f_terms <- c(f_terms, sprintf("s(%s, bs='re')", .env))
  }

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

  if (family_try == "quasibinomial") {
    m_gam <- fit_qb(df)
    fam_used <- "quasibinomial"
  } else {
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
  }

  tmin <- min(df[[.time]], na.rm = TRUE)
  tmax <- max(df[[.time]], na.rm = TRUE)
  t_grid <- seq(tmin, tmax, length.out = grid_n)

  if(!is.null(.env)){
    if(is.null(env_ref)) env_ref <- env_levels0[1]
    if(!(env_ref %in% env_levels0)) stop("`env_ref` not found in environment levels.")
  }

  pred_trt <- purrr::map_dfr(trt_levels0, function(trt_i){
    newd <- tibble::tibble(tmp_time = t_grid)
    names(newd)[1] <- .time
    newd[[.trt]] <- factor(trt_i, levels = trt_levels0)

    if (!is.null(covariates)) {
      cov_vals <- genotype_info |> dplyr::filter(.data[[.trt]] == trt_i)
      for (cov in covariates) {
        newd[[cov]] <- cov_vals[[cov]][1]
      }
    }

    if(!is.null(.env)) newd[[.env]] <- factor(env_levels0[1], levels = env_levels0)

    excl <- character()
    if(!is.null(.env)) excl <- c(excl, sprintf("s(%s)", .env))

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

  out <- list(
    gam = m_gam,
    curves = pred_trt,
    grid = t_grid,
    observed_data = df,
    genotype_info = genotype_info,
    vars = list(time = .time, response = .y, treatment = .trt,
                environment = .env, block = .blk, unit = unit_used,
                covariates = covariates),
    settings = list(
      response_scale = response_scale, eps = eps, min_points = min_points,
      grid_n = grid_n, env_ref = env_ref,
      k_smooth = k_smooth, k_env = k_env, k_trt = k_trt,
      k_smooth_eff = k_smooth_eff,
      k_env_eff = k_env_eff,
      k_trt_eff = k_trt_eff,
      gamma = gamma, discrete = discrete,
      family_try = family_try
    ),
    warnings_betar = w,
    family_used = fam_used,
    plot_mean = plot_mean,
    scores = trt_score
  )
  class(out) <- "functional_curves"
  out
}

#' Print functional_curves
#' @export
print.functional_curves <- function(x, ...) {
  cat("A functional_curves object\n")
  cat("Time variable:", x$vars$time, "\n")
  cat("Treatment variable:", x$vars$treatment, "\n")
  cat("Grid size:", length(x$grid), "points\n")
  invisible(x)
}

#' Plot functional_curves
#' @export
plot.functional_curves <- function(x, ...) {
  x$plot_mean
}
