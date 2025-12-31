#' Diagnostic tools for functional epidemic curve models
#'
#' @description
#' Computes residuals, fitted values, and model diagnostics for
#' GAM-based epidemic curve models fitted with \code{compare_curves()}.
#' Produces diagnostic plots without invoking base graphics.
#'
#' @param x An object of class \code{"r4pde_compare_curves"}.
#' @param grid_n Number of points used for the diagnostic smooth curve.
#'
#' @return An object of class \code{"r4pde_curve_diagnostics"} containing:
#' \itemize{
#'   \item diagnostic data,
#'   \item k-index checks,
#'   \item concurvity estimates,
#'   \item ggplot-based diagnostic panels.
#' }
#'
#' @export
diagnose_curves <- function(x, grid_n = 200){
  stopifnot(inherits(x, "r4pde_compare_curves"))

  req <- c("mgcv","ggplot2","dplyr","tibble")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if(length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))

  df  <- x$data
  tim <- x$vars$time
  trt <- x$vars$treatment

  # fitted values and residuals (ensure same length as df)
  mu_hat   <- as.numeric(mgcv::predict.gam(x$gam, type = "response"))
  res_pear <- as.numeric(stats::residuals(x$gam, type = "pearson"))
  res_dev  <- suppressWarnings(as.numeric(stats::residuals(x$gam, type = "deviance")))

  if(length(mu_hat) != nrow(df) || length(res_pear) != nrow(df)) {
    stop("Diagnostics length mismatch: fitted/residuals not aligned with x$data (nrow).")
  }
  if(length(res_dev) != nrow(df)) res_dev <- rep(NA_real_, nrow(df))

  diag_df <- dplyr::mutate(
    df,
    .mu_hat   = mu_hat,
    .res_pear = res_pear,
    .res_dev  = res_dev
  )

  summ <- list(
    family_used = x$family_used %||% NA_character_,
    n = nrow(df),
    n_treatments = nlevels(df[[trt]]),
    n_units = nlevels(df[[x$vars$unit]]),
    time_range = range(df[[tim]], na.rm = TRUE),
    y_range = range(df$.y, na.rm = TRUE),
    edf = sum(x$gam$edf),
    aic = tryCatch(stats::AIC(x$gam), error = function(e) NA_real_),
    warnings_betar = x$warnings_betar %||% character()
  )

  kcheck <- tryCatch(
    suppressWarnings(mgcv::k.check(x$gam)),
    error = function(e) NULL
  )

  conc <- tryCatch(
    suppressWarnings(mgcv::concurvity(x$gam, estimate = "worst")),
    error = function(e) NULL
  )

  # ---- plots ----
  p_fit <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .mu_hat, y = .y)) +
    ggplot2::geom_point(alpha = 0.35) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(x = "Fitted mean", y = "Observed (transformed proportion)") +
    ggplot2::theme_classic(base_size = 12)

  p_res_time <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .data[[tim]], y = .res_pear)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.35) +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::labs(x = tim, y = "Pearson residual") +
    ggplot2::theme_classic(base_size = 12)

  p_res_fit <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .mu_hat, y = .res_pear)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.35) +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::labs(x = "Fitted mean", y = "Pearson residual") +
    ggplot2::theme_classic(base_size = 12)

  p_qq <- ggplot2::ggplot(
    diag_df |> dplyr::filter(is.finite(.res_dev)),
    ggplot2::aes(sample = .res_dev)
  ) +
    ggplot2::stat_qq(alpha = 0.35) +
    ggplot2::stat_qq_line() +
    ggplot2::labs(x = "Theoretical quantiles", y = "Deviance residuals") +
    ggplot2::theme_classic(base_size = 12)

  # ---- sanity curve ----
  trt0 <- levels(df[[trt]])[1]
  t_grid <- seq(min(df[[tim]]), max(df[[tim]]), length.out = grid_n)

  newd <- tibble::tibble(.tmp = t_grid)
  names(newd)[1] <- tim

  newd[[trt]] <- factor(trt0, levels = levels(df[[trt]]))

  if(!is.null(x$vars$environment)){
    env <- x$vars$environment
    env_ref <- x$settings$env_ref %||% levels(df[[env]])[1]
    newd[[env]] <- factor(env_ref, levels = levels(df[[env]]))
  }

  newd[[x$vars$unit]] <- df[[x$vars$unit]][1]

  excl <- c(sprintf("s(%s)", x$vars$unit))
  if(!is.null(x$vars$environment)){
    excl <- c(excl, sprintf("s(%s,%s)", tim, x$vars$environment))
  }

  mu_grid <- mgcv::predict.gam(
    x$gam, newdata = newd, type = "response", exclude = excl
  )

  p_curve <- ggplot2::ggplot(
    tibble::tibble(time = t_grid, mu = as.numeric(mu_grid)),
    ggplot2::aes(x = time, y = mu)
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::expand_limits(y = c(0, 1)) +
    ggplot2::labs(
      x = tim,
      y = paste0("Adjusted mean (", trt0, ")"),
      title = "Sanity check: one adjusted curve"
    ) +
    ggplot2::theme_classic(base_size = 12)

  out <- list(
    summary = summ,
    kcheck = kcheck,
    concurvity = conc,
    diag_data = diag_df,
    plots = list(
      fitted_vs_observed = p_fit,
      residuals_vs_time = p_res_time,
      residuals_vs_fitted = p_res_fit,
      qq_deviance = p_qq,
      one_adjusted_curve = p_curve
    )
  )
  class(out) <- "r4pde_curve_diagnostics"
  out
}

#' Plot diagnostics for functional epidemic curve models
#'
#' @description
#' Combines multiple diagnostic plots from
#' \code{\link{diagnose_curves}} into a multi-panel layout.
#'
#' @param diag An object of class \code{"r4pde_curve_diagnostics"}.
#' @param layout Layout of diagnostic panels: \code{"2x2"}, \code{"vertical"},
#'   or \code{"horizontal"}.
#' @param show_curve Logical; whether to include the example adjusted curve.
#' @param rel_heights Relative heights for panels when stacking.
#'
#' @return A composite plot object generated using \pkg{cowplot}.
#'
#' @export
plot_diagnostics <- function(
    diag,
    layout = c("2x2", "vertical", "horizontal"),
    show_curve = FALSE,
    rel_heights = c(1, 0.7)
){
  stopifnot(inherits(diag, "r4pde_curve_diagnostics"))

  req <- c("cowplot")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if(length(miss)) stop("Missing package: ", paste(miss, collapse = ", "))

  layout <- match.arg(layout)

  p_fit  <- diag$plots$fitted_vs_observed
  p_rt   <- diag$plots$residuals_vs_time
  p_rf   <- diag$plots$residuals_vs_fitted
  p_qq   <- diag$plots$qq_deviance
  p_curv <- diag$plots$one_adjusted_curve

  subtitle <- paste0(
    "Family: ", diag$summary$family_used,
    " | EDF = ", signif(diag$summary$edf, 4),
    " | AIC = ", signif(diag$summary$aic, 6)
  )

  title_g <- cowplot::ggdraw() +
    cowplot::draw_label("Model diagnostics", x = 0, hjust = 0, fontface = "bold", size = 13) +
    cowplot::draw_label(subtitle, x = 0, y = 0.25, hjust = 0, size = 10)

  fig <- switch(
    layout,

    "2x2" = {
      top <- cowplot::plot_grid(p_fit, p_rt, nrow = 1, align = "hv")
      bot <- cowplot::plot_grid(p_rf,  p_qq, nrow = 1, align = "hv")
      core <- cowplot::plot_grid(top, bot, ncol = 1, align = "hv")

      if(show_curve){
        cowplot::plot_grid(core, p_curv, ncol = 1, rel_heights = rel_heights, align = "hv")
      } else core
    },

    "vertical" = {
      if(show_curve){
        cowplot::plot_grid(p_fit, p_rt, p_rf, p_qq, p_curv, ncol = 1, align = "hv")
      } else {
        cowplot::plot_grid(p_fit, p_rt, p_rf, p_qq, ncol = 1, align = "hv")
      }
    },

    "horizontal" = {
      if(show_curve){
        cowplot::plot_grid(p_fit, p_rt, p_rf, p_qq, p_curv, nrow = 1, align = "hv")
      } else {
        cowplot::plot_grid(p_fit, p_rt, p_rf, p_qq, nrow = 1, align = "hv")
      }
    }
  )

  cowplot::plot_grid(title_g, fig, ncol = 1, rel_heights = c(0.13, 1))
}

#' Print method for curve diagnostics
#'
#' @param x An object of class \code{"r4pde_curve_diagnostics"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.r4pde_curve_diagnostics <- function(x, ...){
  s <- x$summary
  cat("Curve model diagnostics\n")
  cat("  family_used: ", s$family_used, "\n", sep = "")
  cat("  n (rows): ", s$n, "\n", sep = "")
  cat("  treatments: ", s$n_treatments, "\n", sep = "")
  cat("  units: ", s$n_units, "\n", sep = "")
  cat("  time range: ", paste(s$time_range, collapse = " to "), "\n", sep = "")
  cat("  y range: ", paste(signif(s$y_range, 3), collapse = " to "), "\n", sep = "")
  cat("  EDF sum: ", signif(s$edf, 4), "\n", sep = "")
  cat("  AIC: ", signif(s$aic, 6), "\n", sep = "")
  if(length(s$warnings_betar)) {
    cat("  captured warnings (betar attempt):\n")
    for(w in s$warnings_betar) cat("   - ", w, "\n", sep = "")
  }
  invisible(x)
}
