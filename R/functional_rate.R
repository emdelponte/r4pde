#' Instantaneous Disease Progress Rates from Functional Curves
#'
#' @description
#' Estimates instantaneous rates of plant disease progress (\eqn{S'(t) = dS(t)/dt})
#' and associated uncertainty from epidemic trajectories fitted by
#' \code{\link{functional_curves}}. Derivatives are computed using the GAM linear
#' predictor matrix (\code{lpmatrix}) with finite differences, supporting both
#' the original response scale (e.g., severity proportion/percentage points per day)
#' and the link (linear predictor) scale. Also computes epidemic rate phenotypes
#' including maximum rate (\eqn{r_{max}}), time of maximum rate (\eqn{t_{r_{max}}}),
#' growth duration, cumulative positive growth, and distinct growth windows.
#'
#' @param object An object of class \code{"functional_curves"}.
#' @param scale Character specifying the derivative scale: \code{"response"}
#'   (default; e.g., severity change per unit time) or \code{"link"} (linear predictor scale).
#' @param method Character specifying the finite difference method. Currently \code{"central"}
#'   (central difference in the interior with one-sided differences at temporal boundaries).
#' @param n_grid Integer; number of points in the regular temporal evaluation grid per group.
#'   Default is 200.
#' @param eps Numeric step size used in finite differences. If \code{NULL} (default),
#'   automatically chosen as \eqn{10^{-4} \times (t_{max} - t_{min})}.
#' @param interval Numeric confidence level for intervals. Default is 0.95.
#' @param threshold Numeric minimum rate threshold considered biologically relevant.
#'   Default is 0.
#' @param positive_only Logical; whether to restrict summarized growth descriptors
#'   to positive rates. Default is \code{TRUE}. Note that negative rates are always
#'   preserved in the detailed rate table (\code{data}).
#' @param uncertainty Logical; whether to compute standard errors and confidence
#'   intervals using the GAM covariance matrix. Default is \code{TRUE}.
#' @param growth_criterion Character specifying how active epidemic growth is defined:
#'   \code{"rate"} (default; \code{rate > threshold}) or \code{"significant"}
#'   (\code{rate_lower > threshold}, requiring \code{uncertainty = TRUE}).
#' @param newdata Optional user-supplied data frame with custom evaluation points.
#'   Must contain the time and treatment variables.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' \strong{Mathematical Formulation:}
#' Let \eqn{S(t)} denote the fitted disease trajectory and \eqn{\eta(t) = X(t)\beta}
#' the linear predictor from the GAM model, with \eqn{\mu(t) = g^{-1}(\eta(t))}.
#'
#' The derivative on the link scale is estimated by:
#' \deqn{D = \frac{X(t + \epsilon) - X(t - \epsilon)}{2\epsilon}}
#' \deqn{\eta'(t) = D \beta}
#' with variance \eqn{\operatorname{Var}(\eta'(t)) = \operatorname{diag}(D V_p D^T)},
#' where \eqn{V_p} is the Bayesian posterior covariance matrix of the GAM parameters.
#'
#' On the response scale (\code{scale = "response"}), the derivative of the mean trajectory is:
#' \deqn{S'(t) = \frac{\mu(\eta(t + \epsilon)) - \mu(\eta(t - \epsilon))}{2\epsilon}}
#' The approximate gradient vector \eqn{G} with respect to \eqn{\beta} is:
#' \deqn{G = \frac{\mu'(\eta(t + \epsilon)) X(t + \epsilon) - \mu'(\eta(t - \epsilon)) X(t - \epsilon)}{2\epsilon}}
#' where \eqn{\mu'(\eta) = d\mu/d\eta} is given by \code{family$mu.eta()}.
#' The uncertainty is then:
#' \deqn{\operatorname{Var}(S'(t)) = \operatorname{diag}(G V_p G^T)}
#'
#' \strong{Boundary Handling:}
#' At the observed boundaries \eqn{t_{min}} and \eqn{t_{max}}, one-sided forward and backward
#' differences are used to prevent unintended extrapolation beyond the observed domain.
#'
#' \strong{Time-Invariant Terms and Random Effects:}
#' Intercepts, main treatment effects, and terms constant with respect to time cancel
#' analytically in \eqn{X(t + \epsilon) - X(t - \epsilon)}. Environmental and experimental
#' unit random effects are excluded by default to yield population-adjusted treatment rates,
#' consistent with \code{\link{functional_curves}}.
#'
#' \strong{No-Epidemic and Flat Curves:}
#' Curves with all observed response values equal to zero (or estimated rates below
#' numerical tolerance \eqn{10^{-6}}) explicitly return \eqn{r_{max} = 0}, \eqn{t_{r_{max}} = \text{NA}},
#' \eqn{\text{growth\_duration} = 0}, and \eqn{\text{cumulative\_positive\_growth} = 0}.
#'
#' @return An object of class \code{"functional_rate"} containing:
#' \describe{
#'   \item{\code{data}}{A tibble of full evaluation grid with columns including
#'     treatment, time, \code{fitted}, \code{fitted_se}, \code{fitted_lower}, \code{fitted_upper},
#'     \code{rate}, \code{rate_se}, \code{rate_lower}, \code{rate_upper},
#'     \code{significant_increase}, \code{above_threshold}, and \code{significant_above_threshold}.}
#'   \item{\code{summary}}{A tibble of curve-level rate phenotypes:
#'     \code{r_max}, \code{t_r_max}, \code{r_mean}, \code{r_mean_positive},
#'     \code{t_growth_start}, \code{t_growth_end}, \code{growth_duration},
#'     \code{cumulative_positive_growth}, and \code{n_growth_windows}.}
#'   \item{\code{model}}{The underlying GAM model object.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{settings}}{List of analysis settings used.}
#'   \item{\code{vars}}{List of variable names identified in the model.}
#'   \item{\code{family}}{Character string naming the GAM family.}
#' }
#'
#' @references
#' Madden, L. V., Hughes, G., & van den Bosch, F. (2007).
#' The Study of Plant Disease Epidemics. APS Press, St. Paul, MN.
#'
#' Wood, S. N. (2017). Generalized Additive Models: An Introduction with R (2nd ed.).
#' Chapman and Hall/CRC.
#'
#' Simpson, G. L. (2018). Modelling Palaeoecological Time Series Using Generalized Additive Models.
#' Frontiers in Ecology and Evolution, 6, 149.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' df <- data.frame(
#'   time = rep(1:30, 3),
#'   treatment = rep(c("Early", "Late", "None"), each = 30),
#'   y = c(
#'     1 / (1 + exp(-(1:30 - 10) * 0.4)),  # Early
#'     1 / (1 + exp(-(1:30 - 20) * 0.4)),  # Late
#'     rep(0, 30)                          # None
#'   )
#' )
#'
#' curves <- functional_curves(
#'   data = df, time = "time", response = "y", treatment = "treatment",
#'   min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
#' )
#'
#' rates <- functional_rate(curves)
#' print(rates)
#' summary(rates)
#' plot(rates)
#' }
#'
#' @export
functional_rate <- function(
    object,
    scale = c("response", "link"),
    method = c("central"),
    n_grid = 200,
    eps = NULL,
    interval = 0.95,
    threshold = 0,
    positive_only = TRUE,
    uncertainty = TRUE,
    growth_criterion = c("rate", "significant"),
    newdata = NULL,
    ...
) {
  scale <- match.arg(scale)
  method <- match.arg(method)
  growth_criterion <- match.arg(growth_criterion)

  # Validation of input object
  if (!inherits(object, "functional_curves")) {
    stop("`object` must be of class 'functional_curves'.", call. = FALSE)
  }

  m_gam <- object$gam
  if (is.null(m_gam)) {
    stop("The GAM model is missing from `object`.", call. = FALSE)
  }

  .time <- object$vars$time
  .trt  <- object$vars$treatment
  .env  <- object$vars$environment
  .unit <- object$vars$unit
  covariates <- object$vars$covariates

  if (is.null(.time) || !is.character(.time) || length(.time) != 1) {
    stop("Time variable is not identifiable in `object`.", call. = FALSE)
  }
  if (is.null(.trt) || !is.character(.trt) || length(.trt) != 1) {
    stop("Treatment variable is not identifiable in `object`.", call. = FALSE)
  }

  # Validate arguments
  if (!is.numeric(n_grid) || length(n_grid) != 1 || is.na(n_grid) || n_grid < 2) {
    stop("`n_grid` must be a single integer >= 2.", call. = FALSE)
  }
  n_grid <- as.integer(n_grid)

  if (!is.null(eps)) {
    if (!is.numeric(eps) || length(eps) != 1 || is.na(eps) || eps <= 0) {
      stop("`eps` must be a single positive numeric value.", call. = FALSE)
    }
  }

  if (!is.numeric(interval) || length(interval) != 1 || is.na(interval) ||
      interval <= 0 || interval >= 1) {
    stop("`interval` must be a single numeric value in (0, 1).", call. = FALSE)
  }

  obs_df <- object$observed_data
  if (is.null(obs_df) && is.null(newdata)) {
    stop("`observed_data` is missing from `object` and `newdata` was not supplied.", call. = FALSE)
  }

  # Validate time domain in observed data if available
  if (!is.null(obs_df)) {
    times_by_trt <- obs_df |>
      dplyr::group_by(.data[[.trt]]) |>
      dplyr::summarise(
        n_times = dplyr::n_distinct(.data[[.time]]),
        t_min = min(.data[[.time]], na.rm = TRUE),
        t_max = max(.data[[.time]], na.rm = TRUE),
        .groups = "drop"
      )

    if (any(times_by_trt$n_times < 2)) {
      bad_trts <- times_by_trt |> dplyr::filter(.data$n_times < 2) |> dplyr::pull(.data[[.trt]])
      stop("The following groups have fewer than two distinct observed time points: ",
           paste(bad_trts, collapse = ", "), call. = FALSE)
    }

    if (any(times_by_trt$n_times < 4)) {
      warning("Some groups have very few (< 4) distinct observed time points; rate estimates may be unstable.",
              call. = FALSE)
    }
  }

  # Identify treatments with all-zero observed responses
  zero_epidemic_trts <- character(0)
  if (!is.null(obs_df)) {
    resp_col <- if (".y_raw" %in% names(obs_df)) {
      ".y_raw"
    } else if (!is.null(object$vars$response) && object$vars$response %in% names(obs_df)) {
      object$vars$response
    } else {
      NULL
    }

    if (!is.null(resp_col)) {
      zero_check <- obs_df |>
        dplyr::group_by(.data[[.trt]]) |>
        dplyr::summarise(
          all_zero = all(abs(.data[[resp_col]]) < 1e-6, na.rm = TRUE),
          .groups = "drop"
        )
      zero_epidemic_trts <- as.character(zero_check[[.trt]][zero_check$all_zero])
    }
  }

  # Extract model levels
  trt_var_levels <- if (!is.null(obs_df)) {
    levels(factor(obs_df[[.trt]]))
  } else {
    levels(factor(object$curves[[.trt]]))
  }
  env_var_levels <- if (!is.null(.env) && !is.null(obs_df)) levels(factor(obs_df[[.env]])) else NULL
  unit_var_levels <- if (!is.null(.unit) && !is.null(obs_df) && .unit %in% names(obs_df)) {
    levels(factor(obs_df[[.unit]]))
  } else {
    NULL
  }

  # Construct evaluation grid
  if (is.null(newdata)) {
    genotype_info <- object$genotype_info

    grid_list <- lapply(trt_var_levels, function(trt_i) {
      sub_obs <- if (!is.null(obs_df)) {
        obs_df[obs_df[[.trt]] == trt_i, , drop = FALSE]
      } else {
        object$curves[object$curves[[.trt]] == trt_i, , drop = FALSE]
      }

      tmin_i <- min(sub_obs[[.time]], na.rm = TRUE)
      tmax_i <- max(sub_obs[[.time]], na.rm = TRUE)

      if (!is.finite(tmin_i) || !is.finite(tmax_i) || tmax_i - tmin_i < 1e-6) {
        warning(sprintf("Group '%s' has an extremely short time span (span < 1e-6).", trt_i), call. = FALSE)
      }

      t_grid_i <- seq(tmin_i, tmax_i, length.out = n_grid)

      df_i <- tibble::tibble(
        !!.trt := factor(trt_i, levels = trt_var_levels),
        !!.time := t_grid_i,
        .t_min = tmin_i,
        .t_max = tmax_i
      )

      if (!is.null(covariates) && !is.null(genotype_info)) {
        cov_vals <- genotype_info |> dplyr::filter(.data[[.trt]] == trt_i)
        for (cov in covariates) {
          if (cov %in% names(cov_vals)) {
            df_i[[cov]] <- cov_vals[[cov]][1]
          }
        }
      }

      if (!is.null(.env)) {
        ref_env <- if (!is.null(object$settings$env_ref)) object$settings$env_ref else env_var_levels[1]
        df_i[[.env]] <- factor(ref_env, levels = env_var_levels)
      }

      if (!is.null(.unit) && !is.null(unit_var_levels) && length(unit_var_levels) > 0) {
        df_i[[.unit]] <- factor(unit_var_levels[1], levels = unit_var_levels)
      }

      df_i
    })

    grid_df <- dplyr::bind_rows(grid_list)

  } else {
    # Custom newdata provided
    req_cols <- c(.time, .trt)
    if (!is.null(.env)) req_cols <- c(req_cols, .env)
    if (!is.null(covariates)) req_cols <- c(req_cols, covariates)

    missing_cols <- setdiff(req_cols, names(newdata))
    if (length(missing_cols) > 0) {
      stop("`newdata` is missing required variables: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    grid_df <- tibble::as_tibble(newdata)
    grid_df[[.trt]] <- factor(grid_df[[.trt]], levels = trt_var_levels)
    grid_df[[.time]] <- as.numeric(grid_df[[.time]])

    if (!is.null(.env)) {
      grid_df[[.env]] <- factor(grid_df[[.env]], levels = env_var_levels)
    }
    if (!is.null(.unit) && !is.null(unit_var_levels) && length(unit_var_levels) > 0) {
      if (!(.unit %in% names(grid_df))) {
        grid_df[[.unit]] <- factor(unit_var_levels[1], levels = unit_var_levels)
      } else {
        grid_df[[.unit]] <- factor(grid_df[[.unit]], levels = unit_var_levels)
      }
    }

    # Check for extrapolation
    if (!is.null(obs_df)) {
      t_bounds <- obs_df |>
        dplyr::group_by(.data[[.trt]]) |>
        dplyr::summarise(
          t_min = min(.data[[.time]], na.rm = TRUE),
          t_max = max(.data[[.time]], na.rm = TRUE),
          .groups = "drop"
        )
      grid_df <- dplyr::left_join(grid_df, t_bounds, by = .trt)

      # Check bounds
      out_of_bounds <- grid_df[[.time]] < grid_df$t_min - 1e-7 | grid_df[[.time]] > grid_df$t_max + 1e-7
      if (any(out_of_bounds, na.rm = TRUE)) {
        stop("`newdata` contains time points outside the observed domain (extrapolation is not permitted).",
             call. = FALSE)
      }
      grid_df <- dplyr::rename(grid_df, .t_min = "t_min", .t_max = "t_max")
    } else {
      grid_df$.t_min <- min(grid_df[[.time]], na.rm = TRUE)
      grid_df$.t_max <- max(grid_df[[.time]], na.rm = TRUE)
    }
  }

  # Build perturbation grids for finite differences
  t_vals <- grid_df[[.time]]
  t_min_v <- grid_df$.t_min
  t_max_v <- grid_df$.t_max
  span_v <- pmax(t_max_v - t_min_v, 1e-6)

  eps_v <- if (!is.null(eps)) rep(eps, length(t_vals)) else 1e-4 * span_v
  eps_v <- pmin(eps_v, span_v / 4)

  # Boundary identification:
  # Left boundary: t - eps < t_min
  # Right boundary: t + eps > t_max
  # Interior: central difference
  is_left  <- (t_vals - eps_v < t_min_v - 1e-9)
  is_right <- (t_vals + eps_v > t_max_v + 1e-9)

  t1_v <- ifelse(is_left, t_vals, ifelse(is_right, t_vals - eps_v, t_vals - eps_v))
  t2_v <- ifelse(is_left, t_vals + eps_v, ifelse(is_right, t_vals, t_vals + eps_v))
  h_v  <- t2_v - t1_v

  # Data frames for GAM prediction
  df_0 <- grid_df |> dplyr::select(-".t_min", -".t_max")
  df_1 <- df_0
  df_2 <- df_0
  df_1[[.time]] <- t1_v
  df_2[[.time]] <- t2_v

  # Exclude random effect terms from prediction
  excl <- character()
  if (!is.null(.env)) excl <- c(excl, sprintf("s(%s)", .env))

  smooth_labels <- if (!is.null(m_gam$smooth)) {
    vapply(m_gam$smooth, function(s) s$label, character(1))
  } else {
    character(0)
  }
  if (!is.null(.unit) && any(grepl(sprintf("s(%s)", .unit), smooth_labels, fixed = TRUE))) {
    excl <- c(excl, sprintf("s(%s)", .unit))
  }

  # Predict lpmatrix
  X_0 <- mgcv::predict.gam(m_gam, newdata = df_0, type = "lpmatrix", exclude = excl)
  X_1 <- mgcv::predict.gam(m_gam, newdata = df_1, type = "lpmatrix", exclude = excl)
  X_2 <- mgcv::predict.gam(m_gam, newdata = df_2, type = "lpmatrix", exclude = excl)

  beta <- stats::coef(m_gam)
  Vp   <- m_gam$Vp
  if (is.null(Vp)) Vp <- stats::vcov(m_gam)

  fam     <- m_gam$family
  linkinv <- fam$linkinv
  mu.eta  <- fam$mu.eta

  if (is.null(linkinv) || is.null(mu.eta)) {
    stop("The model family does not provide valid `linkinv` or `mu.eta` functions.", call. = FALSE)
  }

  eta_0 <- as.numeric(X_0 %*% beta)
  eta_1 <- as.numeric(X_1 %*% beta)
  eta_2 <- as.numeric(X_2 %*% beta)

  mu_0 <- linkinv(eta_0)
  mu_1 <- linkinv(eta_1)
  mu_2 <- linkinv(eta_2)

  z <- stats::qnorm(1 - (1 - interval) / 2)

  if (scale == "link") {
    D <- (X_2 - X_1) / h_v
    rate <- (eta_2 - eta_1) / h_v
    fitted <- eta_0

    if (uncertainty) {
      fitted_se <- sqrt(pmax(0, rowSums((X_0 %*% Vp) * X_0)))
      fitted_lower <- fitted - z * fitted_se
      fitted_upper <- fitted + z * fitted_se

      var_rate <- pmax(0, rowSums((D %*% Vp) * D))
      rate_se  <- sqrt(var_rate)
      rate_lower <- rate - z * rate_se
      rate_upper <- rate + z * rate_se
    } else {
      fitted_se <- NA_real_
      fitted_lower <- NA_real_
      fitted_upper <- NA_real_
      rate_se <- NA_real_
      rate_lower <- NA_real_
      rate_upper <- NA_real_
    }
  } else {
    # scale == "response"
    rate <- (mu_2 - mu_1) / h_v
    fitted <- mu_0

    if (uncertainty) {
      se_eta_0 <- sqrt(pmax(0, rowSums((X_0 %*% Vp) * X_0)))
      fitted_se <- abs(mu.eta(eta_0)) * se_eta_0
      fitted_lower <- linkinv(eta_0 - z * se_eta_0)
      fitted_upper <- linkinv(eta_0 + z * se_eta_0)

      G <- (mu.eta(eta_2) * X_2 - mu.eta(eta_1) * X_1) / h_v
      var_rate <- pmax(0, rowSums((G %*% Vp) * G))
      rate_se  <- sqrt(var_rate)
      rate_lower <- rate - z * rate_se
      rate_upper <- rate + z * rate_se
    } else {
      fitted_se <- NA_real_
      fitted_lower <- NA_real_
      fitted_upper <- NA_real_
      rate_se <- NA_real_
      rate_lower <- NA_real_
      rate_upper <- NA_real_
    }
  }

  significant_increase <- if (uncertainty) (rate_lower > 0) else NA
  above_threshold      <- (rate > threshold)
  significant_above_threshold <- if (uncertainty) (rate_lower > threshold) else NA

  # Handle zero-epidemic curves explicitly in rate_data
  trt_vec <- as.character(df_0[[.trt]])
  for (zero_trt in zero_epidemic_trts) {
    idx_z <- which(trt_vec == zero_trt)
    if (length(idx_z) > 0) {
      fitted[idx_z] <- 0
      rate[idx_z]   <- 0
      if (uncertainty) {
        fitted_se[idx_z]    <- 0
        fitted_lower[idx_z] <- 0
        fitted_upper[idx_z] <- 0
        rate_se[idx_z]      <- 0
        rate_lower[idx_z]   <- 0
        rate_upper[idx_z]   <- 0
        significant_increase[idx_z] <- FALSE
        significant_above_threshold[idx_z] <- FALSE
      }
      above_threshold[idx_z] <- FALSE
    }
  }

  # Build full rate table
  rate_data <- df_0
  rate_data$fitted <- fitted
  rate_data$fitted_se <- fitted_se
  rate_data$fitted_lower <- fitted_lower
  rate_data$fitted_upper <- fitted_upper
  rate_data$rate <- rate
  rate_data$rate_se <- rate_se
  rate_data$rate_lower <- rate_lower
  rate_data$rate_upper <- rate_upper
  rate_data$significant_increase <- significant_increase
  rate_data$above_threshold <- above_threshold
  rate_data$significant_above_threshold <- significant_above_threshold

  # Warnings for negative rates and extreme uncertainty
  neg_prop <- mean(rate < 0, na.rm = TRUE)
  if (is.finite(neg_prop) && neg_prop > 0.30) {
    warning(sprintf("A large proportion (%.1f%%) of estimated rates are negative. Check model fit or smoothing dimension.",
                    neg_prop * 100), call. = FALSE)
  }

  if (uncertainty) {
    high_se <- which(is.finite(rate_se) & is.finite(rate) & (abs(rate) > 1e-4) & (rate_se / abs(rate) > 50))
    if (length(high_se) > 0) {
      warning("Very high uncertainty detected in some rate estimates (rate_se / |rate| > 50).", call. = FALSE)
    }
  }

  # Compute rate_summary per curve
  group_cols <- c(.trt)
  if (!is.null(.env) && .env %in% names(rate_data) && dplyr::n_distinct(rate_data[[.env]]) > 1) {
    group_cols <- c(group_cols, .env)
  }

  summary_rows <- rate_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(~ {
      sub_t   <- .x[[.time]]
      sub_r   <- .x$rate
      sub_low <- .x$rate_lower
      sub_fit <- .x$fitted
      trt_cur <- as.character(.y[[.trt]][1])

      is_zero_curve <- (trt_cur %in% zero_epidemic_trts) ||
                       all(abs(sub_r) < 1e-6, na.rm = TRUE) ||
                       all(abs(sub_fit) < 1e-6, na.rm = TRUE)

      # Evaluate active growth criterion
      is_growth <- if (growth_criterion == "significant") {
        if (!uncertainty) {
          warning("`growth_criterion = 'significant'` requires `uncertainty = TRUE`; falling back to `rate > threshold`.", call. = FALSE)
          sub_r > threshold
        } else {
          sub_low > threshold
        }
      } else {
        sub_r > threshold
      }

      # Duration and windows
      if (is_zero_curve || !any(is_growth, na.rm = TRUE)) {
        r_max <- if (any(sub_r > 0, na.rm = TRUE) && !is_zero_curve) max(sub_r, na.rm = TRUE) else 0
        if (r_max < 1e-6) r_max <- 0
        t_r_max <- if (r_max > 0) sub_t[which.max(sub_r)] else NA_real_

        span_t <- max(sub_t, na.rm = TRUE) - min(sub_t, na.rm = TRUE)
        r_mean <- if (span_t > 0) trapz_vec_internal(sub_t, sub_r) / span_t else mean(sub_r, na.rm = TRUE)
        r_mean_pos <- 0
        t_growth_start <- NA_real_
        t_growth_end   <- NA_real_
        growth_duration <- 0
        cum_growth <- 0
        n_windows <- 0
      } else {
        r_max <- max(sub_r, na.rm = TRUE)
        if (r_max < 1e-6) {
          r_max <- 0
          t_r_max <- NA_real_
        } else {
          t_r_max <- sub_t[which.max(sub_r)]
        }

        span_t <- max(sub_t, na.rm = TRUE) - min(sub_t, na.rm = TRUE)
        r_mean <- if (span_t > 0) trapz_vec_internal(sub_t, sub_r) / span_t else mean(sub_r, na.rm = TRUE)

        rates_pos <- sub_r[sub_r > threshold]
        r_mean_pos <- if (length(rates_pos) > 0) mean(rates_pos, na.rm = TRUE) else 0

        grow_idx <- which(is_growth)
        t_growth_start <- sub_t[grow_idx[1]]
        t_growth_end   <- sub_t[tail(grow_idx, 1)]

        # Runs of growth
        rle_growth <- rle(is_growth)
        starts <- c(1, cumsum(rle_growth$lengths) + 1)[1:length(rle_growth$lengths)]
        ends   <- cumsum(rle_growth$lengths)

        grow_runs <- which(rle_growth$values)
        n_windows <- length(grow_runs)

        if (n_windows > 0) {
          durations <- vapply(grow_runs, function(idx) {
            s_i <- starts[idx]
            e_i <- ends[idx]
            if (e_i > s_i) sub_t[e_i] - sub_t[s_i] else 0
          }, numeric(1))
          growth_duration <- sum(durations)
        } else {
          growth_duration <- 0
        }

        cum_growth <- trapz_vec_internal(sub_t, pmax(sub_r - threshold, 0))
        if (is.na(cum_growth)) cum_growth <- 0
      }

      tibble::tibble(
        r_max = r_max,
        t_r_max = t_r_max,
        r_mean = r_mean,
        r_mean_positive = r_mean_pos,
        t_growth_start = t_growth_start,
        t_growth_end = t_growth_end,
        growth_duration = growth_duration,
        cumulative_positive_growth = cum_growth,
        n_growth_windows = as.integer(n_windows)
      )
    }) |>
    dplyr::ungroup()

  res <- list(
    data = rate_data,
    summary = summary_rows,
    model = m_gam,
    call = match.call(),
    settings = list(
      scale = scale,
      method = method,
      n_grid = n_grid,
      eps = eps,
      interval = interval,
      threshold = threshold,
      positive_only = positive_only,
      uncertainty = uncertainty,
      growth_criterion = growth_criterion
    ),
    vars = object$vars,
    family = fam$family
  )

  class(res) <- "functional_rate"
  res
}

#' Print a functional_rate object
#'
#' @param x An object of class \code{"functional_rate"}.
#' @param ... Additional arguments.
#' @return The object \code{x}, invisibly.
#' @export
print.functional_rate <- function(x, ...) {
  cat("=========================================================\n")
  cat("Instantaneous Epidemic Rates (functional_rate)\n")
  cat("=========================================================\n")

  .trt  <- x$vars$treatment
  .time <- x$vars$time

  n_curves <- nrow(x$summary)
  t_min <- min(x$data[[.time]], na.rm = TRUE)
  t_max <- max(x$data[[.time]], na.rm = TRUE)

  n_detectable <- sum(x$summary$r_max > x$settings$threshold, na.rm = TRUE)

  cat(sprintf("Number of curves: %d\n", n_curves))
  cat(sprintf("Time domain: [%.2f, %.2f] (%s)\n", t_min, t_max, .time))
  cat(sprintf("Rate scale: %s\n", x$settings$scale))
  cat(sprintf("Method: %s finite difference\n", x$settings$method))
  cat(sprintf("Grid resolution: %d points per curve\n", x$settings$n_grid))
  cat(sprintf("Curves with detectable growth: %d of %d (threshold = %g)\n",
              n_detectable, n_curves, x$settings$threshold))
  cat("=========================================================\n")
  invisible(x)
}

#' Summarize a functional_rate object
#'
#' @param object An object of class \code{"functional_rate"}.
#' @param ... Additional arguments.
#' @return A tibble of curve-level rate phenotypes.
#' @export
summary.functional_rate <- function(object, ...) {
  object$summary
}

#' Plot functional_rate curves and rates
#'
#' @param x An object of class \code{"functional_rate"}.
#' @param which Optional character vector specifying which treatments/groups to display.
#' @param type Character specifying what to plot: \code{"both"} (default; fitted curve
#'   in top panel and derivative in bottom panel), \code{"curve"} (only trajectory),
#'   or \code{"rate"} (only instantaneous rate).
#' @param show_ci Logical; whether to show confidence ribbons. Default is \code{TRUE}.
#' @param show_threshold Logical; whether to show threshold and zero horizontal lines.
#'   Default is \code{TRUE}.
#' @param ncol Optional integer specifying number of facet columns.
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} or combined \code{patchwork} object.
#' @export
plot.functional_rate <- function(
    x,
    which = NULL,
    type = c("both", "curve", "rate"),
    show_ci = TRUE,
    show_threshold = TRUE,
    ncol = NULL,
    ...
) {
  type <- match.arg(type)

  .trt  <- x$vars$treatment
  .time <- x$vars$time

  df <- x$data

  if (!is.null(which)) {
    df <- df |> dplyr::filter(.data[[.trt]] %in% which)
    if (nrow(df) == 0) {
      stop("None of the specified treatments found in `x`.", call. = FALSE)
    }
  }

  y_label_curve <- if (x$settings$scale == "response") "Disease intensity, S(t)" else "Linear predictor, \u03b7(t)"
  y_label_rate  <- if (x$settings$scale == "response") "Instantaneous rate, S'(t)" else "Rate of linear change, \u03b7'(t)"

  # Plot curve S(t)
  p_curve <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[.time]], y = .data$fitted, color = .data[[.trt]], fill = .data[[.trt]]))

  if (show_ci && !all(is.na(df$fitted_lower))) {
    p_curve <- p_curve +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$fitted_lower, ymax = .data$fitted_upper),
                           alpha = 0.2, color = NA)
  }

  p_curve <- p_curve +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = .time, y = y_label_curve, color = .trt, fill = .trt) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  if (type == "curve") return(p_curve)

  # Plot rate S'(t)
  p_rate <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[.time]], y = .data$rate, color = .data[[.trt]], fill = .data[[.trt]]))

  if (show_threshold) {
    p_rate <- p_rate + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
    if (x$settings$threshold != 0) {
      p_rate <- p_rate + ggplot2::geom_hline(yintercept = x$settings$threshold,
                                             linetype = "dotted", color = "firebrick")
    }
  }

  if (show_ci && !all(is.na(df$rate_lower))) {
    p_rate <- p_rate +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$rate_lower, ymax = .data$rate_upper),
                           alpha = 0.2, color = NA)
  }

  p_rate <- p_rate +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = .time, y = y_label_rate, color = .trt, fill = .trt) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  if (type == "rate") return(p_rate)

  # Combine both
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p_curve / p_rate + patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "top")
    return(combined)
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::plot_grid(p_curve, p_rate, ncol = 1, align = "v"))
  } else {
    message("Package 'patchwork' or 'cowplot' recommended to stack plots. Returning list.")
    return(list(curve = p_curve, rate = p_rate))
  }
}

#' Augment generic for functional objects
#'
#' @param x An object to augment.
#' @param ... Additional arguments.
#' @return An augmented data frame or tibble.
#' @export
augment <- function(x, ...) {
  UseMethod("augment")
}

#' Augment functional_rate
#'
#' @param x An object of class \code{"functional_rate"}.
#' @param ... Additional arguments.
#' @return A tibble of the full rate evaluation grid with predictions, derivatives,
#'   and confidence intervals.
#' @export
augment.functional_rate <- function(x, ...) {
  if (!inherits(x, "functional_rate")) {
    stop("`x` must be of class 'functional_rate'.", call. = FALSE)
  }
  x$data
}

#' Augment functional_rate (alias)
#'
#' @param x An object of class \code{"functional_rate"}.
#' @param ... Additional arguments.
#' @return A tibble of the full rate evaluation grid.
#' @export
augment_functional_rate <- function(x, ...) {
  augment.functional_rate(x, ...)
}
