#' Suggest GAM smoothing parameters for epidemic curve models
#'
#' @description
#' A helper function to guide selection of the \code{k_smooth}, \code{k_trt},
#' \code{k_env}, and \code{gamma} parameters used by \code{\link{functional_curves}}.
#' Can infer the effective replication from a data frame or accept scalar values
#' directly when the data are not yet available.
#'
#' For sparse epidemic data it is strongly recommended to use
#' \code{rule = "minimum"} so that k values are based on the least-replicated
#' treatment-by-environment combination, guarding against over-fitting.
#'
#' @param data Optional data frame.  If supplied, \code{time}, \code{treatment},
#'   and optionally \code{environment} are used to compute summaries of the
#'   number of unique time points and environments per treatment.
#' @param time Unquoted column name for the time variable, or a character
#'   string naming the column.
#' @param treatment Unquoted column name for the treatment / cultivar variable,
#'   or a character string naming the column.
#' @param environment Optional unquoted column name for the environment
#'   variable, or a character string naming the column.
#' @param n_time Integer.  Number of unique time points to use directly
#'   (ignored when \code{data} is supplied).
#' @param n_env Integer.  Number of unique environments to use directly
#'   (ignored when \code{data} is supplied, or when there is no environment
#'   variable).
#' @param rule Character string. How to summarise the distribution of
#'   replication counts across treatment-by-environment combinations:
#'   \code{"minimum"} (default, conservative) or \code{"median"}.
#' @param smoothness Character string controlling how liberal the
#'   recommendations are: \code{"conservative"} (default), \code{"moderate"},
#'   or \code{"flexible"}.
#'
#' @details
#' The function computes, for every treatment-by-environment combination, the
#' number of unique time points at which observations are available.  It then
#' summarises these counts using the chosen \code{rule} to obtain
#' \code{effective_n_time}.  Similarly it computes, per treatment, the number
#' of unique environments, and summarises using the same \code{rule} to obtain
#' \code{effective_n_env}.
#'
#' Recommended k values follow these heuristics:
#' \itemize{
#'   \item \code{k_smooth}: global smooth over time; capped below
#'     \code{effective_n_time - 1}.
#'   \item \code{k_trt}: treatment-specific smooth; capped below
#'     \code{k_smooth}.
#'   \item \code{k_env}: environment random effect; capped below
#'     \code{effective_n_env}.
#'   \item \code{gamma}: penalty multiplier; higher values encourage
#'     smoother fits.
#' }
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{time_summary}}{Named numeric vector with minimum, median, and
#'     maximum unique time points per treatment-by-environment combination.}
#'   \item{\code{environment_summary}}{Named numeric vector with minimum, median,
#'     and maximum unique environments per treatment (or \code{NULL} when no
#'     environment variable is given).}
#'   \item{\code{effective_n_time}}{The effective number of time points chosen
#'     by \code{rule}.}
#'   \item{\code{effective_n_env}}{The effective number of environments chosen
#'     by \code{rule}, or \code{NULL}.}
#'   \item{\code{k_smooth}}{Recommended basis dimension for the global time smooth.}
#'   \item{\code{k_trt}}{Recommended basis dimension for the treatment-specific smooth.}
#'   \item{\code{k_env}}{Recommended basis dimension for the environment random effect.}
#'   \item{\code{gamma}}{Recommended penalisation multiplier.}
#'   \item{\code{message}}{A short interpretation message.}
#' }
#'
#' @examples
#' # Using explicit values
#' suggest_k(n_time = 5, n_env = 8)
#'
#' # Inferring from a data frame (unquoted column names)
#' df <- data.frame(
#'   time      = rep(1:6, times = 6),
#'   cultivar  = rep(c("A", "B", "C"), each = 12),
#'   env       = rep(rep(c("E1", "E2"), each = 6), times = 3),
#'   severity  = runif(36, 0, 0.5)
#' )
#' suggest_k(
#'   data        = df,
#'   time        = time,
#'   treatment   = cultivar,
#'   environment = env,
#'   rule        = "minimum",
#'   smoothness  = "conservative"
#' )
#'
#' @export
suggest_k <- function(
    data        = NULL,
    time        = NULL,
    treatment   = NULL,
    environment = NULL,
    n_time      = NULL,
    n_env       = NULL,
    rule        = c("minimum", "median"),
    smoothness  = c("conservative", "moderate", "flexible")
) {
  rule       <- match.arg(rule)
  smoothness <- match.arg(smoothness)

  # ---- resolve tidy-eval column names ----------------------------------------
  time_col <- NULL
  trt_col  <- NULL
  env_col  <- NULL

  if (!is.null(data)) {
    time_quo <- rlang::enquo(time)
    trt_quo  <- rlang::enquo(treatment)
    env_quo  <- rlang::enquo(environment)

    time_col <- if (!rlang::quo_is_null(time_quo)) rlang::as_name(time_quo) else NULL
    trt_col  <- if (!rlang::quo_is_null(trt_quo))  rlang::as_name(trt_quo)  else NULL
    env_col  <- if (!rlang::quo_is_null(env_quo))  rlang::as_name(env_quo)  else NULL

    if (is.null(time_col) || is.null(trt_col)) {
      stop("When `data` is supplied, both `time` and `treatment` must be specified.", call. = FALSE)
    }
    bad <- setdiff(c(time_col, trt_col, env_col), names(data))
    if (length(bad)) stop("Column(s) not found in `data`: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # ---- compute summaries from data -------------------------------------------
  time_summary <- NULL
  env_summary  <- NULL

  if (!is.null(data)) {
    if (!is.null(env_col)) {
      time_counts <- tapply(
        data[[time_col]],
        list(data[[trt_col]], data[[env_col]]),
        function(x) length(unique(x[!is.na(x)]))
      )
      time_vec <- as.numeric(time_counts)
      time_vec <- time_vec[!is.na(time_vec)]

      env_counts <- tapply(
        data[[env_col]],
        data[[trt_col]],
        function(x) length(unique(x[!is.na(x)]))
      )
      env_vec <- as.numeric(env_counts)
    } else {
      time_counts <- tapply(
        data[[time_col]],
        data[[trt_col]],
        function(x) length(unique(x[!is.na(x)]))
      )
      time_vec <- as.numeric(time_counts)
      time_vec <- time_vec[!is.na(time_vec)]
      env_vec  <- NULL
    }

    time_summary <- c(
      minimum = min(time_vec),
      median  = stats::median(time_vec),
      maximum = max(time_vec)
    )

    if (!is.null(env_vec)) {
      env_summary <- c(
        minimum = min(env_vec),
        median  = stats::median(env_vec),
        maximum = max(env_vec)
      )
    }

    n_time <- if (rule == "minimum") time_summary["minimum"] else time_summary["median"]
    n_env  <- if (!is.null(env_summary)) {
      if (rule == "minimum") env_summary["minimum"] else env_summary["median"]
    } else NULL

  } else {
    # use explicit scalars
    if (is.null(n_time)) stop("`n_time` must be provided when `data` is NULL.", call. = FALSE)
  }

  effective_n_time <- as.integer(n_time)
  effective_n_env  <- if (!is.null(n_env)) as.integer(n_env) else NULL

  # ---- recommend k values -----------------------------------------------------
  gamma_tbl <- c(conservative = 1.8, moderate = 1.6, flexible = 1.4)
  gamma_rec <- gamma_tbl[[smoothness]]

  # k_smooth: global smooth over time
  k_smooth_raw <- if (effective_n_time <= 4) 4L
                  else if (effective_n_time <= 6) 5L
                  else if (effective_n_time <= 8) 7L
                  else 10L
  if (smoothness == "flexible")     k_smooth_raw <- min(k_smooth_raw + 2L, effective_n_time - 1L)
  else if (smoothness == "moderate") k_smooth_raw <- min(k_smooth_raw + 1L, effective_n_time - 1L)
  k_smooth_rec <- min(k_smooth_raw, max(3L, effective_n_time - 1L))

  # k_trt: treatment-specific smooth
  k_trt_raw <- if (effective_n_time <= 4) 3L
               else if (effective_n_time <= 6) 4L
               else if (smoothness == "conservative") 5L
               else if (smoothness == "moderate") 5L
               else 6L
  k_trt_rec <- min(k_trt_raw, k_smooth_rec - 1L)
  k_trt_rec <- max(k_trt_rec, 2L)

  # k_env: environment random effect basis dimension
  k_env_rec <- NULL
  if (!is.null(effective_n_env)) {
    k_env_raw <- if (effective_n_env < 6) 3L
                 else if (smoothness == "conservative") 3L
                 else 4L
    k_env_rec <- min(k_env_raw, max(2L, effective_n_env - 1L))
  }

  # ---- interpretation message -------------------------------------------------
  parts <- c(
    sprintf("effective_n_time = %d (%s across combinations)", effective_n_time, rule),
    if (!is.null(effective_n_env)) sprintf("effective_n_env = %d (%s across treatments)", effective_n_env, rule),
    sprintf("smoothness = '%s'", smoothness),
    sprintf("k_smooth = %d, k_trt = %d, k_env = %s, gamma = %.1f",
            k_smooth_rec, k_trt_rec,
            if (!is.null(k_env_rec)) as.character(k_env_rec) else "NA (no environment)",
            gamma_rec)
  )
  msg <- paste(parts, collapse = "\n  ")

  if (effective_n_time <= 4) {
    msg <- paste0(msg, "\n  Warning: very few time points detected. Prefer rule = 'minimum' and smoothness = 'conservative'.")
  }

  # ---- return -----------------------------------------------------------------
  structure(
    list(
      time_summary     = time_summary,
      environment_summary = env_summary,
      effective_n_time = effective_n_time,
      effective_n_env  = effective_n_env,
      k_smooth         = k_smooth_rec,
      k_trt            = k_trt_rec,
      k_env            = k_env_rec,
      gamma            = gamma_rec,
      message          = msg
    ),
    class = "suggest_k"
  )
}

#' Print method for suggest_k
#'
#' @param x An object of class \code{"suggest_k"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.suggest_k <- function(x, ...) {
  cat("── suggest_k: GAM parameter recommendations ──────────────────────\n")
  cat(" ", x$message, "\n")
  cat("\nSuggested call to functional_curves():\n")
  cat(sprintf(
    "  functional_curves(..., k_smooth = %d, k_trt = %d, k_env = %s, gamma = %.1f)\n",
    x$k_smooth, x$k_trt,
    if (!is.null(x$k_env)) as.character(x$k_env) else "NA",
    x$gamma
  ))
  if (!is.null(x$time_summary)) {
    cat("\nTime points per treatment-by-environment:\n")
    cat(sprintf("  min = %g, median = %g, max = %g\n",
                x$time_summary["minimum"], x$time_summary["median"], x$time_summary["maximum"]))
  }
  if (!is.null(x$environment_summary)) {
    cat("Environments per treatment:\n")
    cat(sprintf("  min = %g, median = %g, max = %g\n",
                x$environment_summary["minimum"], x$environment_summary["median"],
                x$environment_summary["maximum"]))
  }
  invisible(x)
}
