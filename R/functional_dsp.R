#' Functional Contrast (Disease Suppression Profiles)
#'
#' @description
#' Calculates functional contrast curves (Disease Suppression Profiles - DSP) between a 
#' reference treatment (e.g., untreated control) and other fungicide treatments over time. 
#' The DSP is defined by default as: \deqn{DSP(t) = y_{control}(t) - y_{treatment}(t)}
#' where \eqn{y(t)} is the disease intensity at time \eqn{t}. The goal is to evaluate
#' the intensity, timing, and persistence of the fungicide protection.
#'
#' @param data A data frame containing disease progress observations, or an object of class \code{"functional_curves"}.
#' @param time Character string naming the time variable.
#' @param response Character string naming the response variable.
#' @param treatment Character string naming the treatment variable.
#' @param reference Character string specifying the reference treatment level. Default is \code{"Control"}.
#' @param group Optional character vector of grouping variables (e.g., \code{c("experiment", "location")}).
#' @param smooth Logical; if \code{TRUE}, the function could apply smoothing (currently placeholder).
#' @param grid Optional numeric vector of common time points for interpolation. If provided, curves are interpolated.
#' @param contrast Character specifying the contrast type. Currently \code{"difference"} is implemented.
#' @param keep_reference Logical; if \code{TRUE}, includes the reference treatment in the output (with DSP = 0). Default is \code{FALSE}.
#'
#' @return An object of class \code{"functional_dsp"} (inheriting from \code{"tbl_df"}), containing:
#' \itemize{
#'   \item grouping variables (if any)
#'   \item treatment variable
#'   \item time variable
#'   \item \code{y_reference}: disease intensity of the reference
#'   \item \code{y_treatment}: disease intensity of the treatment
#'   \item \code{DSP}: the computed contrast
#'   \item \code{contrast_type}: type of contrast used
#'   \item \code{reference}: name of the reference treatment
#' }
#'
#' @seealso \code{\link{functional_summary}}, \code{\link{plot_dsp}}
#'
#' @export
functional_contrast <- function(
    data,
    time = "time",
    response = "response",
    treatment = "treatment",
    reference = "Control",
    group = NULL,
    smooth = FALSE,
    grid = NULL,
    contrast = c("difference", "relative"),
    keep_reference = FALSE
) {
  contrast <- match.arg(contrast)
  
  if (inherits(data, "functional_curves")) {
    time_col <- data$vars$time
    trt_col <- data$vars$treatment
    resp_col <- "mu"
    df <- data$curves
    group_cols <- NULL
  } else {
    time_col <- time
    trt_col <- treatment
    resp_col <- response
    df <- tibble::as_tibble(data)
    group_cols <- group
  }
  
  if (!time_col %in% names(df)) stop(sprintf("Time column '%s' not found.", time_col))
  if (!trt_col %in% names(df)) stop(sprintf("Treatment column '%s' not found.", trt_col))
  if (!resp_col %in% names(df)) stop(sprintf("Response column '%s' not found.", resp_col))
  if (!is.null(group_cols)) {
    missing_groups <- setdiff(group_cols, names(df))
    if (length(missing_groups) > 0) {
      stop("Grouping columns not found: ", paste(missing_groups, collapse = ", "))
    }
  }
  
  # Ensure reference exists
  if (!reference %in% df[[trt_col]]) {
    stop(sprintf("Reference level '%s' not found in treatment column.", reference))
  }
  
  # Interpolation if grid is provided
  if (!is.null(grid)) {
    grid <- sort(unique(grid))
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, trt_col)))) |>
      dplyr::reframe(
        !!time_col := grid,
        !!resp_col := stats::approx(x = .data[[time_col]], y = .data[[resp_col]], xout = grid, rule = 2)$y
      )
  } else {
    # If no grid, ensure same time points by interpolating treatments to match reference time points if needed
    # But usually data is already on a grid if from functional_curves.
    # For raw data, we will just use approx against the union of time points per group.
    
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::group_modify(~ {
        common_times <- sort(unique(.x[[time_col]]))
        
        .x |>
          dplyr::group_by(.data[[trt_col]]) |>
          dplyr::reframe(
            !!time_col := common_times,
            !!resp_col := stats::approx(x = .data[[time_col]], y = .data[[resp_col]], xout = common_times, rule = 2)$y
          )
      }) |>
      dplyr::ungroup()
  }
  
  # Split reference
  ref_df <- df |>
    dplyr::filter(.data[[trt_col]] == reference) |>
    dplyr::select(dplyr::all_of(c(group_cols, time_col)), y_reference = dplyr::all_of(resp_col))
  
  # Check if reference exists for all groups
  if (!is.null(group_cols)) {
    ref_counts <- ref_df |> dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |> dplyr::tally()
    all_groups <- df |> dplyr::distinct(dplyr::across(dplyr::all_of(group_cols)))
    if (nrow(ref_counts) < nrow(all_groups)) {
      stop("Reference treatment is missing in some groups.")
    }
  }
  
  # Calculate contrast
  res <- df |>
    dplyr::left_join(ref_df, by = c(group_cols, time_col)) |>
    dplyr::rename(y_treatment = dplyr::all_of(resp_col)) |>
    dplyr::mutate(
      DSP = if (contrast == "difference") {
        .data$y_reference - .data$y_treatment
      } else {
        (.data$y_reference - .data$y_treatment) / dplyr::na_if(.data$y_reference, 0)
      },
      contrast_type = contrast,
      reference = reference
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(group_cols, trt_col, time_col))))
  
  if (!keep_reference) {
    res <- res |> dplyr::filter(.data[[trt_col]] != reference)
  }
  
  # Create attributes for downstream methods
  attr(res, "dsp_vars") <- list(
    time = time_col,
    treatment = trt_col,
    response = "DSP",
    group = group_cols
  )
  class(res) <- c("functional_dsp", class(res))
  
  res
}

#' @keywords internal
trapz_vec_internal <- function(x, y) {
  ok <- !(is.na(x) | is.na(y))
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

#' Summarize Disease Suppression Profiles
#'
#' @description
#' Extracts functional descriptors from Disease Suppression Profiles (DSPs).
#'
#' @param object An object of class \code{"functional_dsp"} or a data frame with a DSP column.
#' @param time Character string naming the time column.
#' @param response Character string naming the DSP column.
#' @param treatment Character string naming the treatment column.
#' @param group Optional character vector of grouping variables.
#' @param threshold Numeric; minimum suppression threshold for calculating persistence. Default is 0.
#' @param positive_only Logical; if \code{TRUE}, negative DSP values are set to 0 before area/centroid calculations. Default is \code{TRUE}.
#'
#' @return A \code{tibble} with one row per group/treatment containing:
#' \itemize{
#'   \item \code{protected_area}: Trapezoidal integral of DSP(t).
#'   \item \code{mean_suppression}: Temporal mean of DSP(t).
#'   \item \code{max_suppression}: Maximum value of DSP(t).
#'   \item \code{time_max_suppression}: Time at which maximum suppression occurs.
#'   \item \code{persistence}: Total duration where DSP(t) > threshold.
#'   \item \code{centroid}: Temporal centroid defined as \eqn{\int t \cdot DSP(t) dt / \int DSP(t) dt}.
#'   \item \code{energy}: Integral of \eqn{DSP(t)^2}.
#'   \item \code{early_area, mid_area, late_area}: Protected area in the first, second, and third tertiles of the time period.
#'   \item \code{late_decline_slope}: Slope of DSP(t) in the late period.
#'   \item \code{n_time_points}: Number of time points evaluated.
#' }
#' @examples
#' \dontrun{
#' # Using a generic data frame with custom column names
#' dat <- data.frame(
#'   DAA = c(0, 10, 20),
#'   fungicide = c("A", "A", "A"),
#'   contrast = c(0, 10, 0)
#' )
#' functional_summary(
#'   dat, 
#'   time = "DAA", 
#'   response = "contrast", 
#'   treatment = "fungicide"
#' )
#' }
#'
#' @export
functional_summary <- function(
    object,
    time = "time",
    response = "DSP",
    treatment = "treatment",
    group = NULL,
    threshold = 0,
    positive_only = TRUE
) {
  df <- tibble::as_tibble(object)
  
  if (inherits(object, "functional_dsp")) {
    vars <- attr(object, "dsp_vars")
    if (!is.null(vars)) {
      if (missing(time) && !is.null(vars$time)) time <- vars$time
      if (missing(response) && !is.null(vars$response)) response <- vars$response
      if (missing(treatment) && !is.null(vars$treatment)) treatment <- vars$treatment
      if (missing(group) && !is.null(vars$group)) group <- vars$group
    }
  }
  
  if (!all(c(time, response, treatment) %in% names(df))) {
    stop("Missing required columns in object.")
  }
  
  res <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group, treatment)))) |>
    dplyr::summarise(
      n_time_points = dplyr::n(),
      protected_area = {
        t <- .data[[time]]
        y <- .data[[response]]
        y_val <- if (positive_only) pmax(y, 0) else y
        trapz_vec_internal(t, y_val)
      },
      mean_suppression = {
        y <- .data[[response]]
        mean(y, na.rm = TRUE)
      },
      max_suppression = {
        y <- .data[[response]]
        max(y, na.rm = TRUE)
      },
      time_max_suppression = {
        t <- .data[[time]]
        y <- .data[[response]]
        idx <- which.max(y)
        if (length(idx) > 0) t[idx[1]] else NA_real_
      },
      persistence = {
        t <- .data[[time]]
        y <- .data[[response]]
        above <- y > threshold
        if (sum(above) < 2) {
          0
        } else {
          t_above <- t[above]
          max(t_above) - min(t_above)
        }
      },
      centroid = {
        t <- .data[[time]]
        y <- .data[[response]]
        y_val <- if (positive_only) pmax(y, 0) else y
        num <- trapz_vec_internal(t, t * y_val)
        den <- trapz_vec_internal(t, y_val)
        if (is.na(den) || den == 0) NA_real_ else num / den
      },
      energy = {
        t <- .data[[time]]
        y <- .data[[response]]
        trapz_vec_internal(t, y^2)
      },
      early_area = {
        t <- .data[[time]]
        y <- .data[[response]]
        y_val <- if (positive_only) pmax(y, 0) else y
        cuts <- stats::quantile(t, probs = c(1/3, 2/3), na.rm = TRUE)
        idx <- t <= cuts[1]
        trapz_vec_internal(t[idx], y_val[idx])
      },
      mid_area = {
        t <- .data[[time]]
        y <- .data[[response]]
        y_val <- if (positive_only) pmax(y, 0) else y
        cuts <- stats::quantile(t, probs = c(1/3, 2/3), na.rm = TRUE)
        idx <- t > cuts[1] & t <= cuts[2]
        trapz_vec_internal(t[idx], y_val[idx])
      },
      late_area = {
        t <- .data[[time]]
        y <- .data[[response]]
        y_val <- if (positive_only) pmax(y, 0) else y
        cuts <- stats::quantile(t, probs = c(1/3, 2/3), na.rm = TRUE)
        idx <- t > cuts[2]
        trapz_vec_internal(t[idx], y_val[idx])
      },
      late_decline_slope = {
        t <- .data[[time]]
        y <- .data[[response]]
        cuts <- stats::quantile(t, probs = c(1/3, 2/3), na.rm = TRUE)
        idx <- t > cuts[2]
        if (sum(idx) >= 2) {
          stats::coef(stats::lm(y[idx] ~ t[idx]))[2]
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    )
    
  if (inherits(object, "functional_dsp")) {
    class(res) <- unique(c("functional_dsp", class(res)))
    attr(res, "dsp_vars") <- attr(object, "dsp_vars")
  }
  
  res
}

#' Plot Disease Suppression Profiles
#'
#' @description
#' Plots functional contrast curves (DSP) over time.
#'
#' @param object An object of class \code{"functional_dsp"}.
#' @param time Character string for the time column (auto-detected if object is \code{"functional_dsp"}).
#' @param response Character string for the DSP column.
#' @param treatment Character string for the treatment column.
#' @param group Optional character vector of grouping variables.
#' @param facet Optional character string specifying a variable for faceting.
#' @param show_zero Logical; whether to draw a horizontal reference line at y = 0. Default is \code{TRUE}.
#' @param linewidth Numeric; width of the lines.
#' @param alpha Numeric; transparency of the lines.
#' @param ylab Character string; label for the y-axis.
#' @param xlab Character string; label for the x-axis.
#'
#' @return A \code{ggplot} object.
#' @export
plot_dsp <- function(
    object,
    time = NULL,
    response = NULL,
    treatment = NULL,
    group = NULL,
    facet = NULL,
    show_zero = TRUE,
    linewidth = 1,
    alpha = 1,
    ylab = "Disease suppression profile",
    xlab = "Time"
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  
  if (inherits(object, "functional_dsp")) {
    vars <- attr(object, "dsp_vars")
    if (is.null(time)) time <- vars$time
    if (is.null(response)) response <- vars$response
    if (is.null(treatment)) treatment <- vars$treatment
    if (is.null(group) && length(vars$group) > 0) group <- vars$group
  } else {
    if (is.null(time)) time <- "time"
    if (is.null(response)) response <- "DSP"
    if (is.null(treatment)) treatment <- "treatment"
  }
  
  df <- tibble::as_tibble(object)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[time]], y = .data[[response]], color = .data[[treatment]], group = .data[[treatment]]))
  
  if (show_zero) {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  }
  
  p <- p + 
    ggplot2::geom_line(linewidth = linewidth, alpha = alpha) +
    ggplot2::labs(x = xlab, y = ylab, color = "Treatment") +
    ggplot2::theme_classic()
  
  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  } else if (!is.null(group) && length(group) == 1) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", group[1])))
  }
  
  p
}

#' Rank Treatments based on Functional DSP Descriptors
#'
#' @description
#' Generates rankings of treatments based on extracted DSP functional descriptors.
#'
#' @param summary_data The tibble output from \code{\link{functional_summary}}.
#' @param metrics Character vector of metrics to rank.
#' @param higher_is_better Logical vector or named logical vector indicating if higher values are better for each metric. If a single logical, applies to all.
#' @param group Optional character vector of grouping variables.
#' @param ties.method Character string specifying how to handle ties. Default is \code{"average"}.
#'
#' @return A tibble with ranks for each specified metric, and an \code{average_rank}.
#' @export
rank_dsp <- function(
    summary_data,
    metrics = c("protected_area", "max_suppression", "persistence", "centroid"),
    higher_is_better = TRUE,
    group = NULL,
    ties.method = "average"
) {
  df <- tibble::as_tibble(summary_data)
  
  if (length(higher_is_better) == 1) {
    hib <- rep(higher_is_better, length(metrics))
    names(hib) <- metrics
  } else {
    if (is.null(names(higher_is_better))) {
      if (length(higher_is_better) != length(metrics)) stop("higher_is_better must have the same length as metrics if not named.")
      names(higher_is_better) <- metrics
      hib <- higher_is_better
    } else {
      hib <- higher_is_better
    }
  }
  
  if (is.null(group)) {
    grouped_df <- df |> dplyr::group_by()
  } else {
    grouped_df <- df |> dplyr::group_by(dplyr::across(dplyr::all_of(group)))
  }
  
  for (m in metrics) {
    if (!m %in% names(df)) stop(sprintf("Metric '%s' not found in summary_data.", m))
    
    grouped_df <- grouped_df |>
      dplyr::mutate(
        !!paste0("rank_", m) := {
          vals <- .data[[m]]
          if (hib[m]) {
            rank(-vals, ties.method = ties.method, na.last = "keep")
          } else {
            rank(vals, ties.method = ties.method, na.last = "keep")
          }
        }
      )
  }
  
  res <- grouped_df |> dplyr::ungroup()
  
  rank_cols <- paste0("rank_", metrics)
  res <- res |>
    dplyr::mutate(average_rank = rowMeans(dplyr::select(res, dplyr::all_of(rank_cols)), na.rm = TRUE))
  
  res
}

#' Plot DSP Rankings Heatmap
#'
#' @description
#' Creates a heatmap visualization of treatment rankings based on DSP descriptors.
#'
#' @param ranks_data Tibble output from \code{\link{rank_dsp}}.
#' @param treatment Character string naming the treatment column (for the y-axis).
#' @param group Optional character string for faceting by group.
#'
#' @return A \code{ggplot} object.
#' @export
plot_dsp_rank_heatmap <- function(
    ranks_data,
    treatment = "treatment",
    group = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  
  rank_cols <- grep("^rank_", names(ranks_data), value = TRUE)
  if (length(rank_cols) == 0) stop("No rank columns found. Did you run rank_dsp()?")
  
  long_df <- ranks_data |>
    dplyr::select(dplyr::all_of(c(treatment, group, rank_cols, "average_rank"))) |>
    tidyr::pivot_longer(
      cols = c(dplyr::all_of(rank_cols), "average_rank"),
      names_to = "metric",
      values_to = "rank"
    ) |>
    dplyr::mutate(
      metric = gsub("^rank_", "", .data$metric),
      metric = factor(.data$metric, levels = c(gsub("^rank_", "", rank_cols), "average_rank"))
    )
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = .data$metric, y = .data[[treatment]], fill = .data$rank)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(.data$rank, 1)), color = "black") +
    ggplot2::scale_fill_gradient2(low = "forestgreen", mid = "white", high = "red", midpoint = mean(long_df$rank, na.rm = TRUE)) +
    ggplot2::labs(x = "Metric", y = "Treatment", fill = "Rank") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  if (!is.null(group) && length(group) == 1) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", group[1])))
  }
  
  p
}

#' #' Simulate Disease Suppression Profile Data
#'
#' @description
#' Simulates a dataset for demonstrating functional analysis with Disease Suppression Profiles.
#' Returns a data frame with varying treatments and replicates over multiple time points.
#'
#' @return A \code{tibble} with columns \code{treatment}, \code{rep}, \code{time}, \code{y_control}, \code{dsp_true}, \code{severity_true}, and \code{severity}.
#' @export
simulate_dsp_data <- function() {
  set.seed(123)
  
  time_seq <- seq(0, 60, by = 5)
  
  control_curve <- function(t) {
    78 / (1 + exp(-0.13 * (t - 35)))
  }
  
  dsp_fun <- function(t, profile) {
    dplyr::case_when(
      profile == "Early suppression" ~ 
        156.7 * exp(-((t - 30)^2) / (2 * 12^2)),
      
      profile == "Late suppression" ~ 
        78.9 * exp(-((t - 58)^2) / (2 * 15^2)),
      
      profile == "Persistent suppression" ~ 
        43.9 * stats::plogis((t - 20) / 3) * (1 - stats::plogis((t - 65) / 5)),
      
      profile == "Control" ~ 0
    )
  }
  
  treatments <- tibble::tibble(
    treatment = c(
      "Control",
      "Early suppression",
      "Persistent suppression",
      "Late suppression"
    ),
    profile = c(
      "Control",
      "Early suppression",
      "Persistent suppression",
      "Late suppression"
    )
  )
  
  dat_concept <- tidyr::expand_grid(
    treatment = treatments$treatment,
    rep = 1:4,
    time = time_seq
  ) |>
    dplyr::left_join(treatments, by = "treatment") |>
    dplyr::mutate(
      y_control = control_curve(.data$time),
      dsp_true = dsp_fun(.data$time, .data$profile),
      severity_true = pmax(0, .data$y_control - .data$dsp_true)
    ) |>
    dplyr::group_by(.data$treatment, .data$rep) |>
    dplyr::mutate(
      severity_true = cummax(.data$severity_true),
      severity = pmin(100, pmax(0, .data$severity_true + stats::rnorm(dplyr::n(), 0, 2)))
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.data$profile)
    
  return(dat_concept)
}
