#' Summarize functional suppression profiles
#'
#' @param object An object of class \code{"functional_suppression_profiles"}.
#' @param ... Additional arguments.
#'
#' @return A list of class \code{"summary.functional_suppression_profiles"}.
#'
#' @export
summary.functional_suppression_profiles <- function(object, ...) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  
  metrics <- object$parameters$metrics
  available_metrics <- intersect(metrics, names(object$summary))
  
  class_df <- object$classification
  trt_var <- object$parameters$treatment
  
  profile_summary <- class_df |>
    dplyr::group_by(.data$profile) |>
    dplyr::summarise(
      n_treatments = dplyr::n(),
      treatments = paste(.data[[trt_var]], collapse = ", "),
      mean_rank = mean(.data$mean_rank, na.rm = TRUE),
      best_treatment = .data[[trt_var]][which.min(.data$mean_rank)],
      .groups = "drop"
    )
    
  res <- list(
    classification = class_df,
    profile_summary = profile_summary,
    ranking = object$ranking,
    metrics = object$summary,
    silhouette = object$silhouette_avg
  )
  
  class(res) <- "summary.functional_suppression_profiles"
  res
}

#' Print summary of functional suppression profiles
#'
#' @param x An object of class \code{"summary.functional_suppression_profiles"}.
#' @param ... Additional arguments.
#'
#' @export
print.summary.functional_suppression_profiles <- function(x, ...) {
  cat("----------------------------------------\n")
  cat("Functional Suppression Profiles\n")
  cat("----------------------------------------\n\n")
  
  cat("Number of treatments:", nrow(x$classification) + 1, "\n") # +1 for reference
  cat("Estimated functional profiles:", nrow(x$profile_summary), "\n")
  
  if (!is.null(x$silhouette) && !is.na(x$silhouette)) {
    cat("Average silhouette width:", round(x$silhouette, 2), "\n")
  }
  
  cat("\n")
  
  # Print treatments for each profile
  for (i in 1:nrow(x$profile_summary)) {
    row <- x$profile_summary[i, ]
    cat(as.character(row$profile), "|", row$n_treatments, "|", row$treatments, "|", round(row$mean_rank, 1), "|", as.character(row$best_treatment), "\n")
  }
  
  invisible(x)
}
