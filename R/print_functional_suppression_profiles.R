#' Print functional suppression profiles
#'
#' @param x An object of class \code{"functional_suppression_profiles"}.
#' @param ... Additional arguments.
#'
#' @export
print.functional_suppression_profiles <- function(x, ...) {
  n_treatments <- length(unique(x$clusters$treatment))
  n_profiles <- length(unique(x$clusters$cluster))
  
  cat("Functional Suppression Profiles\n")
  cat("--------------------------------\n")
  cat("Reference treatment:", x$reference, "\n")
  cat("Treatments:", n_treatments + 1, "\n") # +1 for reference
  cat("Profiles:", n_profiles, "\n")
  cat("Distance: functional distance among smoothed DSP curves\n")
  cat("Persistence threshold:", x$parameters$threshold, "\n\n")
  
  # Top treatments
  cat("Top treatments by mean rank:\n")
  top_n <- min(3, nrow(x$classification))
  top_df <- x$classification[1:top_n, ]
  
  trt_var <- x$parameters$treatment
  for (i in 1:top_n) {
    cat(sprintf("  %-3s %-12s mean rank = %.1f\n", 
                as.character(top_df[[trt_var]][i]), 
                as.character(top_df$profile[i]), 
                top_df$mean_rank[i]))
  }
  
  invisible(x)
}
