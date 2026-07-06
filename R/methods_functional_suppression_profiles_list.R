#' Print a list of functional suppression profiles
#'
#' @param x An object of class \code{"functional_suppression_profiles_list"}.
#' @param ... Additional arguments.
#'
#' @export
print.functional_suppression_profiles_list <- function(x, ...) {
  cat("List of Functional Suppression Profiles per Environment\n")
  cat("----------------------------------------------------\n")
  cat("Environments represented (", length(x), "):\n", sep = "")
  for (name in names(x)) {
    cat("  - ", name, "\n", sep = "")
  }
  cat("\nUse x$<environment_name> to inspect an individual profile.\n")
  invisible(x)
}

#' Summarize a list of functional suppression profiles
#'
#' @param object An object of class \code{"functional_suppression_profiles_list"}.
#' @param ... Additional arguments.
#'
#' @return A list of summaries, one for each environment.
#' @export
summary.functional_suppression_profiles_list <- function(object, ...) {
  res <- lapply(object, summary, ...)
  class(res) <- "summary.functional_suppression_profiles_list"
  res
}

#' Print summary of a list of functional suppression profiles
#'
#' @param x An object of class \code{"summary.functional_suppression_profiles_list"}.
#' @param ... Additional arguments.
#'
#' @export
print.summary.functional_suppression_profiles_list <- function(x, ...) {
  for (name in names(x)) {
    cat("========================================\n")
    cat("Environment: ", name, "\n", sep = "")
    cat("========================================\n")
    print(x[[name]])
    cat("\n")
  }
  invisible(x)
}

#' Plot a list of functional suppression profiles
#'
#' @param x An object of class \code{"functional_suppression_profiles_list"}.
#' @param ... Additional arguments passed to \code{plot.functional_suppression_profiles}.
#'
#' @return A list of ggplot objects, one for each environment.
#' @export
plot.functional_suppression_profiles_list <- function(x, ...) {
  lapply(x, plot, ...)
}
