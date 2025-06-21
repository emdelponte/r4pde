# zzz.R

# Suppress NOTES from R CMD check about non-standard evaluation (NSE) variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".data", "across", "as_label", "bind_cols", "distinct", "element_rect",
    "enquo", "linearHypothesis", "pull", "rename_with", "row_number", "theme",
    "ungroup", "field", "i", "n", "p", "V", "Vbin", "x", "y", "Y", "mod_x", "focus_id",
    "start_day", "end_day", "column_name", "m", "simes_threshold", "score",
    "mk", "treat", "treatment", "Slow", "Sup", ":="
  ))
}
.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("Icens", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Icens", ask = FALSE, update = FALSE)
  }
}
# Import from base and recommended packages to suppress NOTES on undefined functions
#' @importFrom stats aggregate coef complete.cases cor.test lm median sd var
#' @importFrom utils install.packages
#' @importFrom car linearHypothesis
#' @importFrom ggplot2 theme element_rect
#' @importFrom dplyr across distinct ungroup pull rename_with bind_cols row_number
#' @importFrom rlang enquo as_label
NULL
