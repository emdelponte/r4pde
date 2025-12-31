#' Plot method for compare_curves objects
#'
#' @description
#' Produces a paired visualization of environment-adjusted mean curves
#' and the corresponding functional dendrogram.
#'
#' @param x An object of class \code{"r4pde_compare_curves"}.
#' @param label_fun Optional function to modify treatment labels.
#' @param palette Optional named vector of colors for clusters.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A patchwork object combining curves and dendrogram.
#'
#' @export
plot.r4pde_compare_curves <- function(x, label_fun = NULL, palette = NULL, ...){
  if(!requireNamespace("patchwork", quietly = TRUE)) stop("Need patchwork.")
  (plot_curves(x, label_fun = label_fun, palette = palette) |
      plot_dendrogram(x, label_fun = label_fun, palette = palette)) +
    patchwork::plot_annotation(tag_levels = "A")
}
