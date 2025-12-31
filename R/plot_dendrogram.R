#' Plot functional dendrogram of epidemic curves
#'
#' @description
#' Visualizes the hierarchical clustering of treatments based on
#' functional distances among epidemic curves.
#'
#' @param x An object of class \code{"r4pde_compare_curves"}.
#' @param label_fun Optional function to modify treatment labels.
#' @param palette Optional named vector of colors for clusters.
#' @param show_cut Logical; whether to display the cluster cut height.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_dendrogram <- function(
    x,
    label_fun = NULL,
    palette = NULL,
    show_cut = TRUE
){
  stopifnot(inherits(x, "r4pde_compare_curves"))
  if(!requireNamespace("ggplot2", quietly = TRUE)) stop("Need ggplot2.")
  if(!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr.")
  if(!requireNamespace("ggdendro", quietly = TRUE)) stop("Need ggdendro.")

  trt <- x$vars$treatment
  if(is.null(label_fun)) label_fun <- function(z) z

  dd <- ggdendro::dendro_data(stats::as.dendrogram(x$hc), type = "rectangle")

  lab_map <- x$clusters |>
    dplyr::mutate(
      label = .data[[trt]],
      label_short = label_fun(.data[[trt]]),
      cluster = factor(cluster)
    ) |>
    dplyr::select(label, label_short, cluster)

  labs_df <- dd$labels |>
    dplyr::left_join(lab_map, by = "label")

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dd$segments,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.45
    ) +
    ggplot2::geom_text(
      data = labs_df,
      ggplot2::aes(x = x, y = y, label = label_short, colour = cluster),
      angle = 90, hjust = 1, vjust = 0.5, size = 3
    ) +
    ggplot2::labs(x = NULL, y = "Functional distance", colour = "Cluster") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.line.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 18, l = 5.5)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if(show_cut){
    k <- x$settings$cluster_k
    h_cut <- x$hc$height[length(x$hc$height) - (k - 1)]
    p <- p + ggplot2::geom_hline(yintercept = h_cut, linetype = "dashed", linewidth = 0.5)
  }

  if(!is.null(palette)) p <- p + ggplot2::scale_colour_manual(values = palette)
  p
}
