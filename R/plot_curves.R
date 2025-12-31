#' Plot environment-adjusted epidemic curves by cluster
#'
#' @description
#' Plots environment-adjusted mean epidemic curves for each treatment,
#' colored according to functional cluster membership.
#'
#' The returned object is a \code{ggplot} and can be further modified
#' using standard ggplot2 layers.
#'
#' @param x An object of class \code{"r4pde_compare_curves"}.
#' @param label_fun Optional function to modify treatment labels.
#' @param palette Optional named vector of colors for clusters.
#' @param alpha Line transparency.
#' @param linewidth Line width.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_curves <- function(
    x,
    label_fun = NULL,
    palette = NULL,
    alpha = 0.9,
    linewidth = 1.1
){
  stopifnot(inherits(x, "r4pde_compare_curves"))
  if(!requireNamespace("ggplot2", quietly = TRUE)) stop("Need ggplot2.")
  if(!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr.")

  trt <- x$vars$treatment
  tim <- x$vars$time
  if(is.null(label_fun)) label_fun <- function(z) z

  df <- dplyr::left_join(x$pred, x$clusters, by = trt) |>
    dplyr::mutate(
      trt_lab = label_fun(.data[[trt]]),
      cluster = factor(cluster)
    )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[tim]],
      y = mu,
      group = trt_lab,
      colour = cluster
    )
  ) +
    ggplot2::geom_line(
      linewidth = linewidth,
      alpha = alpha
    ) +
    ggplot2::expand_limits(y = c(0, 1)) +
    ggplot2::labs(
      x = tim,
      y = "Environment-adjusted mean severity",
      colour = "Cluster"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(legend.position = "bottom")

  if(!is.null(palette)){
    p <- p + ggplot2::scale_colour_manual(values = palette)
  }

  p
}
