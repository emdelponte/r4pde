#' Plot ASFD
#'
#' This function creates a tile plot of the foci (cluster) identified by the AFSD function.
#' It colors each cell in a foci and labels the centroid of each cluster with the foci ID.
#' The 'ggplot2' package is used for the plot, and will be automatically installed if not already present.
#'
#' @param df A dataframe containing at least three columns: 'x', 'y', and 'cluster_id'.
#'           'x' and 'y' are spatial coordinates and 'cluster_id' is the cluster identifier
#'           to which each cell belongs.
#'
#' @return A ggplot object with the scatter plot of foci (clusters).
#'
#' @examples
#' df <- data.frame(x = sample(1:100, 500, replace = TRUE),
#'                  y = sample(1:100, 500, replace = TRUE),
#'                  i = sample(0:1, 500, replace = TRUE, prob = c(0.7, 0.3)))
#'
#' # Perform the AFSD
#' result <- AFSD(df)
#' # Plot the foci
#' plot_AFSD(result[[3]])
#'
#' @importFrom ggplot2 ggplot geom_tile theme_minimal geom_text
#' @importFrom dplyr group_by summarise
#' @export
#' @family Spatial analysis
plot_AFSD <- function(df) {
  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required. Please install it with install.packages('ggplot2').")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required. Please install it with install.packages('dplyr').")
  }

  # Calculate the centroid of each cluster
  centroids <- dplyr::group_by(df, focus_id) |>
    dplyr::summarise(x = mean(x), y = mean(y), .groups = "drop")

  # Create scatter plot
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_tile(fill = "gray70", color = "black") +
    ggplot2::geom_text(data = centroids, ggplot2::aes(label = focus_id), size = 5, vjust = -0.5) +
    ggplot2::theme_minimal()
}
