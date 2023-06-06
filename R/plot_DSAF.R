#' Plot Dynamics and Structure Analysif of Foci (plot_DSAF)
#'
#' This function creates a tile plot of the foci (cluster) identified by the DSAF function.
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
#' # Perform the DSAF
#' result <- DSAF(df)
#' # Plot the foci
#' plot_DSAF(result[[3]])
#'
#' @importFrom ggplot2 ggplot geom_tile theme_minimal geom_text
#' @importFrom dplyr group_by summarise
#' @export
plot_DSAF <- function(df) {
  if (!require("ggplot2")) {
    install.packages("ggplot2")
    library("ggplot2")
  }

  # Calculate the centroid of each cluster
  centroids <- df %>%
    group_by(cluster_id) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop")

  # Create scatter plot
  ggplot(df, aes(x = x, y = y)) +
    geom_tile(fill = "green", color =  "black") +  # color points by cluster
    geom_text(data = centroids, aes(label = cluster_id), size = 5, vjust = -0.5) +
    theme_minimal()
}
