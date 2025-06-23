#' Analysis of foci structure and dynamics (AFSD)
#'
#' This function performs the analysis of a simple method introduced by Nelson (1996) and expanded by Laranjeira et al. (1998).
#' The function assumes
#' the dataframe supplied as input has columns 'x', 'y', and 'i', where 'x' and 'y' are spatial coordinates
#' and 'i' is a disease indicator variable (1 if diseased, otherwise 0). The function performs several steps
#' including filtering rows where 'i' is 1, converting to an adjacency matrix, and creating foci using igraph.
#' It then calculates various statistics about the foci and returns these in a list.
#'
#' @param df A dataframe containing at least three columns: 'x', 'y', and 'i'. 'x' and 'y' represent spatial coordinates
#'           and 'i' is a disease indicator (1 if diseased, otherwise 0).
#'
#' @return A list containing:
#'         cluster_summary2: a dataframe summarizing the number and size of foci, and proportions of diseased plants.
#'         cluster_df: a dataframe containing foci information, including size and number of rows and columns in each foci.
#'         df_clustered: the original dataframe with an added 'focus_id' column, showing which foci each row belongs to.
#'
#' @examples
#' # Generate a sample dataframe
#' set.seed(123)
#' df <- data.frame(x = sample(1:100, 500, replace = TRUE),
#'                  y = sample(1:100, 500, replace = TRUE),
#'                  i = sample(0:1, 500, replace = TRUE, prob = c(0.7, 0.3)))
#'
#' # Perform the AFSD
#' result <- AFSD(df)
#'
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom tidyr pivot_longer
#' @importFrom stats dist
#'
#' @export
#' @family Spatial analysis
AFSD <- function(df) {
  # Check that the necessary columns are present in the dataframe
  if (!all(c("x", "y", "i") %in% colnames(df))) {
    stop("Dataframe must contain 'x', 'y', and 'i' columns.")
  }

  # Check for required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required. Please install it with install.packages('igraph').")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("The 'tidyr' package is required. Please install it with install.packages('tidyr').")
  }

  # Filter for rows where i == 1
  df_filtered <- df[df$i == 1, c("x", "y")]

  # Check that filtered dataframe is not empty
  if (nrow(df_filtered) == 0) {
    stop("There are no rows where i == 1.")
  }

  # Convert to adjacency matrix
  adj <- as.matrix(dist(df_filtered))

  # Threshold the adjacency matrix
  adj_thresholded <- adj <= sqrt(2)  # Consider cells touching diagonally as well

  # Create graph and find clusters
  g <- igraph::graph_from_adjacency_matrix(adj_thresholded, mode = "undirected")
  clusters <- igraph::components(g)

  # Calculate proportion of diseased plants
  prop_diseased <- sum(df$i) / nrow(df)

  # Calculate the size of each cluster
  cluster_sizes <- table(clusters$membership)

  # Create a dataframe with cluster information and sizes
  cluster_df <- data.frame(
    focus_id = names(cluster_sizes),
    size = as.integer(cluster_sizes)
  )

  cluster_summary <- data.frame(
    NF = length(cluster_sizes),
    NF1000 = length(cluster_sizes) * 1000 / nrow(df),
    NSF = sum(cluster_sizes == 1),
    NSF1000 = sum(cluster_sizes == 1) * 1000 / nrow(df),
    DIS_INC = prop_diseased
  )

  # Add columns for the number of rows and columns in each cluster
  cluster_df$rows <- sapply(cluster_df$focus_id, function(id) {
    nodes_in_cluster <- which(clusters$membership == as.integer(id))
    max(df_filtered[nodes_in_cluster, 'y']) - min(df_filtered[nodes_in_cluster, 'y']) + 1
  })
  cluster_df$cols <- sapply(cluster_df$focus_id, function(id) {
    nodes_in_cluster <- which(clusters$membership == as.integer(id))
    max(df_filtered[nodes_in_cluster, 'x']) - min(df_filtered[nodes_in_cluster, 'x']) + 1
  })

  cluster_df$SIF <- cluster_df$rows / cluster_df$cols
  cluster_df$CIF <- cluster_df$size / (cluster_df$rows * cluster_df$cols)

  # Add mean indices for shape and compactness
  cluster_summary$mean_SIF <- mean(cluster_df$cols / cluster_df$rows)
  cluster_summary$mean_CIF <- mean(cluster_df$size / (cluster_df$rows * cluster_df$cols))

  # Pivot the summary table
  cluster_summary2 <- tidyr::pivot_longer(cluster_summary, cols = 1:7, names_to = "stats", values_to = "value")

  # Add cluster IDs to filtered df
  df_clustered <- df[df$i == 1, ]
  df_clustered$focus_id <- clusters$membership

  return(list(cluster_summary2, cluster_df, df_clustered))
}




