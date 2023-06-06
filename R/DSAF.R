#' Dynamics and Structure Analysis of Foci (DSAF)
#'
#' This function performs a spatial analysis for the dynamics and structure analysis of foci. The function assumes
#' the dataframe supplied as input has columns 'x', 'y', and 'i', where 'x' and 'y' are spatial coordinates
#' and 'i' is a disease indicator variable (1 if diseased, otherwise 0). The function performs several steps
#' including filtering rows where 'i' is 1, converting to an adjacency matrix, and creating clusters using igraph.
#' It then calculates various statistics about the clusters and returns these in a list.
#'
#' @param df A dataframe containing at least three columns: 'x', 'y', and 'i'. 'x' and 'y' represent spatial coordinates
#'           and 'i' is a disease indicator (1 if diseased, otherwise 0).
#'
#' @return A list containing:
#'         cluster_summary2: a dataframe summarizing the number and size of clusters, and proportions of diseased plants.
#'         cluster_df: a dataframe containing cluster information, including size and number of rows and columns in each cluster.
#'         df_clustered: the original dataframe with an added 'cluster_id' column, showing which cluster each row belongs to.
#'
#' @examples
#' # Generate a sample dataframe
#' set.seed(123)
#' df <- data.frame(x = sample(1:100, 500, replace = TRUE),
#'                  y = sample(1:100, 500, replace = TRUE),
#'                  i = sample(0:1, 500, replace = TRUE, prob = c(0.7, 0.3)))
#'
#' # Perform the disease spread analysis
#' result <- DSAF(df)
#'
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom tidyr pivot_longer
#' @importFrom stats dist
#'
#' @export
DSAF <- function(df) {
  # Check that the necessary columns are present in the dataframe
  if (!all(c("x", "y", "i") %in% colnames(df))) {
    stop("Dataframe must contain 'x', 'y', and 'i' columns.")
  }

  # Import the necessary libraries
  if (!require("igraph")) {
    install.packages("igraph")
    library("igraph")
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
  g <- graph_from_adjacency_matrix(adj_thresholded, mode = "undirected")
  clusters <- components(g)

  # Calculate proportion of diseased plants
  prop_diseased <- sum(df$i) / nrow(df)

  # Calculate the size of each cluster
  cluster_sizes <- table(clusters$membership)

  # Create a dataframe with cluster information and sizes
  cluster_df <- data.frame(cluster_id = names(cluster_sizes),
                           size = as.integer(cluster_sizes))


  cluster_summary <- data.frame(num_clusters = length(cluster_sizes),
                                num_clusters_thousand = length(cluster_sizes) * 1000 / nrow(df),
                                num_single_cell_clusters = sum(cluster_sizes == 1),
                                num_single_cell_clusters_thousand = sum(cluster_sizes == 1) * 1000 / nrow(df),
                                prop_diseased = prop_diseased)

  # Add columns for the number of rows and columns in each cluster
  cluster_df$rows <- sapply(cluster_df$cluster_id, function(id) {
    nodes_in_cluster <- which(clusters$membership == id)
    max(df_filtered[nodes_in_cluster, 'y']) - min(df_filtered[nodes_in_cluster, 'y']) + 1
  })
  cluster_df$cols <- sapply(cluster_df$cluster_id, function(id) {
    nodes_in_cluster <- which(clusters$membership == id)
    max(df_filtered[nodes_in_cluster, 'x']) - min(df_filtered[nodes_in_cluster, 'x']) + 1
  })

  # Add a mean index for the shape of cluster
  cluster_summary$mean_shape_index <- mean(cluster_df$cols / cluster_df$rows)

  # Add a mean index for the compactness of the cluster
  cluster_summary$mean_compactness_index <- mean(cluster_df$size / (cluster_df$rows * cluster_df$cols))

  cluster_summary2 <- cluster_summary |>
    pivot_longer(1:7, names_to = "stats", values_to = "value")

  df_clustered <- df[df$i == 1, ]  # only keep rows where i == 1
  df_clustered$cluster_id <- clusters$membership  # add cluster ID to the dataframe



  return(list(cluster_summary2, cluster_df, df_clustered))
}

