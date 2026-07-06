#' Functional suppression profiles of disease control treatments
#'
#' @description
#' Computes disease suppression trajectories relative to a reference treatment, estimates functional distances among suppression curves, clusters treatments according to their temporal suppression profiles, and summarizes each profile using functional suppression metrics.
#'
#' @details 
#' A Functional Suppression Profile (FSP) is the temporal trajectory of disease 
#' suppression produced by a treatment relative to an untreated or reference control. 
#' Functional distances and hierarchical clustering are computed from smoothed disease 
#' suppression trajectories rather than raw observed suppression values. Functional 
#' summary metrics are subsequently used to interpret the resulting suppression profiles.
#'
#' The workflow includes contrasting treatments against a reference, fitting functional 
#' curves, calculating distances, clustering, and summarizing the suppression using 
#' metrics like protected area, maximum suppression, and persistence.
#' 
#' Treatments belonging to the same functional suppression profile are displayed using 
#' the same colour throughout all graphical summaries. This visual consistency facilitates 
#' interpretation of the relationship between suppression trajectories, functional clustering, 
#' and suppression metrics.
#' 
#' The \code{classification} table combines functional profile membership, suppression metrics, 
#' and mean functional rank. Functional profiles are defined from distances among smoothed 
#' suppression trajectories, whereas mean rank summarizes the overall performance across 
#' functional suppression metrics.
#'
#' @param data A data frame or tibble containing the disease progress data.
#' @param reference Character string specifying the reference treatment name.
#' @param time Character string specifying the time column. Default is \code{"time"}.
#' @param response Character string specifying the response column. Default is \code{"severity"}.
#' @param treatment Character string specifying the treatment column. Default is \code{"treatment"}.
#' @param environment Character string specifying the environment column (optional).
#' @param by_environment Logical; if \code{TRUE}, the analysis is run separately for each environment. Default is \code{FALSE}.
#' @param env_ref Character string specifying the reference environment for joint analysis (optional).
#' @param threshold Numeric threshold for the persistence metric. Default is \code{5}.
#' @param metrics Character vector of metrics to include in the ranking. Default is \code{c("protected_area", "max_suppression", "persistence", "centroid")}.
#' @param dist_method Character string specifying the distance method. Default is \code{"euclidean"}.
#' @param hclust_method Character string specifying the hierarchical clustering method. Default is \code{"average"}.
#' @param k Integer specifying the number of functional profiles to create, or \code{"auto"} (default). When \code{"auto"}, the optimal number of profiles is estimated by maximizing the average silhouette width.
#' @param cut_height Numeric specifying the cut height for clustering. If provided, overrides \code{k}.
#' @param ... Additional arguments passed to \code{\link{functional_contrast}}.
#'
#' @return A list of class \code{"functional_suppression_profiles"} containing:
#' \itemize{
#'   \item \code{call}: The matched call.
#'   \item \code{reference}: The reference treatment name.
#'   \item \code{contrast}: The contrast data from \code{functional_contrast}.
#'   \item \code{dsp_curves}: The fitted curves from \code{functional_curves}.
#'   \item \code{distances}: The distances object from \code{functional_distances}.
#'   \item \code{hclust}: The hierarchical clustering object.
#'   \item \code{clusters}: A tibble with treatment cluster assignments (internal).
#'   \item \code{profiles}: A lightweight tibble mapping treatments to profiles.
#'   \item \code{classification}: A tibble with the final classification combining profiles, mean rank, and functional metrics, ordered by mean rank.
#'   \item \code{summary}: A tibble with functional summary metrics.
#'   \item \code{ranking}: A tibble with metric rankings.
#'   \item \code{cluster_summary}: A tibble combining metrics, ranks, and cluster assignments.
#'   \item \code{parameters}: A list of input parameters.
#' }
#'
#' @export
#' @seealso \code{\link{functional_contrast}}, \code{\link{functional_curves}}, \code{\link{functional_distances}}, \code{\link{functional_summary}}, \code{\link{rank_dsp}}, \code{\link{plot_dendrogram}}
#'
#' @examples
#' \dontrun{
#' sim_dat <- tibble::tibble(
#'   treatment = rep(c("Control", "A", "B", "C"), each = 6),
#'   time = rep(seq(0, 25, by = 5), times = 4),
#'   severity = c(
#'     c(5, 10, 20, 35, 50, 65),
#'     c(3, 5, 10, 18, 30, 40),
#'     c(4, 6, 12, 22, 35, 45),
#'     c(2, 4, 8, 14, 22, 32)
#'   )
#' )
#' fsp <- functional_suppression_profiles(
#'   data = sim_dat,
#'   reference = "Control",
#'   time = "time",
#'   response = "severity",
#'   treatment = "treatment",
#'   threshold = 5,
#'   k = 2
#' )
#'
#' fsp
#' summary(fsp)
#' plot(fsp, type = "dendrogram")
#' plot(fsp, type = "profiles")
#' plot(fsp, type = "heatmap")
#' plot(fsp, type = "rank")
#' }
functional_suppression_profiles <- function(
    data,
    reference,
    time = "time",
    response = "severity",
    treatment = "treatment",
    environment = NULL,
    by_environment = FALSE,
    env_ref = NULL,
    threshold = 5,
    metrics = c("protected_area", "max_suppression", "persistence", "centroid"),
    dist_method = "euclidean",
    hclust_method = "average",
    k = "auto",
    cut_height = NULL,
    ...
) {
  
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Package 'stats' is required.")
  }

  if (isTRUE(by_environment)) {
    if (is.null(environment)) {
      stop("`environment` must be specified if `by_environment = TRUE`.")
    }
    if (!environment %in% names(data)) {
      stop(sprintf("Environment column '%s' not found in data.", environment))
    }
    
    env_levels <- unique(as.character(data[[environment]]))
    res_list <- lapply(env_levels, function(env_level) {
      env_data <- data[data[[environment]] == env_level, , drop = FALSE]
      functional_suppression_profiles(
        data = env_data,
        reference = reference,
        time = time,
        response = response,
        treatment = treatment,
        environment = NULL,
        by_environment = FALSE,
        env_ref = NULL,
        threshold = threshold,
        metrics = metrics,
        dist_method = dist_method,
        hclust_method = hclust_method,
        k = k,
        cut_height = cut_height,
        ...
      )
    })
    names(res_list) <- env_levels
    class(res_list) <- "functional_suppression_profiles_list"
    return(res_list)
  }

  # Ensure treatment is character to avoid issues
  if (is.factor(data[[treatment]])) {
    data[[treatment]] <- as.character(data[[treatment]])
  }
  
  # Step 1: Functional contrast
  contrast <- functional_contrast(
    data = data,
    reference = reference,
    time = time,
    response = response,
    treatment = treatment,
    group = environment,
    ...
  )
  
  # Step 2: Functional curves for DSP
  # functional_curves uses character variables for column names
  dsp_curves <- functional_curves(
    data = contrast,
    time = time,
    response = "DSP",
    treatment = treatment,
    environment = environment,
    env_ref = env_ref,
    ...
  )
  
  # Determine k for functional_distances and cutree
  n_treatments <- length(unique(dsp_curves$curves[[treatment]]))
  
  if (is.null(cut_height)) {
    if (is.null(k) || identical(k, "auto")) {
      if (n_treatments >= 6) {
        pass_k <- 4
      } else {
        pass_k <- min(3, max(1, n_treatments - 1))
      }
    } else {
      pass_k <- min(k, n_treatments - 1)
    }
    pass_k <- max(1, pass_k)
  } else {
    pass_k <- min(4, n_treatments)
  }

  # Step 3: Functional distances
  distances <- functional_distances(
    dsp_curves,
    hc_method = hclust_method,
    cluster_k = pass_k
  )
  
  # Step 4: Hierarchical clustering
  dist_mat <- distances$distance_matrix
  if (is.null(dist_mat) && !is.null(distances$curve_distance)) dist_mat <- distances$curve_distance
  if (is.null(dist_mat)) dist_mat <- as.matrix(stats::dist(dsp_curves$curves_summary))
  
  if (!is.null(distances$hc)) {
    hc <- distances$hc
  } else if (!is.null(distances$hclust)) {
    hc <- distances$hclust
  } else {
    hc <- stats::hclust(stats::as.dist(dist_mat), method = hclust_method)
  }
  
  # Step 5: Cluster assignment
  sil_obj <- NULL
  sil_avg <- NA
  
  if (!is.null(cut_height)) {
    clusters <- stats::cutree(hc, h = cut_height)
    final_k <- length(unique(clusters))
  } else if (identical(k, "auto")) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      warning("Package 'cluster' is required for automatic profile estimation. Defaulting to k = 3.")
      final_k <- min(3, max(1, n_treatments - 1))
      clusters <- stats::cutree(hc, k = final_k)
    } else {
      max_k <- min(10, n_treatments - 1)
      if (max_k < 2) {
        final_k <- 1
        clusters <- stats::cutree(hc, k = final_k)
      } else {
        k_values <- 2:max_k
        avg_widths <- numeric(length(k_values))
        for (i in seq_along(k_values)) {
          clust <- stats::cutree(hc, k = k_values[i])
          sil <- cluster::silhouette(clust, stats::as.dist(dist_mat))
          avg_widths[i] <- mean(sil[, "sil_width"])
        }
        best_idx <- which.max(avg_widths)
        final_k <- k_values[best_idx]
        clusters <- stats::cutree(hc, k = final_k)
        sil_obj <- cluster::silhouette(clusters, stats::as.dist(dist_mat))
        sil_avg <- avg_widths[best_idx]
      }
    }
  } else {
    final_k <- k
    clusters <- stats::cutree(hc, k = final_k)
    if (final_k > 1 && requireNamespace("cluster", quietly = TRUE)) {
      sil_obj <- cluster::silhouette(clusters, stats::as.dist(dist_mat))
      sil_avg <- mean(sil_obj[, "sil_width"])
    }
  }
  
  cluster_table <- tibble::tibble(
    !!treatment := names(clusters),
    cluster = clusters
  )
  
  # Convert numeric clusters to Roman numerals
  k_unique <- length(unique(cluster_table$cluster))
  cluster_labels <- paste("Profile", utils::as.roman(cluster_table$cluster))
  cluster_table$cluster <- factor(cluster_labels, levels = paste("Profile", utils::as.roman(sort(unique(cluster_table$cluster)))))
  
  # Generate palette
  if (requireNamespace("scales", quietly = TRUE)) {
    cluster_palette <- scales::hue_pal()(k_unique)
  } else {
    cluster_palette <- grDevices::rainbow(k_unique)
  }
  names(cluster_palette) <- levels(cluster_table$cluster)
  
  # Update distances object to use roman numerals
  if (!is.null(distances$clusters)) {
    distances$clusters$cluster <- cluster_table$cluster[match(distances$clusters[[treatment]], cluster_table[[treatment]])]
  }
  
  # Ensure treatment is character in cluster table
  cluster_table$treatment <- as.character(cluster_table$treatment)
  
  # Step 6: Functional summary
  if (!is.null(environment)) {
    summary_df <- functional_summary(
      object = dsp_curves$curves,
      time = time,
      response = "mu",
      treatment = treatment,
      threshold = threshold
    )
  } else {
    summary_df <- functional_summary(
      object = contrast,
      time = time,
      response = "DSP",
      treatment = treatment,
      threshold = threshold
    )
  }
  
  # Step 7: Ranking
  # Only use metrics that exist in summary_df
  available_metrics <- intersect(metrics, names(summary_df))
  if (length(available_metrics) < length(metrics)) {
    missing <- setdiff(metrics, names(summary_df))
    warning("Some metrics are missing from the summary and will be dropped: ", paste(missing, collapse = ", "))
  }
  
  if (length(available_metrics) > 0) {
    ranking <- rank_dsp(
      summary_data = summary_df,
      metrics = available_metrics
    )
    
    # Extract only ranking variables and rename average_rank to mean_rank
    rank_cols <- intersect(names(ranking), c(paste0("rank_", available_metrics), "average_rank"))
    ranking <- ranking |>
      dplyr::select(dplyr::all_of(c(treatment, rank_cols)))
    
    if ("average_rank" %in% names(ranking)) {
      ranking <- ranking |> dplyr::rename(mean_rank = average_rank)
    }
  } else {
    ranking <- tibble::tibble(!!treatment := summary_df[[treatment]])
  }
  
  # Ensure treatment column is character in all tables before joining
  summary_df[[treatment]] <- as.character(summary_df[[treatment]])
  if (treatment %in% names(ranking)) {
     ranking[[treatment]] <- as.character(ranking[[treatment]])
  }
  
  # Step 8: Cluster-level summary
  cluster_summary <- summary_df |>
    dplyr::inner_join(ranking, by = treatment) |>
    dplyr::inner_join(cluster_table, by = "treatment")
  
  # Create lightweight profiles table
  profiles_df <- tibble::tibble(
    !!treatment := cluster_table[[treatment]],
    profile = cluster_table$cluster
  )
  
  # Step 9: Final classification table
  # Join profile, mean_rank, and functional metrics
  class_cols <- c(treatment, "profile", "mean_rank", available_metrics)
  
  classification_df <- profiles_df |>
    dplyr::inner_join(ranking |> dplyr::select(dplyr::all_of(c(treatment, "mean_rank"))), by = treatment) |>
    dplyr::inner_join(summary_df |> dplyr::select(dplyr::all_of(c(treatment, available_metrics))), by = treatment) |>
    dplyr::arrange(.data[["mean_rank"]])
  
  # Return object
  res <- structure(
    list(
      call = match.call(),
      reference = reference,
      contrast = contrast,
      dsp_curves = dsp_curves,
      distances = distances,
      hclust = hc,
      k = final_k,
      silhouette_avg = sil_avg,
      silhouette = sil_obj,
      classification = classification_df,
      clusters = cluster_table,
      profiles = profiles_df,
      cluster_palette = cluster_palette,
      summary = summary_df,
      ranking = ranking,
      cluster_summary = cluster_summary,
      parameters = list(
        time = time,
        response = response,
        treatment = treatment,
        environment = environment,
        by_environment = by_environment,
        env_ref = env_ref,
        threshold = threshold,
        metrics = metrics,
        dist_method = dist_method,
        hclust_method = hclust_method,
        k = k,
        cut_height = cut_height
      )
    ),
    class = "functional_suppression_profiles"
  )
}
