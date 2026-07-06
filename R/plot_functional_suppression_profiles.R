#' Plot functional suppression profiles
#'
#' @param x An object of class \code{"functional_suppression_profiles"}.
#' @param type Character string specifying the plot type. One of \code{"dendrogram"}, \code{"profiles"}, \code{"heatmap"}, \code{"rank"}, or \code{"all"}.
#' @param show_cut Logical; whether to display the cluster cut height in dendrogram. Default is \code{TRUE}.
#' @param show_points Logical; whether to overlay observed DSP values as points in the profiles plot. Default is \code{FALSE}.
#' @param show_points Logical; whether to overlay observed DSP values as points in the profiles plot. Default is \code{FALSE}.
#' @param ... Additional arguments passed to specific plot functions.
#'
#' @details 
#' The plot method visualizes the components of a Functional Suppression Profile. 
#' The dendrogram defines the Functional Suppression Profiles based on temporal similarity. 
#' The heatmap is ordered by the dendrogram to help interpret these functional groups 
#' using summary metrics.
#'
#' @return A \code{ggplot} object or a list of \code{ggplot} objects.
#'
#' @export
plot.functional_suppression_profiles <- function(
    x,
    type = c("dendrogram", "profiles", "heatmap", "rank", "all"),
    show_cut = TRUE,
    show_points = FALSE,
    ...
) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  # Get dendrogram order
  hc <- x$hclust
  dend_order <- hc$labels[hc$order]

  # Dendrogram
  p_dendrogram <- function() {
    p <- plot_dendrogram(x$distances, show_cut = show_cut, palette = x$cluster_palette, ...)
    p + ggplot2::labs(y = "Functional distance")
  }
  
  # Profiles
  p_profiles <- function() {
    df_smooth <- x$dsp_curves$curves
    df_obs <- x$contrast
    time_var <- x$parameters$time
    trt_var <- x$parameters$treatment
    
    if (!time_var %in% names(df_smooth) || !trt_var %in% names(df_smooth)) {
      stop("Time or treatment variables not found in smoothed curves data.")
    }
    
      df_smooth[[trt_var]] <- factor(df_smooth[[trt_var]], levels = dend_order)
      df_smooth <- df_smooth |> dplyr::left_join(x$classification |> dplyr::select(dplyr::all_of(c(trt_var, "profile"))), by = trt_var)
      if (show_points) {
        df_obs[[trt_var]] <- factor(df_obs[[trt_var]], levels = dend_order)
        df_obs <- df_obs |> dplyr::left_join(x$classification |> dplyr::select(dplyr::all_of(c(trt_var, "profile"))), by = trt_var)
      }
      
      curve_labels <- df_smooth |>
        dplyr::group_by(.data[[trt_var]]) |>
        dplyr::slice_max(.data[[time_var]], n = 1, with_ties = FALSE) |>
        dplyr::ungroup()
        
      p <- ggplot2::ggplot(df_smooth, ggplot2::aes(
        x = .data[[time_var]], 
        y = .data[["mu"]], 
        color = .data[["profile"]], 
        group = .data[[trt_var]]
      )) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_line(linewidth = 0.8, alpha = 0.85, show.legend = FALSE)
      
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = curve_labels,
        ggplot2::aes(label = .data[[trt_var]]),
        nudge_x = 0.5,
        direction = "y",
        hjust = 0,
        segment.size = 0.2,
        segment.alpha = 0.7,
        min.segment.length = 0,
        box.padding = 0.1,
        point.padding = 0.15,
        show.legend = FALSE,
        size = 2.8
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = curve_labels,
        ggplot2::aes(label = .data[[trt_var]]),
        hjust = -0.1,
        size = 2.8,
        show.legend = FALSE
      )
    }
    
    p <- p +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_manual(values = x$cluster_palette) +
      ggplot2::labs(
        x = "Time",
        y = "Disease suppression",
        color = "Functional Profile"
      )
      
    if (show_points) {
      p <- p + ggplot2::geom_point(
        data = df_obs, 
        ggplot2::aes(x = .data[[time_var]], y = .data[["DSP"]]),
        alpha = 0.5,
        show.legend = FALSE
      )
    }
    p
  }
  
  # Heatmap
  p_heatmap <- function() {
    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("Package 'tidyr' is required for heatmap plot.")
    }
    
    rank_df <- x$ranking
    trt_var <- x$parameters$treatment
    n_treatments <- nrow(rank_df)
    
    rank_df[[trt_var]] <- factor(
      rank_df[[trt_var]], 
      levels = rev(dend_order)
    )
    
    rank_long <- rank_df |>
      tidyr::pivot_longer(
        cols = -dplyr::all_of(trt_var),
        names_to = "metric",
        values_to = "rank"
      )
      
    metric_levels <- c("protected_area", "max_suppression", "persistence", "centroid", "mean_rank")
    metric_labels <- c("Area", "Maximum", "Persistence", "Timing", "Rank")
    
    rank_long <- rank_long |>
      dplyr::mutate(
        metric = gsub("^rank_", "", .data[["metric"]]),
        metric = factor(.data[["metric"]], levels = metric_levels, labels = metric_labels)
      )
      
    y_clusters <- x$classification$profile[match(rev(dend_order), x$classification[[trt_var]])]
    y_colors <- x$cluster_palette[as.character(y_clusters)]
    
    ggplot2::ggplot(rank_long, ggplot2::aes(
      x = .data[["metric"]], 
      y = .data[[trt_var]], 
      fill = .data[["rank"]]
    )) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(1, n_treatments)) + 
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = NULL,
        y = "Treatment",
        fill = "Rank"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = ggplot2::element_text(size = 11, color = y_colors)
      )
  }
  
  # Rank
  p_rank <- function() {
    rank_df <- x$classification
    trt_var <- x$parameters$treatment
    
    if (!"mean_rank" %in% names(rank_df)) {
      stop("mean_rank not found in classification data.")
    }
    
    rank_df[[trt_var]] <- factor(
      rank_df[[trt_var]], 
      levels = rank_df[[trt_var]][order(rank_df[["mean_rank"]], decreasing = TRUE)]
    )
    
    ggplot2::ggplot(rank_df, ggplot2::aes(
      x = .data[["mean_rank"]], 
      y = .data[[trt_var]],
      fill = .data[["profile"]]
    )) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = x$cluster_palette) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Mean Rank",
        y = "Treatment"
      )
  }
  
  if (type == "dendrogram") return(p_dendrogram())
  if (type == "profiles") return(p_profiles())
  if (type == "heatmap") return(p_heatmap())
  if (type == "rank") return(p_rank())
  
  if (type == "all") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      message("Package 'patchwork' is required to combine plots. Returning a list of plots instead.")
      out <- list(
        dendrogram = p_dendrogram(),
        profiles = p_profiles(),
        heatmap = p_heatmap()
      )
      return(out)
    } else {
      p1 <- p_profiles()
      p2 <- p_dendrogram()
      p3 <- p_heatmap()
      p4 <- p_rank()
      
      fig <- (p1 | p2) / (p3 | p4)
      
      return(fig + patchwork::plot_annotation(
        tag_levels = "A"
      ))
    }
  }
}
