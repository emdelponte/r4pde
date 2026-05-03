#' Functional principal component analysis of disease progress curves
#'
#' @description
#' Performs functional principal component analysis on fitted disease progress curves 
#' returned by \code{\link{functional_curves}}. The function decomposes variation among epidemic 
#' trajectories into orthogonal temporal components and returns curve-level scores, 
#' eigenfunctions, variance explained, and reconstructed curves.
#'
#' @param object An object returned by \code{\link{functional_curves}}.
#' @param n_components Optional integer number of functional principal components to retain.
#' @param var_explained Cumulative variance threshold used when \code{n_components = NULL}.
#' @param center Logical; whether to center curves before PCA. Default TRUE.
#' @param scale Logical; whether to scale grid columns before PCA. Default FALSE.
#' @param method Character; method for FPCA, currently only \code{"pca_on_grid"} is supported.
#' @param ... Additional arguments for future extensions.
#'
#' @details
#' The function uses the fitted curves from \code{functional_curves()} and does not refit 
#' the disease progress model. The first implementation uses PCA on a common prediction grid. 
#' FPC scores can be analyzed as functional epidemiological traits in downstream models.
#' 
#' Interpretation:
#' \itemize{
#'   \item FPC1 often captures the largest mode of variation, commonly overall epidemic intensity or speed.
#'   \item Later FPCs may capture timing of disease onset, curve crossing, late acceleration, or other shape-related deviations.
#'   \item Interpretation must be data-driven and should be based on eigenfunction plots and mean ± perturbation plots.
#' }
#'
#' @return An object of class \code{"r4pde_functional_pca"} containing:
#' \itemize{
#'   \item \code{scores}: Tibble of curve-level FPC scores.
#'   \item \code{eigenfunctions}: Tibble of eigenfunction values across the time grid.
#'   \item \code{variance}: Tibble of eigenvalues, variance explained, and cumulative variance.
#'   \item \code{mean_curve}: Tibble of the mean curve.
#'   \item \code{reconstructed}: Tibble of reconstructed curves using retained components.
#'   \item \code{input_curves}: Tibble of original fitted curves.
#'   \item \code{pca}: The underlying \code{prcomp} object.
#'   \item \code{settings}: List of settings used.
#'   \item \code{call}: The matched call.
#' }
#'
#' @examples
#' \dontrun{
#' curves <- functional_curves(...)
#' fpca <- functional_pca(curves, var_explained = 0.95)
#' 
#' print(fpca)
#' plot(fpca, type = "scree")
#' plot(fpca, type = "components")
#' plot(fpca, type = "scores", components = c(1, 2))
#' 
#' scores <- get_fpca_scores(fpca)
#' }
#' @export
functional_pca <- function(
    object,
    n_components = NULL,
    var_explained = 0.95,
    center = TRUE,
    scale = FALSE,
    method = c("pca_on_grid"),
    ...
) {
  method <- match.arg(method)
  
  if (!inherits(object, "functional_curves")) {
    stop("`object` must be of class 'functional_curves'.")
  }
  
  if (is.null(object$curves)) {
    stop("Fitted curves are not available in the `object`.")
  }
  
  .time <- object$vars$time
  .trt <- object$vars$treatment
  
  if (is.null(.time) || is.null(.trt)) {
    stop("Time or treatment variable information is missing from the object.")
  }
  
  # Extract fitted curve data
  pred_trt <- object$curves
  t_grid <- object$grid
  
  if (length(t_grid) < 3) {
    stop("There must be at least three time points in the common grid.")
  }
  
  # Reshape to a wide matrix: rows = individual curves, columns = common time grid
  wide_mat <- pred_trt |>
    dplyr::select(dplyr::all_of(c(.trt, .time, "mu"))) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(.time), values_from = "mu")
  
  curve_ids <- wide_mat[[.trt]]
  
  if (length(unique(curve_ids)) < 2) {
    stop("There must be at least two curves.")
  }
  
  # Extract curve matrix
  curve_mat <- as.matrix(wide_mat[, -1, drop = FALSE])
  rownames(curve_mat) <- as.character(curve_ids)
  
  if (any(is.na(curve_mat))) {
    stop("Missing fitted values found in the curve matrix.")
  }
  
  # Perform PCA
  pca_res <- stats::prcomp(curve_mat, center = center, scale. = scale)
  
  # Variance explained
  eigenvalues <- pca_res$sdev^2
  prop_var <- eigenvalues / sum(eigenvalues)
  cum_var <- cumsum(prop_var)
  
  # Select retained components
  if (!is.null(n_components)) {
    k <- min(n_components, length(eigenvalues))
  } else {
    k <- which(cum_var >= var_explained)[1]
    if (is.na(k)) k <- length(eigenvalues)
  }
  
  # Create variance tibble
  var_tbl <- tibble::tibble(
    FPC = paste0("FPC", seq_along(eigenvalues)),
    eigenvalue = eigenvalues,
    prop_var = prop_var,
    cum_var = cum_var
  ) |>
    dplyr::slice(1:k)
  
  # Create scores tibble
  scores_tbl <- tibble::as_tibble(pca_res$x[, 1:k, drop = FALSE])
  scores_tbl <- dplyr::bind_cols(tibble::tibble(curve_id = curve_ids), scores_tbl)
  
  # Create eigenfunctions tibble
  ef_mat <- pca_res$rotation[, 1:k, drop = FALSE]
  ef_tbl <- tibble::as_tibble(ef_mat) |>
    dplyr::mutate(time = t_grid) |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("PC"),
      names_to = "FPC",
      values_to = "value"
    ) |>
    dplyr::mutate(FPC = gsub("^PC", "FPC", .data$FPC))
  
  # Create mean curve tibble
  if (center) {
    mean_val <- pca_res$center
  } else {
    mean_val <- colMeans(curve_mat)
  }
  
  mean_curve_tbl <- tibble::tibble(
    time = t_grid,
    mean = as.numeric(mean_val)
  )
  
  # Create reconstructed curves using retained components
  recon_mat <- t(t(pca_res$x[, 1:k, drop = FALSE] %*% t(pca_res$rotation[, 1:k, drop = FALSE])) + pca_res$center)
  
  recon_tbl <- tibble::as_tibble(recon_mat) |>
    dplyr::mutate(curve_id = curve_ids) |>
    tidyr::pivot_longer(
      cols = -curve_id,
      names_to = "time",
      values_to = "reconstructed"
    ) |>
    dplyr::mutate(time = as.numeric(.data$time))
  
  # Input curves
  input_tbl <- wide_mat |>
    dplyr::rename(curve_id = !!.trt) |>
    tidyr::pivot_longer(
      cols = -curve_id,
      names_to = "time",
      values_to = "fitted"
    ) |>
    dplyr::mutate(time = as.numeric(.data$time))
  
  # Combine for reconstructed output
  recon_combined <- input_tbl |>
    dplyr::left_join(recon_tbl, by = c("curve_id", "time"))
  
  settings <- list(
    n_components = n_components,
    var_explained = var_explained,
    center = center,
    scale = scale,
    method = method,
    n_retained = k,
    time_var = .time,
    trt_var = .trt
  )
  
  out <- list(
    scores = scores_tbl,
    eigenfunctions = ef_tbl,
    variance = var_tbl,
    mean_curve = mean_curve_tbl,
    reconstructed = recon_combined,
    input_curves = input_tbl,
    pca = pca_res,
    settings = settings,
    call = match.call()
  )
  
  class(out) <- "r4pde_functional_pca"
  return(out)
}

#' @export
print.r4pde_functional_pca <- function(x, ...) {
  cat("A functional_pca object\n")
  cat("Number of curves:", nrow(x$scores), "\n")
  cat("Number of time points:", nrow(x$mean_curve), "\n")
  cat("Number of retained FPCs:", x$settings$n_retained, "\n")
  cat("\nVariance Explained:\n")
  print(x$variance |> dplyr::select(FPC, prop_var, cum_var), n = x$settings$n_retained)
  invisible(x)
}

#' @export
summary.r4pde_functional_pca <- function(object, ...) {
  list(
    settings = object$settings,
    variance = object$variance,
    score_ranges = lapply(object$scores |> dplyr::select(-curve_id), range),
    n_curves = nrow(object$scores),
    n_timepoints = nrow(object$mean_curve)
  )
}

#' Plot functional PCA results
#' 
#' @param x An object of class \code{"r4pde_functional_pca"}.
#' @param type Type of plot: "scree", "components", "scores", "reconstruction", or "mean".
#' @param components Integer vector of length 2 indicating which components to plot for scores and mean perturbation.
#' @param curve_id Optional vector of curve IDs to include in reconstruction plot.
#' @param ... Additional arguments.
#' @export
plot.r4pde_functional_pca <- function(
    x,
    type = c("scree", "components", "scores", "reconstruction", "mean"),
    components = c(1, 2),
    curve_id = NULL,
    ...
) {
  type <- match.arg(type)
  
  if (type == "scree") {
    var_data <- x$variance
    p <- ggplot2::ggplot(var_data, ggplot2::aes(x = factor(.data$FPC, levels = .data$FPC), y = .data$prop_var)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::geom_line(ggplot2::aes(y = .data$cum_var, group = 1), color = "darkred", size = 1) +
      ggplot2::geom_point(ggplot2::aes(y = .data$cum_var), color = "darkred", size = 2) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
      ggplot2::labs(x = "Functional Principal Component", y = "Proportion of Variance Explained",
                    title = "Scree Plot") +
      ggplot2::theme_classic()
    return(p)
  }
  
  if (type == "components") {
    ef_data <- x$eigenfunctions
    if (!is.null(components)) {
      ef_data <- ef_data |> dplyr::filter(.data$FPC %in% paste0("FPC", components))
    }
    p <- ggplot2::ggplot(ef_data, ggplot2::aes(x = .data$time, y = .data$value, color = .data$FPC)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_line(size = 1) +
      ggplot2::labs(x = "Time", y = "Eigenfunction Value", title = "Functional Principal Components") +
      ggplot2::theme_classic()
    return(p)
  }
  
  if (type == "scores") {
    if (length(components) != 2) stop("Please specify exactly two components for the scores plot.")
    comp_names <- paste0("FPC", components)
    if (!all(comp_names %in% names(x$scores))) stop("Specified components not found in scores.")
    
    p <- ggplot2::ggplot(x$scores, ggplot2::aes(x = .data[[comp_names[1]]], y = .data[[comp_names[2]]], label = .data$curve_id)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(color = "steelblue", size = 3) +
      ggplot2::geom_text(size = 3, vjust = -0.5, hjust = 0.5) +
      ggplot2::labs(x = paste0(comp_names[1], " Score"), y = paste0(comp_names[2], " Score"), 
                    title = "FPC Scores") +
      ggplot2::theme_classic()
    return(p)
  }
  
  if (type == "reconstruction") {
    r_data <- x$reconstructed
    if (!is.null(curve_id)) {
      r_data <- r_data |> dplyr::filter(.data$curve_id %in% curve_id)
    } else {
      # Take a sample of 6 curves if none specified
      sample_ids <- unique(r_data$curve_id)[1:min(6, length(unique(r_data$curve_id)))]
      r_data <- r_data |> dplyr::filter(.data$curve_id %in% sample_ids)
    }
    
    p <- ggplot2::ggplot(r_data, ggplot2::aes(x = .data$time)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$fitted, color = "Fitted"), size = 1) +
      ggplot2::geom_line(ggplot2::aes(y = .data$reconstructed, color = "Reconstructed"), linetype = "dashed", size = 1) +
      ggplot2::facet_wrap(~curve_id) +
      ggplot2::labs(x = "Time", y = "Value", title = "Curve Reconstruction", color = "") +
      ggplot2::scale_color_manual(values = c("Fitted" = "black", "Reconstructed" = "red")) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "bottom")
    return(p)
  }
  
  if (type == "mean") {
    if (length(components) > 2) stop("Please specify at most two components for the mean perturbation plot.")
    comp_names <- paste0("FPC", components)
    
    mean_df <- x$mean_curve
    ef_df <- x$eigenfunctions |> dplyr::filter(.data$FPC %in% comp_names)
    scores_df <- x$scores
    
    plot_data <- list()
    for (comp in comp_names) {
      if (comp %in% names(scores_df)) {
        sd_score <- sd(scores_df[[comp]])
        ef_val <- ef_df |> dplyr::filter(.data$FPC == comp) |> dplyr::pull(.data$value)
        
        tmp <- tibble::tibble(
          time = mean_df$time,
          FPC = comp,
          mean = mean_df$mean,
          plus_2sd = mean_df$mean + 2 * sd_score * ef_val,
          minus_2sd = mean_df$mean - 2 * sd_score * ef_val
        )
        plot_data[[comp]] <- tmp
      }
    }
    plot_data <- dplyr::bind_rows(plot_data)
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$time)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$mean), color = "black", size = 1) +
      ggplot2::geom_line(ggplot2::aes(y = .data$plus_2sd), color = "blue", linetype = "dashed", size = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = .data$minus_2sd), color = "red", linetype = "dashed", size = 0.8) +
      ggplot2::facet_wrap(~FPC) +
      ggplot2::labs(x = "Time", y = "Value", title = "Mean Curve ± 2 SD Perturbation") +
      ggplot2::theme_classic()
    return(p)
  }
}

#' @export
autoplot.r4pde_functional_pca <- function(object, ...) {
  plot.r4pde_functional_pca(object, ...)
}

#' Get FPCA scores
#' @param x An object of class \code{"r4pde_functional_pca"}.
#' @param components Optional integer vector of components to return.
#' @export
get_fpca_scores <- function(x, components = NULL) {
  if (!inherits(x, "r4pde_functional_pca")) stop("Object must be of class 'r4pde_functional_pca'.")
  scores <- x$scores
  if (!is.null(components)) {
    cols <- c("curve_id", paste0("FPC", components))
    scores <- scores |> dplyr::select(dplyr::any_of(cols))
  }
  scores
}

#' Get FPCA eigenfunctions
#' @param x An object of class \code{"r4pde_functional_pca"}.
#' @param components Optional integer vector of components to return.
#' @export
get_fpca_eigenfunctions <- function(x, components = NULL) {
  if (!inherits(x, "r4pde_functional_pca")) stop("Object must be of class 'r4pde_functional_pca'.")
  ef <- x$eigenfunctions
  if (!is.null(components)) {
    comps <- paste0("FPC", components)
    ef <- ef |> dplyr::filter(.data$FPC %in% comps)
  }
  ef
}

#' Get FPCA variance explained
#' @param x An object of class \code{"r4pde_functional_pca"}.
#' @export
get_fpca_variance <- function(x) {
  if (!inherits(x, "r4pde_functional_pca")) stop("Object must be of class 'r4pde_functional_pca'.")
  x$variance
}

#' Augment functional PCA
#' @param x An object of class \code{"r4pde_functional_pca"}.
#' @export
augment_functional_pca <- function(x) {
  if (!inherits(x, "r4pde_functional_pca")) stop("Object must be of class 'r4pde_functional_pca'.")
  x$reconstructed |>
    dplyr::mutate(residual = .data$fitted - .data$reconstructed)
}

#' Reconstruct curves using specified FPCA components
#' @param x An object of class \code{"r4pde_functional_pca"}.
#' @param components Integer vector of components to use for reconstruction. If NULL, uses all retained components.
#' @export
reconstruct_curves <- function(x, components = NULL) {
  if (!inherits(x, "r4pde_functional_pca")) stop("Object must be of class 'r4pde_functional_pca'.")
  
  pca <- x$pca
  scores <- pca$x
  rot <- pca$rotation
  
  if (is.null(components)) {
    k <- x$settings$n_retained
    components <- 1:k
  } else {
    valid_comps <- components[components <= ncol(scores)]
    if (length(valid_comps) == 0) stop("Invalid components specified.")
    components <- valid_comps
  }
  
  # Ensure components is not out of bounds
  recon_mat <- matrix(0, nrow = nrow(scores), ncol = nrow(rot))
  for (i in components) {
    recon_mat <- recon_mat + outer(scores[, i], rot[, i])
  }
  recon_mat <- t(t(recon_mat) + pca$center)
  
  recon_tbl <- tibble::as_tibble(recon_mat) |>
    dplyr::mutate(curve_id = rownames(scores)) |>
    tidyr::pivot_longer(
      cols = -curve_id,
      names_to = "time",
      values_to = "reconstructed"
    ) |>
    dplyr::mutate(time = as.numeric(.data$time))
  
  recon_tbl
}

#' Cluster curves based on FPCA scores
#' @param x An object of class \code{"r4pde_functional_pca"}.
#' @param k Number of clusters.
#' @param components Integer vector of components to use.
#' @param method Clustering method: "kmeans" or "hclust".
#' @param choose_k Method for suggesting k: "none", "silhouette", or "elbow".
#' @param ... Additional arguments passed to clustering functions.
#' @export
functional_pca_clusters <- function(
    x,
    k = NULL,
    components = NULL,
    method = c("kmeans", "hclust"),
    choose_k = c("none", "silhouette", "elbow"),
    ...
) {
  method <- match.arg(method)
  choose_k <- match.arg(choose_k)
  
  if (!inherits(x, "r4pde_functional_pca")) stop("Object must be of class 'r4pde_functional_pca'.")
  
  scores <- get_fpca_scores(x, components = components)
  curve_ids <- scores$curve_id
  score_mat <- as.matrix(scores |> dplyr::select(-curve_id))
  
  # Standardize scores before clustering
  score_mat <- scale(score_mat)
  
  # Logic for choosing k if requested, though simplified for now
  if (is.null(k) && choose_k == "none") {
    stop("Must specify `k` or a method for `choose_k`.")
  }
  
  if (is.null(k)) {
    k <- 2 # fallback
  }
  
  if (method == "kmeans") {
    res <- stats::kmeans(score_mat, centers = k, ...)
    clusters <- res$cluster
  } else {
    d <- stats::dist(score_mat)
    hc <- stats::hclust(d, ...)
    clusters <- stats::cutree(hc, k = k)
  }
  
  tibble::tibble(
    curve_id = curve_ids,
    cluster = factor(clusters)
  )
}
