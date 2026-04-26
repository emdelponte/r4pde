#' Compute functional resistance and stability-adjusted functional resistance
#'
#' @description
#' Computes a functional resistance index (FRI) relative to a reference genotype curve.
#' Optionally adjusts FRI by a functional instability penalty to calculate a
#' stability-adjusted functional resistance index (SAFRI). Supports grouping
#' genotypes into classes based on quantiles, clustering, or bootstrap-supported differences.
#'
#' @param object An object of class \code{functional_curves}.
#' @param instability Optional output from \code{functional_instability()}.
#' @param reference Character string naming the reference genotype.
#' @param lambda Numeric penalty weight for instability. Default is 1.
#' @param method Character string for distance method.
#' @param scale_nfi Logical; if \code{TRUE}, normalizes nFI to a 0-1 scale before applying penalty.
#' @param n_groups Integer number of classes to group into.
#' @param group_method Character string for the grouping method.
#' @param time_var Optional variable name overrides.
#' @param genotype_var Optional variable name overrides.
#' @param n_boot Integer number of bootstrap iterations.
#' @param ci_level Numeric confidence level for bootstrap intervals.
#' @param clustering_method Character string for clustering method if \code{group_method = "clustering"}.
#' @param ... Additional arguments.
#'
#' @return A list with a \code{table} containing genotype, FRI, nFI, SAFRI, rank, and classes,
#' and optionally bootstrap results.
#'
#' @examples
#' \dontrun{
#' fr <- functional_resistance(
#'   fc,
#'   reference = "Susceptible",
#'   group_method = "bootstrap"
#' )
#' print(fr)
#' }
#' @export
functional_resistance <- function(
    object,
    instability = NULL,
    reference,
    lambda = 1,
    method = c("positive_area", "l2_difference"),
    scale_nfi = FALSE,
    n_groups = 4,
    group_method = c("none", "quantile", "bootstrap", "clustering"),
    time_var = NULL,
    genotype_var = NULL,
    n_boot = 1000,
    ci_level = 0.95,
    clustering_method = c("hclust", "kmeans"),
    ...
) {
  if (!inherits(object, "functional_curves")) {
    stop("`object` must be of class 'functional_curves'")
  }

  method <- match.arg(method)
  group_method <- match.arg(group_method)
  clustering_method <- match.arg(clustering_method)

  .time <- if(is.null(time_var)) object$vars$time else time_var
  .trt <- if(is.null(genotype_var)) object$vars$treatment else genotype_var

  pred_trt <- object$curves
  t_grid <- object$grid
  dt_grid <- mean(diff(t_grid))

  if (!reference %in% pred_trt[[.trt]]) {
    stop("`reference` not found in the treatments of `object`.")
  }

  ref_curve <- pred_trt |>
    dplyr::filter(.data[[.trt]] == reference) |>
    dplyr::arrange(.data[[.time]]) |>
    dplyr::pull("mu")

  # Compute FRI
  mat <- pred_trt |>
    dplyr::select(.data[[.time]], .data[[.trt]], .data[["mu"]]) |>
    tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["mu"]]) |>
    dplyr::arrange(.data[[.time]])

  X <- as.matrix(mat[, setdiff(names(mat), .time), drop = FALSE])
  trts <- colnames(X)

  fri_res <- numeric(length(trts))
  names(fri_res) <- trts

  for (i in seq_along(trts)) {
    y_g <- X[, i]
    if (method == "positive_area") {
      # area of positive difference
      diff_g <- pmax(ref_curve - y_g, 0)
      fri_res[i] <- sum(diff_g) * dt_grid
    } else if (method == "l2_difference") {
      fri_res[i] <- sqrt(sum((ref_curve - y_g)^2) * dt_grid)
    }
  }

  res_df <- tibble::tibble(
    !!.trt := names(fri_res),
    FRI = fri_res
  )

  if (!is.null(instability)) {
    if (is.data.frame(instability)) {
      inst_df <- instability
    } else if ("nfi_table" %in% names(instability)) {
      inst_df <- instability$nfi_table
    } else {
      stop("Could not extract instability data frame.")
    }
    
    trt_col_inst <- NULL
    for (nm in names(inst_df)) {
      if (all(res_df[[.trt]] %in% inst_df[[nm]])) {
        trt_col_inst <- nm
        break
      }
    }
    if (is.null(trt_col_inst)) {
      stop("Could not match genotype names in `instability`.")
    }

    inst_sub <- inst_df |> dplyr::select(!!trt_col_inst, nFI = dplyr::any_of("nFI"))
    res_df <- res_df |>
      dplyr::left_join(inst_sub, by = stats::setNames(trt_col_inst, .trt))
    
    if(scale_nfi && max(res_df$nFI, na.rm=TRUE) > 0) {
      res_df$nFI <- res_df$nFI / max(res_df$nFI, na.rm=TRUE)
    }
    
    res_df$SAFRI <- res_df$FRI / (1 + lambda * res_df$nFI)
  } else {
    res_df$SAFRI <- res_df$FRI
  }

  res_df <- res_df |>
    dplyr::arrange(dplyr::desc(.data[["SAFRI"]])) |>
    dplyr::mutate(rank = dplyr::row_number())

  # Grouping
  default_labels <- c("high stable resistance", "moderate stable resistance", 
                      "low stable resistance", "susceptible or unstable")
  if (n_groups == 4) {
    grp_labels <- default_labels
  } else {
    grp_labels <- paste("Class", seq_len(n_groups))
  }

  if (group_method == "quantile") {
    res_df$resistance_class <- as.character(cut(
      res_df$SAFRI,
      breaks = stats::quantile(res_df$SAFRI, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE),
      labels = rev(grp_labels),
      include.lowest = TRUE
    ))
    res_df$resistance_class <- factor(res_df$resistance_class, levels = grp_labels)
  } else if (group_method == "clustering") {
    if (clustering_method == "hclust") {
      hc <- stats::hclust(stats::dist(res_df$SAFRI), method = "ward.D2")
      cl <- stats::cutree(hc, k = n_groups)
    } else {
      km <- stats::kmeans(res_df$SAFRI, centers = n_groups)
      cl <- km$cluster
    }
    
    cl_means <- tapply(res_df$SAFRI, cl, mean)
    cl_order <- order(cl_means, decreasing = TRUE)
    map_cl <- match(cl, cl_order)
    
    res_df$resistance_class <- factor(grp_labels[map_cl], levels = grp_labels)
  }

  boot_out <- NULL

  if (group_method == "bootstrap") {
    # Perform bootstrap of predicted curves
    df_boot <- object$observed_data
    m_gam <- object$gam
    .env <- object$vars$environment
    .blk <- object$vars$block
    unit_used <- object$vars$unit
    
    ids0 <- if(!is.null(.blk)) levels(df_boot[[.blk]]) else levels(df_boot[[unit_used]])
    id_col <- if(!is.null(.blk)) .blk else unit_used
    
    # We need curve grid
    curve_meta <- df_boot |>
      dplyr::distinct(.data[[unit_used]], .data[[.trt]], .data[[id_col]]) |>
      dplyr::mutate(
        curve_id_chr = as.character(.data[[unit_used]]),
        trt_chr      = as.character(.data[[.trt]])
      )
      
    new_curve_grid <- tidyr::crossing(
      curve_id_chr = curve_meta$curve_id_chr,
      !!.time := t_grid
    ) |>
      dplyr::left_join(curve_meta, by = "curve_id_chr") |>
      dplyr::mutate(
        !!unit_used := factor(curve_id_chr, levels = levels(df_boot[[unit_used]])),
        !!.trt      := factor(trt_chr, levels = levels(df_boot[[.trt]]))
      )
      
    excl_c <- character()
    if(any(grepl(sprintf("s\\(%s\\)", unit_used), names(m_gam$smooth), fixed = TRUE))) {
       excl_c <- c(excl_c, sprintf("s(%s)", unit_used))
    }
    
    new_curve_grid$y_hat <- as.numeric(
      mgcv::predict.gam(m_gam, newdata = new_curve_grid, type = "response", exclude = excl_c)
    )
    
    boot_safri_list <- vector("list", n_boot)
    
    for (b in seq_len(n_boot)) {
      ids_b <- sample(ids0, replace = TRUE)
      w_draw <- tibble::tibble(id = factor(ids_b, levels = ids0)) |> dplyr::count(id, name = "w")
      grid_b <- new_curve_grid |>
        dplyr::mutate(id = factor(.data[[id_col]], levels = ids0)) |>
        dplyr::left_join(w_draw, by = "id") |>
        dplyr::filter(!is.na(.data[["w"]]))
        
      pred_b <- grid_b |>
        dplyr::group_by(.data[[.trt]], .data[[.time]]) |>
        dplyr::summarise(mu = stats::weighted.mean(.data[["y_hat"]], w = .data[["w"]], na.rm = TRUE), .groups = "drop")
        
      # Compute FRI for b
      mat_b <- pred_b |>
        tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["mu"]]) |>
        dplyr::arrange(.data[[.time]])
      X_b <- as.matrix(mat_b[, setdiff(names(mat_b), .time), drop = FALSE])
      
      if (!reference %in% colnames(X_b)) next
      ref_b <- X_b[, reference]
      
      fri_b <- numeric(ncol(X_b))
      names(fri_b) <- colnames(X_b)
      for (j in seq_len(ncol(X_b))) {
        if (method == "positive_area") {
          fri_b[j] <- sum(pmax(ref_b - X_b[, j], 0)) * dt_grid
        } else {
          fri_b[j] <- sqrt(sum((ref_b - X_b[, j])^2) * dt_grid)
        }
      }
      
      safri_b <- fri_b
      if (!is.null(instability)) {
        for(nm in names(fri_b)) {
          nfi_val <- res_df$nFI[res_df[[.trt]] == nm]
          if(length(nfi_val)>0) {
            safri_b[nm] <- fri_b[nm] / (1 + lambda * nfi_val[1])
          }
        }
      }
      boot_safri_list[[b]] <- tibble::tibble(!!.trt := names(safri_b), SAFRI_b = safri_b, iter = b)
    }
    
    boot_safri_df <- dplyr::bind_rows(boot_safri_list)
    
    # intervals
    ci_low <- (1 - ci_level)/2
    ci_high <- 1 - ci_low
    
    boot_ci_tbl <- boot_safri_df |>
      dplyr::group_by(.data[[.trt]]) |>
      dplyr::summarise(
        SAFRI_lower = stats::quantile(.data[["SAFRI_b"]], probs = ci_low, na.rm = TRUE),
        SAFRI_upper = stats::quantile(.data[["SAFRI_b"]], probs = ci_high, na.rm = TRUE),
        .groups = "drop"
      )
      
    res_df <- res_df |> dplyr::left_join(boot_ci_tbl, by = .trt)
    
    # pairwise comparisons
    trts_ordered <- res_df[[.trt]]
    n_t <- length(trts_ordered)
    group_labels <- rep(NA, n_t)
    
    # Simple grouping algorithm based on bootstrap CIs of differences
    # We assign classes A, B, C, D (or 1, 2, 3...)
    class_names <- LETTERS[1:n_groups]
    if (n_groups == 4) {
      class_names <- c("A", "B", "C", "D")
      # A: high stability-adjusted resistance, B: intermediate-high, C: intermediate-low, D: susceptible or unstable
    }
    
    # Convert boot_safri_df to matrix
    mat_boot <- boot_safri_df |> tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["SAFRI_b"]])
    
    current_class <- 1
    group_assign <- numeric(n_t)
    
    for (i in 1:n_t) {
      if (group_assign[i] > 0) next
      group_assign[i] <- current_class
      for (j in (i+1):n_t) {
        if (j > n_t) break
        if (group_assign[j] > 0) next
        
        diff_ij <- mat_boot[[trts_ordered[i]]] - mat_boot[[trts_ordered[j]]]
        ci_diff <- stats::quantile(diff_ij, probs = c(ci_low, ci_high), na.rm = TRUE)
        
        if (ci_diff[1] <= 0 && ci_diff[2] >= 0) {
          # Includes zero -> cannot separate, assign same class
          group_assign[j] <- current_class
        }
      }
      current_class <- current_class + 1
    }
    
    # Normalize class indices to max n_groups
    group_assign <- pmin(group_assign, n_groups)
    res_df$bootstrap_group <- class_names[group_assign]
    res_df$n_boot <- n_boot
    
    res_df$resistance_class <- factor(res_df$bootstrap_group, levels = class_names)
    boot_out <- list(boot_distributions = boot_safri_df)
  }

  out <- list(
    table = res_df,
    settings = list(lambda = lambda, method = method, scale_nfi = scale_nfi, group_method = group_method),
    bootstrap = boot_out
  )
  class(out) <- "functional_resistance"
  out
}

#' Print functional_resistance
#' @export
print.functional_resistance <- function(x, ...) {
  cat("A functional_resistance object\n")
  print(x$table)
  invisible(x)
}
