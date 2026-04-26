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
#' @param reference Character string naming the reference genotype, or a named vector if \code{reference_within_group = TRUE}.
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
#' @param adjust_by Optional character string naming a genotype-level covariate to stratify by.
#' @param reference_within_group Logical; if \code{TRUE}, the reference is evaluated within each \code{adjust_by} group.
#' @param group_within_adjust_by Logical; if \code{TRUE}, resistance classes are assigned within each \code{adjust_by} group.
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
    adjust_by = NULL,
    reference_within_group = FALSE,
    group_within_adjust_by = TRUE,
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

  if (!is.null(adjust_by)) {
    if (!adjust_by %in% names(object$genotype_info)) {
      stop(sprintf("`adjust_by` variable '%s' not found in object$genotype_info.", adjust_by))
    }
    geno_groups <- object$genotype_info |> dplyr::select(dplyr::all_of(c(.trt, adjust_by)))
    names(geno_groups)[2] <- "adjust_group"
  } else {
    geno_groups <- tibble::tibble(!!.trt := unique(pred_trt[[.trt]]), adjust_group = "All")
    adjust_by <- "All"
  }

  groups <- unique(geno_groups$adjust_group)
  ref_curves <- list()

  if (reference_within_group && adjust_by != "All") {
    for (g in groups) {
      ref_g <- if (length(reference) > 1 && !is.null(names(reference))) reference[as.character(g)] else reference[1]
      if (is.na(ref_g) || !ref_g %in% geno_groups[[.trt]][geno_groups$adjust_group == g]) {
        stop(sprintf("Reference genotype '%s' not found in group '%s'. When reference_within_group = TRUE, each group must have a valid reference genotype.", ref_g, g))
      }
      ref_curves[[as.character(g)]] <- pred_trt |> dplyr::filter(.data[[.trt]] == ref_g) |> dplyr::arrange(.data[[.time]]) |> dplyr::pull("mu")
    }
  } else {
    ref_g <- reference[1]
    if (!ref_g %in% pred_trt[[.trt]]) stop(sprintf("Global reference genotype '%s' not found in the treatments of `object`.", ref_g))
    global_ref <- pred_trt |> dplyr::filter(.data[[.trt]] == ref_g) |> dplyr::arrange(.data[[.time]]) |> dplyr::pull("mu")
    for (g in groups) {
      ref_curves[[as.character(g)]] <- global_ref
    }
  }

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
    g_val <- as.character(geno_groups$adjust_group[geno_groups[[.trt]] == trts[i]])
    ref_curve <- ref_curves[[g_val]]
    
    if (method == "positive_area") {
      fri_res[i] <- sum(pmax(ref_curve - y_g, 0)) * dt_grid
    } else {
      fri_res[i] <- sqrt(sum((ref_curve - y_g)^2) * dt_grid)
    }
  }

  res_df <- tibble::tibble(!!.trt := names(fri_res), FRI = fri_res)

  if (!is.null(instability)) {
    inst_df <- if(is.data.frame(instability)) instability else instability$nfi_table
    trt_col_inst <- names(inst_df)[apply(inst_df, 2, function(col) all(res_df[[.trt]] %in% col))][1]
    if (is.na(trt_col_inst)) stop("Could not match genotype names in `instability`.")

    inst_sub <- inst_df |> dplyr::select(!!trt_col_inst, nFI = dplyr::any_of("nFI"))
    res_df <- res_df |> dplyr::left_join(inst_sub, by = stats::setNames(trt_col_inst, .trt))
    
    if(scale_nfi && max(res_df$nFI, na.rm=TRUE) > 0) res_df$nFI <- res_df$nFI / max(res_df$nFI, na.rm=TRUE)
    res_df$SAFRI <- res_df$FRI / (1 + lambda * res_df$nFI)
  } else {
    res_df$SAFRI <- res_df$FRI
  }

  res_df <- res_df |>
    dplyr::arrange(dplyr::desc(.data[["SAFRI"]])) |>
    dplyr::mutate(rank_global = dplyr::row_number())

  if (adjust_by != "All") {
    res_df <- res_df |>
      dplyr::left_join(geno_groups |> dplyr::rename(!!adjust_by := adjust_group), by = .trt) |>
      dplyr::group_by(.data[[adjust_by]]) |>
      dplyr::arrange(dplyr::desc(.data[["SAFRI"]])) |>
      dplyr::mutate(rank_within_group = dplyr::row_number()) |>
      dplyr::ungroup()
  }

  default_labels <- c("high stable resistance", "moderate stable resistance", "low stable resistance", "susceptible or unstable")
  grp_labels <- if (n_groups == 4) default_labels else paste("Class", seq_len(n_groups))
  class_names <- if (n_groups == 4) c("A", "B", "C", "D") else LETTERS[1:n_groups]

  assign_classes <- function(sub_df) {
    if (group_method == "quantile") {
      c_vec <- as.character(cut(sub_df$SAFRI, breaks = stats::quantile(sub_df$SAFRI, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE), labels = rev(grp_labels), include.lowest = TRUE))
      sub_df$resistance_class <- factor(c_vec, levels = grp_labels)
    } else if (group_method == "clustering") {
      if (nrow(sub_df) < n_groups) {
         warning("Not enough genotypes for clustering. Skipping class assignment.")
         sub_df$resistance_class <- NA
         return(sub_df)
      }
      cl <- if (clustering_method == "hclust") {
        stats::cutree(stats::hclust(stats::dist(sub_df$SAFRI), method = "ward.D2"), k = n_groups)
      } else {
        stats::kmeans(sub_df$SAFRI, centers = n_groups)$cluster
      }
      cl_means <- tapply(sub_df$SAFRI, cl, mean)
      map_cl <- match(cl, order(cl_means, decreasing = TRUE))
      sub_df$resistance_class <- factor(grp_labels[map_cl], levels = grp_labels)
    }
    return(sub_df)
  }

  if (group_method %in% c("quantile", "clustering")) {
    if (adjust_by != "All" && group_within_adjust_by) {
      res_list <- split(res_df, res_df[[adjust_by]])
      res_list <- lapply(res_list, assign_classes)
      res_df <- dplyr::bind_rows(res_list) |> dplyr::arrange(dplyr::desc(.data[["SAFRI"]]))
    } else {
      res_df <- assign_classes(res_df)
    }
  }

  boot_out <- NULL
  if (group_method == "bootstrap") {
    df_boot <- object$observed_data
    m_gam <- object$gam
    .blk <- object$vars$block
    unit_used <- object$vars$unit
    ids0 <- if(!is.null(.blk)) levels(df_boot[[.blk]]) else levels(df_boot[[unit_used]])
    id_col <- if(!is.null(.blk)) .blk else unit_used
    
    curve_meta <- df_boot |> dplyr::distinct(.data[[unit_used]], .data[[.trt]], .data[[id_col]]) |>
      dplyr::mutate(curve_id_chr = as.character(.data[[unit_used]]), trt_chr = as.character(.data[[.trt]]))
      
    new_curve_grid <- tidyr::crossing(curve_id_chr = curve_meta$curve_id_chr, !!.time := t_grid) |>
      dplyr::left_join(curve_meta, by = "curve_id_chr") |>
      dplyr::mutate(!!unit_used := factor(curve_id_chr, levels = levels(df_boot[[unit_used]])), !!.trt := factor(trt_chr, levels = levels(df_boot[[.trt]])))
      
    if (!is.null(object$vars$covariates)) {
      for (cov in object$vars$covariates) {
        cov_vals <- object$genotype_info[[cov]][match(new_curve_grid[[.trt]], object$genotype_info[[.trt]])]
        new_curve_grid[[cov]] <- cov_vals
      }
    }
      
    excl_c <- character()
    if(any(grepl(sprintf("s\\(%s\\)", unit_used), names(m_gam$smooth), fixed = TRUE))) excl_c <- sprintf("s(%s)", unit_used)
    
    new_curve_grid$y_hat <- as.numeric(mgcv::predict.gam(m_gam, newdata = new_curve_grid, type = "response", exclude = excl_c))
    boot_safri_list <- vector("list", n_boot)
    
    for (b in seq_len(n_boot)) {
      ids_b <- sample(ids0, replace = TRUE)
      w_draw <- tibble::tibble(id = factor(ids_b, levels = ids0)) |> dplyr::count(id, name = "w")
      pred_b <- new_curve_grid |> dplyr::mutate(id = factor(.data[[id_col]], levels = ids0)) |>
        dplyr::left_join(w_draw, by = "id") |> dplyr::filter(!is.na(.data[["w"]])) |>
        dplyr::group_by(.data[[.trt]], .data[[.time]]) |>
        dplyr::summarise(mu = stats::weighted.mean(.data[["y_hat"]], w = .data[["w"]], na.rm = TRUE), .groups = "drop")
        
      mat_b <- pred_b |> tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["mu"]]) |> dplyr::arrange(.data[[.time]])
      X_b <- as.matrix(mat_b[, setdiff(names(mat_b), .time), drop = FALSE])
      
      safri_b <- numeric(ncol(X_b))
      names(safri_b) <- colnames(X_b)
      
      for (j in seq_len(ncol(X_b))) {
        trt_j <- colnames(X_b)[j]
        g_val <- as.character(geno_groups$adjust_group[geno_groups[[.trt]] == trt_j])
        ref_b_nm <- if (reference_within_group && adjust_by != "All") {
           if (length(reference) > 1 && !is.null(names(reference))) reference[as.character(g_val)] else reference[1]
        } else reference[1]
        
        if (!ref_b_nm %in% colnames(X_b)) {
           safri_b[trt_j] <- NA
           next
        }
        ref_b <- X_b[, ref_b_nm]
        
        fri_b <- if(method == "positive_area") sum(pmax(ref_b - X_b[, j], 0)) * dt_grid else sqrt(sum((ref_b - X_b[, j])^2) * dt_grid)
        safri_b[trt_j] <- fri_b
        
        if (!is.null(instability)) {
          nfi_val <- res_df$nFI[res_df[[.trt]] == trt_j]
          if(length(nfi_val)>0) safri_b[trt_j] <- fri_b / (1 + lambda * nfi_val[1])
        }
      }
      
      safri_b <- safri_b[!is.na(safri_b)]
      
      if (adjust_by != "All") {
         cov_b <- geno_groups$adjust_group[match(names(safri_b), geno_groups[[.trt]])]
         boot_safri_list[[b]] <- tibble::tibble(!!.trt := names(safri_b), SAFRI_b = safri_b, !!adjust_by := cov_b, iter = b)
      } else {
         boot_safri_list[[b]] <- tibble::tibble(!!.trt := names(safri_b), SAFRI_b = safri_b, iter = b)
      }
    }
    
    boot_safri_df <- dplyr::bind_rows(boot_safri_list)
    boot_ci_tbl <- boot_safri_df |> dplyr::group_by(.data[[.trt]]) |>
      dplyr::summarise(SAFRI_lower = stats::quantile(.data[["SAFRI_b"]], probs = (1-ci_level)/2, na.rm = TRUE),
                       SAFRI_upper = stats::quantile(.data[["SAFRI_b"]], probs = 1-(1-ci_level)/2, na.rm = TRUE), .groups = "drop")
    res_df <- res_df |> dplyr::left_join(boot_ci_tbl, by = .trt)
    
    assign_boot_classes <- function(sub_res, sub_boot) {
      trts_ordered <- sub_res[[.trt]]
      n_t <- length(trts_ordered)
      mat_boot <- sub_boot |> tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["SAFRI_b"]], id_cols = "iter")
      
      group_assign <- numeric(n_t)
      current_class <- 1
      for (i in seq_len(n_t)) {
        if (group_assign[i] > 0) next
        group_assign[i] <- current_class
        for (j in seq_len(n_t)[-(1:i)]) {
          if (group_assign[j] > 0) next
          diff_ij <- mat_boot[[trts_ordered[i]]] - mat_boot[[trts_ordered[j]]]
          ci_diff <- stats::quantile(diff_ij, probs = c((1-ci_level)/2, 1-(1-ci_level)/2), na.rm = TRUE)
          if (ci_diff[1] <= 0 && ci_diff[2] >= 0) group_assign[j] <- current_class
        }
        current_class <- current_class + 1
      }
      group_assign <- pmin(group_assign, n_groups)
      sub_res$bootstrap_group <- class_names[group_assign]
      sub_res$resistance_class <- factor(sub_res$bootstrap_group, levels = class_names)
      sub_res
    }

    if (adjust_by != "All" && group_within_adjust_by) {
      res_list <- split(res_df, res_df[[adjust_by]])
      boot_list <- split(boot_safri_df, boot_safri_df[[adjust_by]])
      for (g in names(res_list)) res_list[[g]] <- assign_boot_classes(res_list[[g]], boot_list[[g]])
      res_df <- dplyr::bind_rows(res_list) |> dplyr::arrange(dplyr::desc(.data[["SAFRI"]]))
    } else {
      res_df <- assign_boot_classes(res_df, boot_safri_df)
    }
    
    res_df$n_boot <- n_boot
    boot_out <- list(boot_distributions = boot_safri_df)
  }

  out <- list(table = res_df, settings = list(lambda = lambda, method = method, scale_nfi = scale_nfi, group_method = group_method), bootstrap = boot_out)
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
