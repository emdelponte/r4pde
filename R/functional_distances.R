#' Compute pairwise functional distances and cluster curves
#'
#' @description
#' Computes L2 functional distances among treatment mean trajectories evaluated on a
#' common time grid from a \code{functional_curves} object. Also computes curve-level
#' distances, performs hierarchical clustering, and optionally performs permutation testing.
#'
#' @param object An object of class \code{functional_curves}.
#' @param cluster_k Number of clusters used to cut the hierarchical tree.
#' @param hc_method Clustering linkage method passed to \code{stats::hclust()}.
#' @param test_factor Optional a priori grouping used for a distance-based permutation test.
#' @param n_perm Number of permutations for the distance-based test.
#' @param test_mode Character; permutation test mode.
#' @param min_strata_for_pairwise Minimum number of strata levels required for pairwise tests.
#' @param perm_unit Optional character naming the unit at which group labels are permuted.
#' @param perm_strata Optional character naming a factor defining restricted permutation strata.
#' @param bootstrap Logical; if \code{TRUE}, compute bootstrap envelopes and distance summaries.
#' @param boot_B Integer number of bootstrap replicates.
#' @param boot_seed Integer seed for bootstrap resampling.
#' @param boot_ci Length-2 numeric vector of quantiles for bootstrap intervals.
#' @param show_progress Logical; whether to show progress.
#' @param curve_level Logical; whether to compute curve-level distance matrix.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"functional_distances"} containing:
#' \itemize{
#'   \item \code{distance}: treatment-by-treatment functional distance matrix;
#'   \item \code{distance_table}: pairwise distances as a data frame;
#'   \item \code{hc}: \code{hclust} object;
#'   \item \code{clusters}: treatment cluster assignments;
#'   \item \code{curve_distance}: curve-by-curve distance matrix;
#'   \item \code{test}: permutation test results;
#'   \item \code{bootstrap}: bootstrap results if requested.
#' }
#'
#' @examples
#' \dontrun{
#' fd <- functional_distances(fc, cluster_k = 4)
#' print(fd)
#' plot_dendrogram(fd)
#' plot_curves(fd)
#' }
#' @export
functional_distances <- function(
    object,
    cluster_k = 4,
    hc_method = "ward.D2",
    test_factor = NULL,
    n_perm = 999,
    test_mode = c("auto","none","global","global_pairwise"),
    min_strata_for_pairwise = 6,
    perm_unit = NULL,
    perm_strata = NULL,
    bootstrap = FALSE,
    boot_B = 399,
    boot_seed = 1,
    boot_ci = c(0.025, 0.975),
    show_progress = TRUE,
    curve_level = NULL,
    ...
) {
  if (!inherits(object, "functional_curves")) {
    stop("`object` must be of class 'functional_curves'")
  }

  test_mode <- match.arg(test_mode)

  df <- object$observed_data
  m_gam <- object$gam
  t_grid <- object$grid
  pred_trt <- object$curves
  trt_score <- object$scores

  .time <- object$vars$time
  .trt <- object$vars$treatment
  .env <- object$vars$environment
  .blk <- object$vars$block
  unit_used <- object$vars$unit

  trt_levels0 <- levels(df[[.trt]])
  unit_levels0 <- levels(df[[unit_used]])
  env_levels0 <- if (!is.null(.env)) levels(df[[.env]]) else NULL
  blk_levels0 <- if (!is.null(.blk)) levels(df[[.blk]]) else NULL

  # ---- functional distance + clustering (treatment means) ----
  mat <- pred_trt |>
    dplyr::select(.data[[.time]], .data[[.trt]], .data[["mu"]]) |>
    tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["mu"]]) |>
    dplyr::arrange(.data[[.time]])

  dt_grid <- mean(diff(mat[[.time]]))
  X <- as.matrix(mat[, setdiff(names(mat), .time), drop = FALSE])

  # L2 distance via trapezoidal (dt_grid scaling is proportional to trapezoidal for regular grid)
  D <- as.matrix(stats::dist(t(X), method = "euclidean")) * sqrt(dt_grid)
  
  # Distance table
  D_dist <- stats::as.dist(D)
  D_mat <- as.matrix(D_dist)
  dist_df <- as.data.frame(as.table(D_mat))
  names(dist_df) <- c("group1", "group2", "distance")
  dist_df <- dist_df[as.character(dist_df$group1) < as.character(dist_df$group2), ]
  rownames(dist_df) <- NULL

  hc <- stats::hclust(stats::as.dist(D), method = hc_method)
  
  cl_raw <- stats::cutree(hc, k = cluster_k)

  cluster_tbl <- tibble::tibble(
    !!.trt := names(cl_raw),
    cluster_raw = as.integer(cl_raw)
  ) |>
    dplyr::left_join(trt_score, by = .trt) |>
    dplyr::group_by(.data[["cluster_raw"]]) |>
    dplyr::summarise(cluster_mean = mean(.data[["mean_mu"]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data[["cluster_mean"]]) |>
    dplyr::mutate(cluster = dplyr::row_number()) |>
    dplyr::select(.data[["cluster_raw"]], .data[["cluster"]])

  trt_cluster <- tibble::tibble(
    !!.trt := names(cl_raw),
    cluster_raw = as.integer(cl_raw)
  ) |>
    dplyr::left_join(cluster_tbl, by = "cluster_raw") |>
    dplyr::left_join(trt_score, by = .trt)

  # ---- Curve-level fitted trajectories (exclude unit RE) + distance matrix ----
  D_curve <- NULL
  new_curve_grid <- NULL
  
  if (is.null(curve_level)) {
    curve_level <- !is.null(test_factor) || isTRUE(bootstrap)
  }

  if (isTRUE(curve_level)) {
    col_name1 <- function(x){
      if (is.null(x)) return(NULL)
      if (is.character(x) && length(x) == 1L) return(x)
      if (rlang::is_quosure(x)) return(rlang::as_name(rlang::quo_get_expr(x)))
      if (rlang::is_symbol(x))  return(rlang::as_name(x))
      if (is.list(x) && length(x) == 1L) return(col_name1(x[[1]]))
      stop("Could not coerce column spec to a single name. Got: ", paste(class(x), collapse = "/"))
    }

    env_col <- col_name1(.env)
    blk_col <- col_name1(.blk)
    ab_cols <- c(env_col, blk_col)
    ab_cols <- ab_cols[!is.na(ab_cols) & nzchar(ab_cols)]

    curve_meta <- df |>
      dplyr::distinct(.data[[unit_used]], .data[[.trt]], dplyr::across(dplyr::all_of(ab_cols))) |>
      dplyr::mutate(
        curve_id_chr = as.character(.data[[unit_used]]),
        trt_chr      = as.character(.data[[.trt]])
      )

    dup <- curve_meta |>
      dplyr::count(curve_id_chr) |>
      dplyr::filter(n > 1)
    if (nrow(dup) > 0) {
      stop("Each curve (unit) must map to one treatment/block/env. Duplicates for: ",
           paste(dup$curve_id_chr, collapse = ", "))
    }

    new_curve_grid <- tidyr::crossing(
      curve_id_chr = curve_meta$curve_id_chr,
      !!.time := t_grid
    ) |>
      dplyr::left_join(curve_meta, by = "curve_id_chr") |>
      dplyr::mutate(
        !!unit_used := factor(curve_id_chr, levels = unit_levels0),
        !!.trt      := factor(trt_chr, levels = trt_levels0)
      )

    if (!is.null(.env) && .env %in% names(new_curve_grid)) {
      new_curve_grid[[.env]] <- factor(new_curve_grid[[.env]], levels = env_levels0)
    }
    if (!is.null(.blk) && .blk %in% names(new_curve_grid)) {
      new_curve_grid[[.blk]] <- factor(new_curve_grid[[.blk]], levels = blk_levels0)
    }

    excl_c <- character()
    if(any(grepl(sprintf("s\\(%s\\)", unit_used), names(m_gam$smooth), fixed = TRUE))) {
       excl_c <- c(excl_c, sprintf("s(%s)", unit_used))
    }

    new_curve_grid$y_hat <- as.numeric(
      mgcv::predict.gam(
        m_gam,
        newdata = new_curve_grid,
        type = "response",
        exclude = excl_c
      )
    )

    Yhat <- new_curve_grid |>
      dplyr::select(curve_id_chr, !!.time, y_hat) |>
      tidyr::pivot_wider(names_from = !!.time, values_from = y_hat) |>
      dplyr::arrange(curve_id_chr)

    curve_ids <- Yhat$curve_id_chr
    Ymat <- as.matrix(Yhat[, -1, drop = FALSE])

    D_curve <- as.matrix(stats::dist(Ymat, method = "euclidean")) * sqrt(dt_grid)
    rownames(D_curve) <- curve_ids
    colnames(D_curve) <- curve_ids

    D_curve <- D_curve / stats::median(D_curve[D_curve > 0])
  }

  # ---- Permutation test ----
  test_res <- NULL
  if (!is.null(test_factor) && test_mode != "none") {
    if (is.null(perm_strata)) perm_strata <- .blk
    perm_unit_eff <- if (is.null(perm_unit)) unit_used else perm_unit

    if (!perm_unit_eff %in% names(df)) {
      stop("`perm_unit` (or default unit) not found in data: ", perm_unit_eff)
    }

    permute_within_strata <- function(g, strata) {
      gperm <- g
      for (s in unique(strata)) {
        idx <- which(strata == s)
        gperm[idx] <- sample(g[idx], length(idx), replace = FALSE)
      }
      gperm
    }

    permanova_F <- function(Dmat, grp) {
      n <- nrow(Dmat)
      grp <- droplevels(factor(grp))
      k <- nlevels(grp)
      if (n < 3 || k < 2) return(list(F = NA_real_, df1 = NA_integer_, df2 = NA_integer_))
      df1 <- k - 1
      df2 <- n - k
      if (df2 <= 0) return(list(F = NA_real_, df1 = df1, df2 = df2))

      A <- -0.5 * (Dmat^2)
      H <- diag(n) - matrix(1, n, n)/n
      G <- H %*% A %*% H
      Xg <- stats::model.matrix(~ grp)
      P  <- Xg %*% solve(t(Xg) %*% Xg) %*% t(Xg)

      SS_between <- sum(diag(P %*% G))
      SS_total   <- sum(diag(G))
      SS_within  <- SS_total - SS_between
      Fstat <- (SS_between/df1) / (SS_within/df2)
      list(F = as.numeric(Fstat), df1 = df1, df2 = df2)
    }

    if (is.character(test_factor) && length(test_factor) == 1L) {
      tf_col <- test_factor
      if (!tf_col %in% names(df)) stop("`test_factor` not found in data: ", tf_col)

      map_u <- df |>
        dplyr::distinct(.data[[perm_unit_eff]], .data[[tf_col]], .keep_all = FALSE) |>
        dplyr::filter(!is.na(.data[[tf_col]])) |>
        dplyr::mutate(u_chr = as.character(.data[[perm_unit_eff]]))

      dup <- map_u |> dplyr::count(u_chr) |> dplyr::filter(n > 1)
      if (nrow(dup) > 0) {
        stop("`test_factor` must be unique per `perm_unit`. Duplicates for: ",
             paste(dup$u_chr, collapse = ", "))
      }

      g_u <- map_u[[tf_col]]
      names(g_u) <- map_u$u_chr

    } else {
      g_u <- test_factor
      if (is.null(names(g_u))) stop("If `test_factor` is a vector, it must be named by `perm_unit` ids.")
    }

    map_curve_to_u <- df |>
      dplyr::distinct(.data[[unit_used]], .data[[perm_unit_eff]], .keep_all = FALSE) |>
      dplyr::mutate(curve_id_chr = as.character(.data[[unit_used]]),
                    u_chr = as.character(.data[[perm_unit_eff]]))

    dup2 <- map_curve_to_u |> dplyr::count(curve_id_chr) |> dplyr::filter(n > 1)
    if (nrow(dup2) > 0) stop("Each curve must map to one perm_unit. Duplicates for: ", paste(dup2$curve_id_chr, collapse = ", "))

    curve_ids2 <- rownames(D_curve)
    u_curve <- map_curve_to_u$u_chr
    names(u_curve) <- map_curve_to_u$curve_id_chr
    u_curve <- u_curve[curve_ids2]

    g_curve <- g_u[u_curve]
    g_curve <- droplevels(factor(g_curve))

    strata_curve <- NULL
    if (!is.null(perm_strata)) {
      if (!perm_strata %in% names(df)) stop("`perm_strata` not found in data: ", perm_strata)
      map_curve_to_s <- df |>
        dplyr::distinct(.data[[unit_used]], .data[[perm_strata]], .keep_all = FALSE) |>
        dplyr::mutate(curve_id_chr = as.character(.data[[unit_used]]))
      
      strata_curve <- map_curve_to_s[[perm_strata]]
      names(strata_curve) <- map_curve_to_s$curve_id_chr
      strata_curve <- droplevels(factor(strata_curve[curve_ids2]))
    }

    pairwise_allowed <- TRUE
    if (test_mode == "auto") {
      if (!is.null(strata_curve)) {
        n_strata <- nlevels(strata_curve)
        if (n_strata < min_strata_for_pairwise) pairwise_allowed <- FALSE
      }
    }
    if (test_mode == "global") pairwise_allowed <- FALSE
    if (test_mode == "global_pairwise") pairwise_allowed <- TRUE

    obs <- permanova_F(D_curve, g_curve)
    set.seed(1)
    
    if (show_progress) {
      pb_glob <- progress::progress_bar$new(
        format = "  Global permutation test [:bar] :percent eta: :eta",
        total = n_perm, clear = FALSE, width = 60, force = TRUE, show_after = 0)
    }
    
    F_perm <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      if (show_progress) pb_glob$tick()
      gperm <- if (is.null(strata_curve)) {
        sample(g_curve, replace = FALSE)
      } else {
        permute_within_strata(g_curve, strata_curve)
      }
      F_perm[i] <- permanova_F(D_curve, gperm)$F
    }
    
    p_global <- (1 + sum(F_perm >= obs$F, na.rm = TRUE)) / (n_perm + 1)

    pairwise <- NULL
    lev <- levels(g_curve)
    if (pairwise_allowed && length(lev) > 2) {
      pairs <- utils::combn(lev, 2, simplify = FALSE)
      
      if (show_progress) {
        pb_pair <- progress::progress_bar$new(
          format = "  Pairwise permutation tests [:bar] :percent eta: :eta",
          total = length(pairs) * n_perm, clear = FALSE, width = 60, force = TRUE, show_after = 0)
      }

      pairwise <- purrr::map_dfr(pairs, function(pp) {
        keep2 <- g_curve %in% pp
        D2 <- D_curve[keep2, keep2, drop = FALSE]
        g2 <- droplevels(g_curve[keep2])
        s2 <- if (is.null(strata_curve)) NULL else droplevels(strata_curve[keep2])

        ob2 <- permanova_F(D2, g2)
        
        Fp2 <- numeric(n_perm)
        for (j in seq_len(n_perm)) {
          if (show_progress) pb_pair$tick()
          gp <- if (is.null(s2)) {
            sample(g2, replace = FALSE)
          } else {
            permute_within_strata(g2, s2)
          }
          Fp2[j] <- permanova_F(D2, gp)$F
        }
        
        p2 <- (1 + sum(Fp2 >= ob2$F, na.rm = TRUE)) / (n_perm + 1)
        tibble::tibble(group1 = pp[1], group2 = pp[2], F = ob2$F, p = p2)
      }) |>
        dplyr::mutate(p_adj = stats::p.adjust(p, method = "holm"))
    }

    test_res <- list(
      type = "permutation_distance_test",
      mode = test_mode,
      perm_unit = perm_unit_eff,
      perm_strata = perm_strata,
      n_perm = n_perm,
      groups = g_curve,
      strata = strata_curve,
      global = c(F = obs$F, p = p_global, df1 = obs$df1, df2 = obs$df2),
      pairwise = pairwise,
      perm_F = F_perm
    )
  }

  # ---- Bootstrap ----
  boot_res <- NULL
  if (isTRUE(bootstrap)) {
    set.seed(boot_seed)
    if(!is.null(.blk)) {
      ids0 <- blk_levels0
      id_col <- .blk
      boot_unit_label <- "block"
    } else {
      ids0 <- unit_levels0
      id_col <- unit_used
      boot_unit_label <- "unit"
    }

    if (show_progress) {
      pb_boot <- progress::progress_bar$new(
        format = "  Bootstrapping predictions [:bar] :percent eta: :eta",
        total = boot_B, clear = FALSE, width = 60, force = TRUE, show_after = 0)
    }

    boot_pred_list <- vector("list", boot_B)

    for (b in seq_len(boot_B)) {
      if (show_progress) pb_boot$tick()
      ids_b <- sample(ids0, replace = TRUE)

      w_draw <- tibble::tibble(id = factor(ids_b, levels = ids0)) |>
        dplyr::count(id, name = "w")

      grid_b <- new_curve_grid |>
        dplyr::mutate(id = factor(.data[[id_col]], levels = ids0)) |>
        dplyr::left_join(w_draw, by = "id") |>
        dplyr::filter(!is.na(.data[["w"]]))

      boot_pred_list[[b]] <- grid_b |>
        dplyr::group_by(.data[[.trt]], .data[[.time]]) |>
        dplyr::summarise(
          mu = stats::weighted.mean(.data[["y_hat"]], w = .data[["w"]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(.iter = b)
    }

    boot_pred <- dplyr::bind_rows(boot_pred_list)
    trt_order <- object$scores[[.trt]]
    env_tbl <- boot_pred |>
      dplyr::group_by(.data[[.trt]], .data[[.time]]) |>
      dplyr::summarise(
        mu_boot_mean = mean(.data[["mu"]], na.rm = TRUE),
        mu_lo = stats::quantile(.data[["mu"]], probs = boot_ci[1], na.rm = TRUE),
        mu_hi = stats::quantile(.data[["mu"]], probs = boot_ci[2], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(!!.trt := factor(.data[[.trt]], levels = trt_order))

    boot_dist_one_iter <- function(df_iter) {
      mat_i <- df_iter |>
        dplyr::select(.data[[.time]], .data[[.trt]], .data[["mu"]]) |>
        tidyr::pivot_wider(names_from = .data[[.trt]], values_from = .data[["mu"]]) |>
        dplyr::arrange(.data[[.time]])

      X_i <- as.matrix(mat_i[, setdiff(names(mat_i), .time), drop = FALSE])
      coln <- colnames(X_i)
      if (ncol(X_i) < 2) return(NULL)

      pairs <- utils::combn(coln, 2, simplify = FALSE)
      purrr::map_dfr(pairs, function(pp){
        d <- sqrt(sum((X_i[, pp[1]] - X_i[, pp[2]])^2) * dt_grid)
        tibble::tibble(group1 = pp[1], group2 = pp[2], dist = as.numeric(d))
      })
    }

    if (show_progress) {
      pb_dist <- progress::progress_bar$new(
        format = "  Calculating bootstrap distances [:bar] :percent eta: :eta",
        total = boot_B, clear = FALSE, width = 60, force = TRUE, show_after = 0)
    }

    dist_boot <- boot_pred |>
      dplyr::group_by(.iter) |>
      dplyr::group_modify(~ {
        if (show_progress) pb_dist$tick()
        boot_dist_one_iter(.x)
      }) |>
      dplyr::ungroup()

    dist_summary <- dist_boot |>
      dplyr::group_by(.data[["group1"]], .data[["group2"]]) |>
      dplyr::summarise(
        dist_mean = mean(.data[["dist"]], na.rm = TRUE),
        dist_med  = stats::median(.data[["dist"]], na.rm = TRUE),
        dist_lo   = stats::quantile(.data[["dist"]], probs = boot_ci[1], na.rm = TRUE),
        dist_hi   = stats::quantile(.data[["dist"]], probs = boot_ci[2], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(.data[["dist_mean"]]))

    boot_res <- list(
      method = "bootstrap_on_predicted_curves",
      unit = boot_unit_label,
      id_col = id_col,
      B = boot_B,
      seed = boot_seed,
      ci = boot_ci,
      envelopes = env_tbl,
      distances = dist_summary,
      distances_boot = dist_boot
    )
  }

  out <- list(
    distance_matrix = D,
    distance_table = dist_df,
    hc = hc,
    clusters = trt_cluster,
    curve_distance = D_curve,
    test = test_res,
    bootstrap = boot_res,
    functional_curves = object
  )
  class(out) <- "functional_distances"
  out
}

#' Print functional_distances
#' @export
print.functional_distances <- function(x, ...) {
  cat("A functional_distances object\n")
  cat("Dimensions of distance matrix:", nrow(x$distance_matrix), "x", ncol(x$distance_matrix), "\n")
  cat("Number of clusters:", length(unique(x$clusters$cluster)), "\n")
  invisible(x)
}
