# Validation of standard area diagrams (SADs): conventional two-stage workflow
# sad_data(), sad_ccc(), sad_compare(), sad_reliability(), sad_plots()

utils::globalVariables(c(
  "rater", "image", "actual", "estimate", "method", "group", "condition",
  "error", "emmean", "SE", "lower.CL", "upper.CL", ".group", "stat", "cell",
  "N", "Cb", "u", "r", "rho_c", "v", "value", "rho_other", "gain", "panel",
  "other", "k", "n"
))

.sad_need <- function(pkgs, what) {
  rlang::check_installed(pkgs, reason = paste0("to use `", what, "()`."))
}

#' Prepare severity estimates for SAD validation
#'
#' @description
#' Organises visual severity estimates made with and without a standard area
#' diagram (SAD), or with two or more SADs, into the long format used by the
#' `sad_*` family, and checks the data for common problems.
#'
#' @param data A data frame in long or wide format.
#' @param rater,image,actual Names of the columns with the rater identifier,
#'   the image (leaf, fruit, ...) identifier and the actual severity (%).
#' @param estimate Name of the column with the estimates (long format).
#' @param method Name of the column with the assessment method (long format),
#'   for example `"unaided"` and `"aided"`.
#' @param wide Character vector with the estimate columns in wide format. Each
#'   column becomes one level of `method`. If the vector is named, the names are
#'   used as the method labels, e.g. `c(UN = "unaided", SAD = "aided")`.
#' @param group Optional name of a column that splits raters into groups, for
#'   example which SAD each rater used. When given, `condition` is
#'   `method` and `group` pasted together.
#' @param sep Separator used to build `condition` from `method` and `group`.
#'
#' @details
#' Each combination of method and group is a *condition*. The function warns
#' when an image has more than one actual severity value, when estimates fall
#' outside 0 to 100%, when a rater x image x condition combination is repeated,
#' and when a rater lacks one of the methods. The design is labelled `"paired"`
#' when raters assessed the images under more than one condition (the usual
#' unaided then aided study), and `"independent"` otherwise.
#'
#' @return A tibble of class `sad_data` with columns `rater`, `image`,
#'   `actual`, `estimate`, `method`, `group`, `condition` and `error`
#'   (estimate minus actual, in percentage points), and attribute `design`.
#'
#' @seealso [sad_ccc()], [sad_compare()], [sad_fit()]
#' @examples
#' x <- sad_data(SADRustSoybean, rater = "rater", image = "leaf",
#'               wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#' x
#' @export
sad_data <- function(data, rater = "rater", image = "leaf", actual = "actual",
                     estimate = "estimate", method = "method", wide = NULL,
                     group = NULL, sep = "_") {
  d <- tibble::as_tibble(data)
  needed <- c(rater, image, actual, group, if (is.null(wide)) c(estimate, method) else wide)
  miss <- setdiff(needed, names(d))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!is.null(wide)) {
    labs <- if (is.null(names(wide))) wide else ifelse(names(wide) == "", wide, names(wide))
    d <- tidyr::pivot_longer(d, dplyr::all_of(unname(wide)), names_to = ".method",
                             values_to = ".estimate")
    d$.method <- factor(labs[match(d$.method, wide)], levels = labs)
    method <- ".method"; estimate <- ".estimate"
  }
  mth <- d[[method]]
  out <- tibble::tibble(
    rater    = factor(d[[rater]]),
    image    = factor(d[[image]]),
    actual   = as.numeric(d[[actual]]),
    estimate = as.numeric(d[[estimate]]),
    method   = if (is.factor(mth)) droplevels(mth) else factor(mth, levels = unique(mth)),
    group    = if (is.null(group)) factor("all") else factor(d[[group]], levels = unique(d[[group]])))
  out <- out[!is.na(out$estimate) & !is.na(out$actual), ]
  if (is.null(group)) {
    out$condition <- out$method
  } else {
    lv <- as.vector(t(outer(levels(out$method), levels(out$group), paste, sep = sep)))
    out$condition <- factor(paste(out$method, out$group, sep = sep), levels = lv)
    out$condition <- droplevels(out$condition)
  }
  out$error <- out$estimate - out$actual

  msg <- character()
  inc <- tapply(out$actual, out$image, function(a) length(unique(a)))
  if (any(inc > 1, na.rm = TRUE))
    msg <- c(msg, sprintf("%d images have more than one actual severity value.", sum(inc > 1, na.rm = TRUE)))
  if (any(out$estimate < 0 | out$estimate > 100))
    msg <- c(msg, "Some estimates are outside 0-100%.")
  dup <- dplyr::count(out, rater, image, condition)
  if (any(dup$n > 1))
    msg <- c(msg, sprintf("%d rater x image x condition combinations are repeated.", sum(dup$n > 1)))
  tab <- table(out$rater, out$method)
  if (ncol(tab) > 1) {
    lack <- rownames(tab)[rowSums(tab == 0) > 0 & rowSums(tab > 0) > 1]
    if (length(lack)) msg <- c(msg, sprintf("Raters missing a method: %s.", paste(lack, collapse = ", ")))
  }
  if (length(msg)) warning(paste(msg, collapse = "\n"), call. = FALSE)

  nc <- tapply(out$condition, out$rater, function(z) length(unique(z)))
  attr(out, "design") <- if (all(nc > 1)) "paired" else "independent"
  class(out) <- c("sad_data", class(out))
  out
}

#' @export
print.sad_data <- function(x, ...) {
  cat(sprintf("SAD validation data (%s design): %d raters, %d images, %d estimates\n",
              attr(x, "design"), length(unique(x$rater)), length(unique(x$image)), nrow(x)))
  cat("Conditions:", paste(levels(droplevels(x$condition)), collapse = ", "), "\n")
  print(tibble::as_tibble(unclass(x)), ...)
  invisible(x)
}

# Lin's concordance statistics, population variances as in Lin (1989).
# u is negative when estimates are below the actual values.
.lin_stats <- function(est, act) {
  ok <- is.finite(est) & is.finite(act); est <- est[ok]; act <- act[ok]
  n <- length(est)
  s1 <- sqrt(sum((est - mean(est))^2) / n); s2 <- sqrt(sum((act - mean(act))^2) / n)
  r <- stats::cor(est, act); v <- s1 / s2; u <- (mean(est) - mean(act)) / sqrt(s1 * s2)
  Cb <- 2 / (v + 1 / v + u^2)
  f <- stats::lm(est ~ act)
  tibble::tibble(n = n, rho_c = r * Cb, r = r, Cb = Cb, v = v, u = u,
                 intercept = unname(stats::coef(f)[1]), slope = unname(stats::coef(f)[2]))
}

#' Lin's concordance statistics per rater and condition
#'
#' @description
#' First stage of the conventional SAD validation: computes Lin's concordance
#' correlation coefficient and its components for each rater under each
#' condition.
#'
#' @param x An object returned by [sad_data()].
#'
#' @details
#' \eqn{\rho_c = r \times C_b}, where \eqn{r} is Pearson's correlation
#' (precision) and \eqn{C_b = 2 / (v + 1/v + u^2)} is the bias correction
#' factor (accuracy). \eqn{v = \sigma_{est}/\sigma_{act}} is the scale shift
#' and \eqn{u = (\mu_{est} - \mu_{act}) / \sqrt{\sigma_{est}\sigma_{act}}} is
#' the location shift, negative when the rater underestimates. Note that
#' `epiR::epi.ccc()` reports the location shift with the opposite sign.
#' The intercept and slope of the regression of estimates on actual values are
#' also returned.
#'
#' @return A tibble with one row per rater and condition and columns `n`,
#'   `rho_c`, `r`, `Cb`, `v`, `u`, `intercept` and `slope`.
#' @references Lin, L. I. (1989). A concordance correlation coefficient to
#'   evaluate reproducibility. Biometrics 45:255-268.
#' @seealso [sad_compare()]
#' @examples
#' x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#' head(sad_ccc(x))
#' @export
sad_ccc <- function(x) {
  if (!inherits(x, "sad_data")) stop("`x` must be created with sad_data().", call. = FALSE)
  keys <- unique(as.data.frame(x)[c("rater", "group", "method", "condition")])
  res <- lapply(seq_len(nrow(keys)), function(i) {
    s <- x$rater == keys$rater[i] & x$condition == keys$condition[i]
    .lin_stats(x$estimate[s], x$actual[s])
  })
  out <- dplyr::bind_cols(tibble::as_tibble(keys), dplyr::bind_rows(res))
  out <- out[order(out$condition, out$rater), ]
  out$condition <- droplevels(out$condition)
  out
}

#' Compare concordance statistics between conditions
#'
#' @description
#' Second stage of the conventional SAD validation: compares the mean of each
#' Lin statistic between conditions (for example unaided versus aided, or two
#' or more SADs) and returns the summary table with letters of a Tukey test,
#' as in Table 2 of Franceschi et al. (2020).
#'
#' @param x An object returned by [sad_data()].
#' @param ccc Statistics per rater, by default `sad_ccc(x)`.
#' @param stats Statistics to compare.
#' @param adjust Multiplicity adjustment passed to `emmeans::contrast()`.
#' @param conditions Optional character vector to select and order conditions.
#' @param alpha Significance level for the letters.
#' @param digits Digits used in the formatted table.
#'
#' @details
#' When the same raters appear in more than one condition, each statistic is
#' analysed with a linear mixed model with rater as a random intercept
#' (`lmerTest::lmer(y ~ condition + (1 | rater))`), so the comparison is paired.
#' Otherwise an ordinary linear model is used. Letters are assigned so that `a`
#' is the condition with the highest mean; means that share a letter do not
#' differ at level `alpha`.
#'
#' @return An object of class `sad_compare`: a list with `table` (formatted
#'   table), `estimates` (means, standard errors and confidence limits),
#'   `pairwise` (all pairwise contrasts), `models` and `ccc`.
#' @references Franceschi, V. T., Alves, K. S., Mazaro, S. M., Godoy, C. V.,
#'   Duarte, H. S. S., and Del Ponte, E. M. (2020). A new standard area diagram
#'   set for assessment of severity of soybean rust improves accuracy of
#'   estimates and optimizes resource use. Plant Pathology 69:495-505.
#' @seealso [sad_ccc()], [sad_reliability()], [sad_fit()]
#' @examples
#' \donttest{
#' x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#' sad_compare(x)
#' }
#' @export
sad_compare <- function(x, ccc = sad_ccc(x),
                        stats = c("u", "v", "Cb", "r", "rho_c"),
                        adjust = "tukey", conditions = NULL, alpha = 0.05, digits = 2) {
  .sad_need(c("emmeans", "multcompView", "lmerTest"), "sad_compare")
  if (!is.null(conditions)) {
    ccc <- ccc[ccc$condition %in% conditions, ]
    ccc$condition <- factor(ccc$condition, levels = conditions)
  } else ccc$condition <- droplevels(ccc$condition)
  if (nlevels(ccc$condition) < 2) stop("At least two conditions are needed.", call. = FALSE)
  paired <- any(table(ccc$rater) > 1)
  res <- lapply(stats, function(s) {
    d <- ccc; d$y <- d[[s]]
    m <- if (paired) suppressMessages(lmerTest::lmer(y ~ condition + (1 | rater), data = d)) else
      stats::lm(y ~ condition, data = d)
    em <- emmeans::emmeans(m, ~ condition)
    pw <- as.data.frame(emmeans::contrast(em, "pairwise", adjust = adjust))
    emd <- as.data.frame(em)
    lv <- .sad_letters(pw, emd, alpha)
    tab <- tibble::tibble(condition = as.character(emd$condition), stat = s,
                          mean = emd$emmean, se = emd$SE,
                          lower = emd[[grep("lower", names(emd))[1]]],
                          upper = emd[[grep("upper", names(emd))[1]]],
                          letter = unname(lv[as.character(emd$condition)]))
    pw$stat <- s
    list(tab = tab, pw = tibble::as_tibble(pw), model = m)
  })
  names(res) <- stats
  long <- dplyr::bind_rows(lapply(res, `[[`, "tab"))
  long$condition <- factor(long$condition, levels = levels(ccc$condition))
  long$cell <- paste(formatC(long$mean, format = "f", digits = digits), long$letter)
  table2 <- tidyr::pivot_wider(long[c("condition", "stat", "cell")],
                               names_from = "stat", values_from = "cell")
  nn <- as.data.frame(table(condition = ccc$condition), responseName = "N")
  table2 <- dplyr::left_join(table2, nn, by = "condition")
  table2 <- table2[c("condition", "N", stats)]
  out <- list(table = table2, estimates = long[setdiff(names(long), "cell")],
              pairwise = dplyr::bind_rows(lapply(res, `[[`, "pw")),
              models = lapply(res, `[[`, "model"), ccc = ccc,
              design = if (paired) "linear mixed model, rater random" else "linear model",
              adjust = adjust, alpha = alpha)
  class(out) <- "sad_compare"
  out
}

# Compact letter display with "a" for the highest mean
.sad_letters <- function(pw, emd, alpha) {
  lv_names <- as.character(emd$condition)
  if (length(lv_names) == 1) return(stats::setNames("a", lv_names))
  # contrast labels are "A - B"; rebuild them from the levels so that level
  # names with spaces or dashes are handled
  pairs <- utils::combn(lv_names, 2)
  lab <- paste(pairs[1, ], pairs[2, ], sep = " - ")
  lab_par <- paste0("(", pairs[1, ], ") - (", pairs[2, ], ")")
  p <- pw$p.value[match(lab, as.character(pw$contrast))]
  p2 <- pw$p.value[match(lab_par, as.character(pw$contrast))]
  p[is.na(p)] <- p2[is.na(p)]
  safe <- make.names(lv_names, unique = TRUE)
  pv <- stats::setNames(p, paste(safe[match(pairs[1, ], lv_names)],
                                 safe[match(pairs[2, ], lv_names)], sep = "-"))
  L <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
  L <- stats::setNames(L[safe], lv_names)
  ord <- lv_names[order(emd$emmean, decreasing = TRUE)]
  seen <- unique(unlist(strsplit(L[ord], "")))
  remap <- stats::setNames(letters[seq_along(seen)], seen)
  vapply(L, function(z) paste(sort(remap[strsplit(z, "")[[1]]]), collapse = ""), "")
}

#' @export
print.sad_compare <- function(x, ...) {
  cat("Lin's statistics by condition -", x$design, "\n")
  cat(sprintf("Means followed by the same letter do not differ (%s, alpha = %g).\n\n",
              x$adjust, x$alpha))
  print(as.data.frame(x$table), row.names = FALSE)
  invisible(x)
}

# ICC(A,1): two-way random effects, absolute agreement, single rater
# (McGraw and Wong 1996, case 2A), with confidence interval
.icc_a1 <- function(M, conf = 0.95) {
  M <- M[stats::complete.cases(M), , drop = FALSE]
  n <- nrow(M); k <- ncol(M); gm <- mean(M)
  SSR <- k * sum((rowMeans(M) - gm)^2); SSC <- n * sum((colMeans(M) - gm)^2)
  SSE <- sum((M - gm)^2) - SSR - SSC
  MSR <- SSR / (n - 1); MSC <- SSC / (k - 1); MSE <- SSE / ((n - 1) * (k - 1))
  icc <- (MSR - MSE) / (MSR + (k - 1) * MSE + k / n * (MSC - MSE))
  a <- k * icc / (n * (1 - icc)); b <- 1 + k * icc * (n - 1) / (n * (1 - icc))
  v <- (a * MSC + b * MSE)^2 / ((a * MSC)^2 / (k - 1) + (b * MSE)^2 / ((n - 1) * (k - 1)))
  al <- 1 - conf
  FL <- stats::qf(1 - al / 2, n - 1, v); FU <- stats::qf(1 - al / 2, v, n - 1)
  lb <- n * (MSR - FL * MSE) / (FL * (k * MSC + (k * n - k - n) * MSE) + n * MSR)
  ub <- n * (FU * MSR - MSE) / (k * MSC + (k * n - k - n) * MSE + n * FU * MSR)
  c(icc = icc, lower = lb, upper = ub, n_images = n, n_raters = k)
}

# Overall concordance correlation coefficient (Barnhart et al. 2002): the
# weighted average of pairwise CCCs, with sample variances as in epiR::epi.occc()
.occc <- function(M) {
  k <- ncol(M); mu <- colMeans(M); C <- stats::cov(M)
  num <- 0; den <- 0
  for (j in 1:(k - 1)) for (l in (j + 1):k) {
    num <- num + 2 * C[j, l]
    den <- den + C[j, j] + C[l, l] + (mu[j] - mu[l])^2
  }
  unname(num / den)
}

#' Inter-rater reliability by condition
#'
#' @description
#' Computes the intraclass correlation coefficient ICC(A,1) with its
#' confidence interval and the overall concordance correlation coefficient
#' (OCC) among raters for each condition, as in Table 3 of Franceschi et al.
#' (2020).
#'
#' @param x An object returned by [sad_data()].
#' @param conf Confidence level for the ICC interval.
#'
#' @details
#' ICC(A,1) is the two-way random effects, absolute agreement, single rater
#' coefficient of McGraw and Wong (1996). OCC is the overall concordance
#' correlation coefficient of Barnhart et al. (2002), computed as in
#' `epiR::epi.occc()`. Only raters that assessed
#' all images of a condition are used.
#'
#' @return A tibble with one row per condition and columns `ICC`,
#'   `ICC_lower`, `ICC_upper`, `OCC`, `n_raters` and `n_images`.
#' @references McGraw, K. O., and Wong, S. P. (1996). Forming inferences about
#'   some intraclass correlation coefficients. Psychological Methods 1:30-46.
#'
#'   Barnhart, H. X., Haber, M., and Song, J. (2002). Overall concordance
#'   correlation coefficient for evaluating agreement among multiple observers.
#'   Biometrics 58:1020-1027.
#' @examples
#' x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#' sad_reliability(x)
#' @export
sad_reliability <- function(x, conf = 0.95) {
  if (!inherits(x, "sad_data")) stop("`x` must be created with sad_data().", call. = FALSE)
  conds <- levels(droplevels(x$condition))
  res <- lapply(conds, function(cn) {
    d <- x[x$condition == cn, c("image", "rater", "estimate")]
    w <- tidyr::pivot_wider(as.data.frame(d), names_from = "rater", values_from = "estimate",
                            values_fn = mean)
    M <- as.matrix(w[, -1, drop = FALSE])
    M <- M[, colSums(is.na(M)) == 0, drop = FALSE]
    if (ncol(M) < 2) return(tibble::tibble(condition = cn, ICC = NA_real_, ICC_lower = NA_real_,
                                           ICC_upper = NA_real_, OCC = NA_real_,
                                           n_raters = ncol(M), n_images = nrow(M)))
    ic <- .icc_a1(M, conf)
    tibble::tibble(condition = cn, ICC = unname(ic["icc"]), ICC_lower = unname(ic["lower"]),
                   ICC_upper = unname(ic["upper"]), OCC = .occc(M),
                   n_raters = ncol(M), n_images = nrow(M))
  })
  out <- dplyr::bind_rows(res)
  out$condition <- factor(out$condition, levels = conds)
  out
}

#' Exploratory plots for SAD validation
#'
#' @description
#' Produces the exploratory figures of a conventional SAD validation study:
#' estimates against actual severity, Cb against r, boxplots of Lin's
#' statistics, errors along the severity range, error densities, the standard
#' deviation among raters along the severity range and the gain in accuracy
#' relative to unaided estimates.
#'
#' @param x An object returned by [sad_data()].
#' @param ccc Statistics per rater, by default `sad_ccc(x)`.
#' @param baseline Optional name of the method used as baseline (for example
#'   `"UN"`); when given, the plot of gain in \eqn{\rho_c} against baseline
#'   accuracy is added.
#' @param palette Optional named vector of colours, one per condition.
#'
#' @return A named list of ggplot objects: `scatter`, `cb_vs_r`, `box_ccc`,
#'   `error`, `error_density` (requires ggridges), `sd_vs_actual` and, with
#'   `baseline`, `gain_vs_baseline`.
#' @examples
#' x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#' p <- sad_plots(x, baseline = "UN")
#' p$scatter
#' @export
sad_plots <- function(x, ccc = sad_ccc(x), baseline = NULL, palette = NULL) {
  force(ccc)
  conds <- levels(droplevels(x$condition))
  base_pal <- c("#4D4D4D", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442")
  pal <- if (is.null(palette)) stats::setNames(rep_len(base_pal, length(conds)), conds) else palette
  th <- ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")
  x <- tibble::as_tibble(unclass(x)); x$condition <- droplevels(x$condition)
  p <- list()
  p$scatter <- ggplot2::ggplot(x, ggplot2::aes(actual, estimate, colour = condition)) +
    ggplot2::geom_abline(colour = "grey50", linewidth = 1) +
    ggplot2::geom_point(alpha = .4, size = 1) +
    ggplot2::facet_wrap(~ condition, nrow = 1) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::labs(x = "Actual severity (%)", y = "Estimated severity (%)") + th
  p$cb_vs_r <- ggplot2::ggplot(ccc, ggplot2::aes(r, Cb, colour = condition)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = "Pearson's r (precision)", y = "Bias correction factor (Cb)", colour = NULL) +
    ggplot2::theme_bw()
  lab_stats <- c(rho_c = "\u03c1c (agreement)", r = "r (precision)", Cb = "Cb (accuracy)",
                 v = "v (scale shift)", u = "u (location shift)")
  bl <- tidyr::pivot_longer(ccc, c("Cb", "u", "r", "rho_c", "v"), names_to = "stat")
  bl$stat <- factor(lab_stats[bl$stat], levels = lab_stats)
  p$box_ccc <- ggplot2::ggplot(bl, ggplot2::aes(condition, value, colour = condition)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = .15, height = 0, size = 1.2, alpha = .7) +
    ggplot2::facet_wrap(~ stat, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = NULL, y = "Value") + th +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  p$error <- ggplot2::ggplot(x, ggplot2::aes(actual, error, colour = condition)) +
    ggplot2::geom_hline(yintercept = 0) + ggplot2::geom_point(alpha = .3, size = 1) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, colour = "red", se = FALSE) +
    ggplot2::facet_wrap(~ condition) + ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = "Actual severity (%)", y = "Error (percentage points)") + th
  if (requireNamespace("ggridges", quietly = TRUE)) {
    p$error_density <- ggplot2::ggplot(x, ggplot2::aes(error, condition, fill = condition)) +
      ggridges::geom_density_ridges(alpha = .8, scale = 1.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = "Error (percentage points)", y = NULL) + th
  } else {
    p$error_density <- ggplot2::ggplot(x, ggplot2::aes(error, fill = condition)) +
      ggplot2::geom_density(alpha = .5) + ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = "Error (percentage points)", y = "Density", fill = NULL) + ggplot2::theme_bw()
  }
  sdd <- dplyr::summarise(dplyr::group_by(x, condition, image, actual),
                          sd = stats::sd(estimate), .groups = "drop")
  p$sd_vs_actual <- ggplot2::ggplot(sdd, ggplot2::aes(actual, sd, colour = condition)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
    ggplot2::facet_wrap(~ condition, nrow = 1) + ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = "Actual severity (%)", y = "SD among raters (percentage points)") + th
  if (!is.null(baseline)) {
    if (!baseline %in% levels(x$method)) stop("`baseline` must be one of the methods.", call. = FALSE)
    gl <- tidyr::pivot_wider(ccc[c("rater", "group", "method", "rho_c")],
                             names_from = "method", values_from = "rho_c")
    others <- setdiff(levels(x$method), baseline)
    gl <- tidyr::pivot_longer(gl, dplyr::all_of(others), names_to = "other", values_to = "rho_other")
    gl$base <- gl[[baseline]]
    gl <- gl[!is.na(gl$rho_other) & !is.na(gl$base), ]
    gl$gain <- gl$rho_other - gl$base
    gl$panel <- if (length(unique(gl$group)) > 1) paste(gl$other, gl$group, sep = "_") else gl$other
    p$gain_vs_baseline <- ggplot2::ggplot(gl, ggplot2::aes(.data$base, gain)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_point(size = 2.5, colour = "#007f80") +
      ggplot2::facet_wrap(~ panel) +
      ggplot2::labs(x = sprintf("\u03c1c without aid (%s)", baseline),
                    y = "Gain or loss in \u03c1c") + ggplot2::theme_bw()
  }
  p
}
