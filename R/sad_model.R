# Validation of standard area diagrams (SADs): one-stage hierarchical model
# and functional comparison of error curves
# sad_fit(), sad_curve(), sad_regions(), sad_raters(), sad_images(), sad_simulate()

utils::globalVariables(c(
  "z", "bias", "lower", "upper", "expected_abs_error", "contrast", "significant",
  "from", "to", "bias_at_center", "slope", "diff", "difficulty", "se"
))

#' Fit the one-stage hierarchical model for SAD validation
#'
#' @description
#' Fits a single mixed model to all severity estimates of a SAD validation
#' study, instead of first summarising each rater by Lin's statistics. The model
#' describes how the mean estimate, the variation among raters and the residual
#' variation change along the severity range in each condition.
#'
#' @param x An object returned by [sad_data()].
#' @param structure Random-effects structure for raters. `"auto"` (default)
#'   tries `"full"`, `"noim"`, `"sep"` and `"base"` in this order and keeps the
#'   first that converges with a positive definite Hessian and finite standard
#'   errors. See Details.
#' @param dispersion Model for the residual variance: `"quad"` (default), a
#'   quadratic function of severity in each condition, which allows larger
#'   errors at intermediate severity; or `"lin"`, a log-linear function.
#' @param family `"gaussian"` (estimates in %) or `"ordbeta"` (ordered beta
#'   regression on the 0-1 scale, which handles the bounds at 0 and 100% and
#'   accepts exact zeros and ones). `"ordbeta"` is useful as a sensitivity
#'   analysis.
#' @param center,scale Severity is entered as `z = (actual - center) / scale`.
#'   The defaults (10% and 10 percentage points) give intercepts that refer to
#'   10% severity and help convergence.
#'
#' @details
#' The mean model is `estimate ~ z * condition`, so each condition has its own
#' intercept (bias at `center`) and slope (calibration). Random-effects
#' structures:
#' \describe{
#'   \item{full}{`(0 + method + method:z | rater) + (1 | image) + (1 | image:condition)`:
#'     a random intercept and slope per rater and method, all correlated, an
#'     image effect (intrinsic difficulty) and an image x condition effect.}
#'   \item{noim}{as `full` without the image x condition effect.}
#'   \item{sep}{`(0 + method | rater) + (0 + method:z | rater) + (1 | image)`:
#'     intercepts and slopes uncorrelated.}
#'   \item{base}{`(1 + z | rater) + (1 | image)`: rater effects common to all methods.}
#' }
#' The dispersion formula is `~ condition * (z + I(z^2))` or `~ condition * z`.
#' Models are fitted with `glmmTMB::glmmTMB()`.
#'
#' @return An object of class `sad_fit`: a list with the `model`, the
#'   `structure` used, the `data` (with `z`), `center`, `scale`, `family`,
#'   `dispersion`, `conditions` and `methods`.
#' @seealso [sad_curve()], [sad_raters()], [sad_images()], [sad_compare()]
#' @examples
#' \donttest{
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#'   fit <- sad_fit(x)
#'   fit
#' }
#' }
#' @export
sad_fit <- function(x, structure = c("auto", "full", "noim", "sep", "base"),
                    dispersion = c("quad", "lin"), family = c("gaussian", "ordbeta"),
                    center = 10, scale = 10) {
  .sad_need("glmmTMB", "sad_fit")
  if (!inherits(x, "sad_data")) stop("`x` must be created with sad_data().", call. = FALSE)
  structure <- match.arg(structure); dispersion <- match.arg(dispersion); family <- match.arg(family)
  d <- tibble::as_tibble(unclass(x))
  d$z <- (d$actual - center) / scale
  d$condition <- droplevels(d$condition); d$method <- droplevels(d$method)
  if (nlevels(d$condition) < 2) stop("At least two conditions are needed.", call. = FALSE)
  d$y <- if (family == "ordbeta") pmin(pmax(d$estimate / 100, 0), 1) else d$estimate
  one_method <- nlevels(d$method) == 1
  re <- list(
    full = "(0 + method + method:z | rater) + (1 | image) + (1 | image:condition)",
    noim = "(0 + method + method:z | rater) + (1 | image)",
    sep  = "(0 + method | rater) + (0 + method:z | rater) + (1 | image)",
    base = "(1 + z | rater) + (1 | image)")
  if (one_method) re$full <- re$noim <- re$sep <- NULL
  disp_f <- if (dispersion == "quad") ~ condition * (z + I(z^2)) else ~ condition * z
  fam <- if (family == "ordbeta") glmmTMB::ordbeta() else stats::gaussian()
  ok <- function(m) {
    if (inherits(m, "try-error") || !isTRUE(m$sdr$pdHess)) return(FALSE)
    all(is.finite(sqrt(diag(stats::vcov(m)$cond))))
  }
  try_specs <- if (structure == "auto") names(re) else structure
  if (!all(try_specs %in% names(re))) stop("Structure not available with a single method.", call. = FALSE)
  for (s in try_specs) {
    f <- stats::as.formula(paste("y ~ z * condition +", re[[s]]))
    m <- suppressWarnings(try(glmmTMB::glmmTMB(f, dispformula = disp_f, family = fam, data = d),
                              silent = TRUE))
    if (ok(m)) {
      out <- list(model = m, structure = s, data = d, center = center, scale = scale,
                  family = family, dispersion = dispersion,
                  conditions = levels(d$condition), methods = levels(d$method))
      class(out) <- "sad_fit"
      return(out)
    }
    if (structure == "auto") message("Structure '", s, "' did not converge; trying the next one.")
  }
  stop("The model did not converge.", call. = FALSE)
}

#' @export
print.sad_fit <- function(x, ...) {
  cat("One-stage model for SAD validation\n")
  cat(sprintf("  family: %s | dispersion: %s | rater structure: %s\n",
              x$family, x$dispersion, x$structure))
  cat(sprintf("  severity centred at %g%%, in units of %g percentage points\n", x$center, x$scale))
  cat(sprintf("  %d raters, %d images, %d estimates; conditions: %s\n",
              length(unique(x$data$rater)), length(unique(x$data$image)), nrow(x$data),
              paste(x$conditions, collapse = ", ")))
  cat(sprintf("  AIC = %.1f\n", stats::AIC(x$model)))
  invisible(x)
}

#' @export
summary.sad_fit <- function(object, ...) summary(object$model)

# Covariance matrix of rater intercepts and slopes by method
.rater_cov <- function(fit) {
  vc <- glmmTMB::VarCorr(fit$model)$cond; lv <- fit$methods
  nm <- c(paste0("method", lv), paste0("method", lv, ":z"))
  S <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
  if (fit$structure == "base") {
    r <- vc$rater
    for (a in lv) for (b in lv) {
      S[paste0("method", a), paste0("method", b)] <- r[1, 1]
      S[paste0("method", a, ":z"), paste0("method", b, ":z")] <- r[2, 2]
      S[paste0("method", a), paste0("method", b, ":z")] <- r[1, 2]
      S[paste0("method", b, ":z"), paste0("method", a)] <- r[1, 2]
    }
  } else {
    for (r in vc[grepl("^rater", names(vc))]) S[rownames(r), colnames(r)] <- r
  }
  S
}

# Pairs of conditions to compare: each condition against the reference within
# the same group of raters
.sad_pairs <- function(d, reference) {
  cm <- unique(as.data.frame(d)[c("condition", "method", "group")])
  cm[] <- lapply(cm, as.character)
  conds <- levels(droplevels(d$condition))
  ref <- if (is.null(reference)) conds[1] else reference
  if (ref %in% conds) {
    p <- data.frame(a = setdiff(conds, ref), b = ref)
  } else if (ref %in% cm$method) {
    refs <- cm[cm$method == ref, ]
    p <- merge(cm[cm$method != ref, c("condition", "group")], refs[c("condition", "group")],
               by = "group", suffixes = c(".a", ".b"))
    p <- data.frame(a = p$condition.a, b = p$condition.b)
    p <- p[order(match(p$a, conds)), ]
  } else stop("`reference` must be a condition or a method.", call. = FALSE)
  if (!nrow(p)) stop("No pairs of conditions to compare.", call. = FALSE)
  p
}

#' Bias, variance components and expected error along the severity range
#'
#' @description
#' Functional part of the one-stage SAD validation. Describes, for each
#' condition, how the bias, the variation among raters, the image variation, the
#' residual variation and the expected absolute error change along the severity
#' range, and compares the bias curves between conditions with simultaneous
#' confidence bands, so that the severity ranges where a SAD helps or hurts are
#' identified.
#'
#' @param fit An object returned by [sad_fit()].
#' @param type `"linear"` (default) uses the fitted `sad_fit` model, so bias is
#'   linear in severity within each condition. `"smooth"` fits a functional
#'   generalised additive mixed model with `mgcv::bam()`: a smooth error curve
#'   per condition, a random smooth curve per rater and condition and a random
#'   image effect; bias can then take any shape.
#' @param reference Reference for the differences: a condition (e.g.
#'   `"UN_Old"`) or a method (e.g. `"UN"`). With a method, every other
#'   condition is compared with the reference method in the same group of
#'   raters. Defaults to the first condition.
#' @param pairs Optional data frame with columns `a` and `b`; differences
#'   `a - b` are computed for each row. Overrides `reference`.
#' @param n Number of grid points along the severity range.
#' @param level Confidence level of the bands.
#' @param k Basis dimension of the condition smooths (`type = "smooth"`).
#' @param nsim Number of draws used for the simultaneous bands.
#' @param seed Random seed for the simultaneous bands. The user's random number
#'   stream is restored afterwards.
#'
#' @details
#' With `type = "linear"`, the variance among raters at severity \eqn{z} is
#' \eqn{\sigma^2_{b0} + 2z\sigma_{b0b1} + z^2\sigma^2_{b1}}, from the rater
#' intercept and slope covariance of the method. The total standard deviation
#' \eqn{\sigma_T} adds image, image x condition and residual variances, and the
#' expected absolute error is the mean of a folded normal distribution,
#' \eqn{\sigma_T\sqrt{2/\pi}\exp(-b^2/2\sigma_T^2) + b(1 - 2\Phi(-b/\sigma_T))},
#' where \eqn{b} is the bias. Differences are evaluated only where the severity
#' ranges of the two conditions overlap. Simultaneous bands use the
#' distribution of the maximum absolute t statistic over the grid, from draws
#' of the fixed effects.
#'
#' @return An object of class `sad_curve`: a list with `curves` (one row per
#'   condition and grid point; for `"linear"` also `sd_rater`, `sd_image`,
#'   `sd_resid`, `sd_total` and `expected_abs_error`), `differences` (with
#'   `significant`), `type` and `level`. For `"smooth"` also `gam` and `edf`.
#' @seealso [sad_regions()], [plot.sad_curve()]
#' @examples
#' \donttest{
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#'   fit <- sad_fit(x)
#'   cv <- sad_curve(fit, reference = "UN")
#'   sad_regions(cv)
#'   plot(cv, "difference")
#' }
#' }
#' @export
sad_curve <- function(fit, type = c("linear", "smooth"), reference = NULL, pairs = NULL,
                      n = 100, level = 0.95, k = 5, nsim = 5000, seed = 1) {
  if (!inherits(fit, "sad_fit")) stop("`fit` must be created with sad_fit().", call. = FALSE)
  .sad_need("MASS", "sad_curve")
  type <- match.arg(type)
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  if (is.null(pairs)) pairs <- .sad_pairs(fit$data, reference)
  pairs <- data.frame(a = as.character(pairs$a), b = as.character(pairs$b))
  if (type == "smooth") return(.sad_curve_smooth(fit, pairs, n, level, k, nsim, seed))
  if (fit$family != "gaussian") stop("type = 'linear' needs family = 'gaussian' in sad_fit().", call. = FALSE)
  m <- fit$model; d <- fit$data
  cmeth <- unique(as.data.frame(d)[c("condition", "method")])
  rng <- dplyr::summarise(dplyr::group_by(d, condition), lo = min(actual), hi = max(actual))
  grid <- dplyr::bind_rows(lapply(seq_len(nrow(rng)), function(i)
    tibble::tibble(condition = as.character(rng$condition[i]),
                   actual = seq(rng$lo[i], rng$hi[i], length.out = n))))
  grid$method <- as.character(cmeth$method[match(grid$condition, cmeth$condition)])
  grid$condition <- factor(grid$condition, levels = fit$conditions)
  grid$method <- factor(grid$method, levels = fit$methods)
  grid$z <- (grid$actual - fit$center) / fit$scale
  grid$rater <- NA; grid$image <- NA
  beta <- glmmTMB::fixef(m)$cond; V <- stats::vcov(m)$cond
  Xfor <- function(cond, z) {
    nd <- data.frame(z = z, condition = factor(cond, levels = fit$conditions))
    stats::model.matrix(~ z * condition, data = nd)[, names(beta), drop = FALSE]
  }
  X <- Xfor(as.character(grid$condition), grid$z)
  fitv <- drop(X %*% beta); se <- sqrt(rowSums((X %*% V) * X))
  vc <- glmmTMB::VarCorr(m)$cond
  v_img <- vc$image[1, 1] + if (!is.null(vc$`image:condition`)) vc$`image:condition`[1, 1] else 0
  S <- .rater_cov(fit)
  sigma <- sqrt(stats::predict(m, newdata = grid, type = "disp", re.form = NA,
                               allow.new.levels = TRUE))
  i <- paste0("method", grid$method); j <- paste0(i, ":z")
  v_r <- S[cbind(i, i)] + 2 * grid$z * S[cbind(i, j)] + grid$z^2 * S[cbind(j, j)]
  b <- fitv - grid$actual
  sd_tot <- sqrt(v_r + v_img + sigma^2)
  fold <- sd_tot * sqrt(2 / pi) * exp(-b^2 / (2 * sd_tot^2)) + b * (1 - 2 * stats::pnorm(-b / sd_tot))
  zq <- stats::qnorm(1 - (1 - level) / 2)
  curves <- tibble::tibble(condition = grid$condition, method = grid$method, actual = grid$actual,
                           bias = b, lower = b - zq * se, upper = b + zq * se,
                           sd_rater = sqrt(v_r), sd_image = sqrt(v_img), sd_resid = sigma,
                           sd_total = sd_tot, expected_abs_error = fold)
  set.seed(seed)
  bu <- MASS::mvrnorm(nsim, rep(0, length(beta)), V)
  diffs <- dplyr::bind_rows(lapply(seq_len(nrow(pairs)), function(r) {
    a <- pairs$a[r]; bb <- pairs$b[r]
    lo <- max(rng$lo[rng$condition %in% c(a, bb)]); hi <- min(rng$hi[rng$condition %in% c(a, bb)])
    act <- seq(lo, hi, length.out = n); zz <- (act - fit$center) / fit$scale
    Xd <- Xfor(a, zz) - Xfor(bb, zz)
    .sad_band(Xd, beta, V, bu, level, paste(a, "-", bb), act)
  }))
  out <- list(curves = curves, differences = diffs, type = "linear", level = level)
  class(out) <- "sad_curve"
  out
}

.sad_band <- function(Xd, beta, V, bu, level, label, act) {
  est <- drop(Xd %*% beta); s <- sqrt(rowSums((Xd %*% V) * Xd))
  s[s < 1e-12] <- 1e-12
  crit <- stats::quantile(apply(abs((Xd %*% t(bu)) / s), 2, max), level, names = FALSE)
  tibble::tibble(contrast = label, actual = act, diff = est,
                 lower = est - crit * s, upper = est + crit * s,
                 significant = (est - crit * s) > 0 | (est + crit * s) < 0)
}

# Functional GAMM: smooth error curve per condition, random curve per rater
.sad_curve_smooth <- function(fit, pairs, n, level, k, nsim, seed) {
  d <- fit$data
  d$err <- d$estimate - d$actual
  d$rc <- interaction(d$rater, d$condition, drop = TRUE)
  g <- mgcv::bam(err ~ condition + s(actual, by = condition, k = k) +
                   s(actual, rc, bs = "fs", k = 5, m = 1) + s(image, bs = "re"),
                 data = d, method = "fREML", discrete = TRUE)
  rng <- dplyr::summarise(dplyr::group_by(d, condition), lo = min(actual), hi = max(actual))
  bcoef <- stats::coef(g); V <- stats::vcov(g, unconditional = TRUE)
  re_cols <- grep("s\\(actual,rc\\)|s\\(image\\)", names(bcoef))
  Xfor <- function(cond, act) {
    nd <- data.frame(actual = act, condition = factor(cond, levels = levels(d$condition)),
                     rc = d$rc[1], image = d$image[1])
    X <- stats::predict(g, newdata = nd, type = "lpmatrix"); X[, re_cols] <- 0; X
  }
  set.seed(seed); bu <- MASS::mvrnorm(nsim, rep(0, length(bcoef)), V)
  curves <- dplyr::bind_rows(lapply(seq_len(nrow(rng)), function(i) {
    act <- seq(rng$lo[i], rng$hi[i], length.out = n)
    bb <- .sad_band(Xfor(as.character(rng$condition[i]), act), bcoef, V, bu, level, "", act)
    tibble::tibble(condition = rng$condition[i], actual = act, bias = bb$diff,
                   lower = bb$lower, upper = bb$upper)
  }))
  diffs <- dplyr::bind_rows(lapply(seq_len(nrow(pairs)), function(r) {
    a <- pairs$a[r]; b <- pairs$b[r]
    lo <- max(rng$lo[rng$condition %in% c(a, b)]); hi <- min(rng$hi[rng$condition %in% c(a, b)])
    act <- seq(lo, hi, length.out = n)
    .sad_band(Xfor(a, act) - Xfor(b, act), bcoef, V, bu, level, paste(a, "-", b), act)
  }))
  st <- summary(g)$s.table
  out <- list(curves = curves, differences = diffs, type = "smooth", level = level,
              gam = g, edf = st[grepl("^s\\(actual\\):", rownames(st)), "edf"])
  class(out) <- "sad_curve"
  out
}

#' Severity ranges where conditions differ
#'
#' @description
#' Summarises the differences of a [sad_curve()] object into the severity
#' ranges where the simultaneous band excludes zero.
#'
#' @param curve An object returned by [sad_curve()].
#' @return A tibble with `contrast`, `from`, `to` (actual severity, %) and
#'   `direction` of the difference in bias (`"a > b"` or `"a < b"`).
#' @export
sad_regions <- function(curve) {
  dd <- curve$differences
  res <- lapply(split(dd, factor(dd$contrast, levels = unique(dd$contrast))), function(z) {
    run <- cumsum(c(1, diff(z$significant) != 0))
    out <- lapply(split(z, run), function(w) {
      if (!w$significant[1]) return(NULL)
      tibble::tibble(contrast = w$contrast[1], from = min(w$actual), to = max(w$actual),
                     direction = ifelse(mean(w$diff) < 0, "a < b", "a > b"))
    })
    dplyr::bind_rows(out)
  })
  out <- dplyr::bind_rows(res)
  if (!nrow(out)) out <- tibble::tibble(contrast = character(), from = numeric(),
                                        to = numeric(), direction = character())
  out
}

#' Plot curves from sad_curve()
#'
#' @param x An object returned by [sad_curve()].
#' @param what `"bias"` (bias curves with pointwise intervals),
#'   `"expected_abs_error"` (only for `type = "linear"`) or `"difference"`
#'   (differences with simultaneous bands; shaded ranges are significant).
#' @param ... Not used.
#' @return A ggplot object.
#' @export
plot.sad_curve <- function(x, what = c("bias", "expected_abs_error", "difference"), ...) {
  what <- match.arg(what)
  if (what == "difference") {
    reg <- sad_regions(x)
    p <- ggplot2::ggplot(x$differences, ggplot2::aes(actual, diff))
    if (nrow(reg)) p <- p + ggplot2::geom_rect(
      data = reg, ggplot2::aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey88")
    return(p + ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .3) +
      ggplot2::geom_line(linewidth = 1) + ggplot2::facet_wrap(~ contrast) +
      ggplot2::labs(x = "Actual severity (%)", y = "Difference in bias (percentage points)",
                    subtitle = sprintf("%d%% simultaneous band; shaded = significant difference",
                                       round(100 * x$level))) + ggplot2::theme_bw())
  }
  if (what == "expected_abs_error") {
    if (x$type != "linear") stop("Expected absolute error is available only for type = 'linear'.", call. = FALSE)
    return(ggplot2::ggplot(x$curves, ggplot2::aes(actual, expected_abs_error, colour = condition)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(x = "Actual severity (%)", y = "Expected absolute error (percentage points)",
                    colour = NULL) + ggplot2::theme_bw())
  }
  ggplot2::ggplot(x$curves, ggplot2::aes(actual, bias, colour = condition, fill = condition)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .2, colour = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = "Actual severity (%)", y = "Bias (percentage points)",
                  colour = NULL, fill = NULL) + ggplot2::theme_bw()
}

#' Bias and calibration of each rater
#'
#' @description
#' Combines fixed and random effects of a [sad_fit()] model into the bias of
#' each rater at the centring severity and the slope of the rater's estimates
#' on actual severity (1 = perfect calibration), per condition.
#'
#' @param fit An object returned by [sad_fit()] with `family = "gaussian"`.
#' @param x An object returned by `sad_raters()`.
#' @param ... Not used.
#' @return `sad_raters()` returns a tibble with `rater`, `condition`,
#'   `method`, `bias_at_center` (percentage points) and `slope`.
#'   `plot_sad_raters()` returns a ggplot object, with lines linking the same
#'   rater across conditions.
#' @examples
#' \donttest{
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   x <- sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad")
#'   rt <- sad_raters(sad_fit(x))
#'   plot_sad_raters(rt)
#' }
#' }
#' @export
sad_raters <- function(fit) {
  if (fit$family != "gaussian") stop("sad_raters() needs family = 'gaussian'.", call. = FALSE)
  m <- fit$model; b <- glmmTMB::fixef(m)$cond
  re <- glmmTMB::ranef(m)$cond
  R <- lapply(re[grepl("^rater", names(re))], function(z) {
    z <- as.data.frame(z); z$rater <- rownames(z); z
  })
  R <- Reduce(function(a, bb) merge(a, bb, by = "rater", all = TRUE), R)
  d <- unique(as.data.frame(fit$data)[c("rater", "condition", "method")])
  d$rater <- as.character(d$rater)
  d <- merge(d, R, by = "rater", all.x = TRUE, sort = FALSE)
  get_b <- function(nm) ifelse(nm %in% names(b), b[nm], 0)
  if (fit$structure == "base") {
    int_re <- d[["(Intercept)"]]; slope_re <- d[["z"]]
  } else {
    int_re <- d[cbind(seq_len(nrow(d)), match(paste0("method", d$method), names(d)))]
    slope_re <- d[cbind(seq_len(nrow(d)), match(paste0("method", d$method, ":z"), names(d)))]
    int_re <- as.numeric(int_re); slope_re <- as.numeric(slope_re)
  }
  int_fix <- b["(Intercept)"] + get_b(paste0("condition", d$condition))
  slope_fix <- b["z"] + get_b(paste0("z:condition", d$condition))
  out <- tibble::tibble(rater = d$rater, condition = factor(d$condition, levels = fit$conditions),
                        method = factor(d$method, levels = fit$methods),
                        bias_at_center = unname(int_fix + int_re - fit$center),
                        slope = unname((slope_fix + slope_re) / fit$scale))
  out <- out[order(out$condition, suppressWarnings(as.numeric(out$rater)), out$rater), ]
  attr(out, "center") <- fit$center
  out
}

#' @rdname sad_raters
#' @export
plot_sad_raters <- function(x, ...) {
  center <- if (is.null(attr(x, "center"))) 10 else attr(x, "center")
  ggplot2::ggplot(x, ggplot2::aes(bias_at_center, slope, colour = condition)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(group = rater), colour = "grey75") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = sprintf("Rater bias at %g%% severity (percentage points)", center),
                  y = "Slope (calibration; ideal = 1)", colour = NULL) +
    ggplot2::theme_bw()
}

#' Intrinsic difficulty of each image
#'
#' @description
#' Extracts the image random effects of a [sad_fit()] model: the shift, common
#' to all raters and conditions, of the estimates of each image relative to
#' what its actual severity predicts. Large absolute values flag images that
#' are hard to assess, or whose actual severity may be wrong.
#'
#' @param fit An object returned by [sad_fit()].
#' @return A tibble with `image`, `actual`, `difficulty` (percentage points for
#'   the gaussian family) and its conditional standard error `se`, ordered by
#'   absolute difficulty.
#' @export
sad_images <- function(fit) {
  re <- glmmTMB::ranef(fit$model, condVar = TRUE)$cond$image
  se <- sqrt(attr(re, "condVar")[1, 1, ])
  act <- unique(as.data.frame(fit$data)[c("image", "actual")])
  act <- act[!duplicated(act$image), ]
  out <- tibble::tibble(image = rownames(re), difficulty = re[, 1], se = se)
  out$actual <- act$actual[match(out$image, as.character(act$image))]
  out <- out[order(-abs(out$difficulty)), c("image", "actual", "difficulty", "se")]
  out
}

#' Simulate power for planning a SAD validation study
#'
#' @description
#' Simulates paired studies (each rater assesses the same images without and
#' with a SAD) and compares the power of the conventional two-stage approach,
#' a paired t test on \eqn{\rho_c}, with the likelihood ratio test of the
#' one-stage model for any effect of the SAD on mean or dispersion (4 df).
#'
#' @param n_rater,n_image Vectors with the numbers of raters and images; all
#'   combinations are simulated.
#' @param nsim Number of simulated studies per scenario.
#' @param shift Change in bias at 10% severity caused by the SAD (percentage points).
#' @param slope_change Change in slope caused by the SAD.
#' @param precision Multiplier of the rater and residual standard deviations
#'   with the SAD (< 1 means the SAD makes estimates more precise).
#' @param bias0,slope0 Bias at 10% severity and slope without the SAD.
#' @param sd_rater,sd_slope Standard deviations of rater intercepts and slopes.
#' @param sd_image Standard deviation of image difficulty.
#' @param sigma0 Residual standard deviation at 10% severity without the SAD.
#' @param rho_rater Correlation of a rater's effects with and without the SAD.
#' @param severity Function of `n` that draws the actual severities.
#' @param cores Number of cores (`parallel::mclapply()`; ignored on Windows).
#' @param seed Random seed; the user's random number stream is restored.
#'
#' @details
#' Default parameters are typical values estimated from ten published SAD
#' validation studies. Estimates are truncated to 0-100%. The one-stage model
#' uses the `"sep"` structure. With fewer than about 20 raters the likelihood
#' ratio test can be liberal; a parametric bootstrap is then safer.
#'
#' @return A tibble with one row per scenario and the rejection rates at 5%
#'   (`power_two_stage`, `power_one_stage`) and the number of converged fits.
#' @examples
#' \donttest{
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   sad_simulate(n_rater = 10, n_image = 20, nsim = 5)
#' }
#' }
#' @export
sad_simulate <- function(n_rater = c(10, 20, 40), n_image = c(20, 40), nsim = 200,
                         shift = -3, slope_change = 0, precision = 1,
                         bias0 = 3, slope0 = 1.10, sd_rater = 5, sd_slope = 0.18,
                         sd_image = 3.5, sigma0 = 8, rho_rater = 0.5,
                         severity = function(n) stats::runif(n, 0.5, 80), cores = 1, seed = 1) {
  .sad_need("glmmTMB", "sad_simulate")
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  one <- function(nr, ni) {
    a0 <- severity(ni); zc <- (a0 - 10) / 10
    v <- stats::rnorm(ni, 0, sd_image)
    b0u <- stats::rnorm(nr, 0, sd_rater); b1u <- stats::rnorm(nr, 0, sd_slope)
    b0a <- precision * (rho_rater * b0u + sqrt(1 - rho_rater^2) * stats::rnorm(nr, 0, sd_rater))
    b1a <- precision * (rho_rater * b1u + sqrt(1 - rho_rater^2) * stats::rnorm(nr, 0, sd_slope))
    d <- expand.grid(s = 0:1, i = seq_len(ni), r = seq_len(nr))
    d$actual <- a0[d$i]; d$z <- zc[d$i]
    sg <- sigma0 * exp(0.1 * d$z) * ifelse(d$s == 1, precision, 1)
    mu <- 10 + bias0 + slope0 * 10 * d$z + d$s * (shift + slope_change * 10 * d$z) +
      ifelse(d$s == 1, b0a[d$r] + b1a[d$r] * 10 * d$z, b0u[d$r] + b1u[d$r] * 10 * d$z) + v[d$i]
    d$estimate <- pmin(pmax(mu + stats::rnorm(nrow(d), 0, sg), 0), 100)
    d$method <- factor(ifelse(d$s == 1, "aided", "unaided"), levels = c("unaided", "aided"))
    d$rater <- factor(d$r); d$image <- factor(d$i)
    rc <- tapply(seq_len(nrow(d)), list(d$rater, d$method),
                 function(ix) .lin_stats(d$estimate[ix], d$actual[ix])$rho_c)
    p_two <- stats::t.test(rc[, "aided"], rc[, "unaided"], paired = TRUE)$p.value
    re <- "(0 + method | rater) + (0 + method:z | rater) + (1 | image)"
    f1 <- stats::as.formula(paste("estimate ~ z * method +", re))
    f0 <- stats::as.formula(paste("estimate ~ z +", re))
    m1 <- suppressWarnings(try(glmmTMB::glmmTMB(f1, dispformula = ~ method * z, data = d), silent = TRUE))
    m0 <- suppressWarnings(try(glmmTMB::glmmTMB(f0, dispformula = ~ z, data = d), silent = TRUE))
    p_one <- if (inherits(m1, "try-error") || inherits(m0, "try-error")) NA else
      stats::pchisq(max(2 * (m0$fit$objective - m1$fit$objective), 0), 4, lower.tail = FALSE)
    c(p_two = p_two, p_one = p_one)
  }
  sc <- expand.grid(n_image = n_image, n_rater = n_rater)
  set.seed(seed)
  use_mc <- cores > 1 && .Platform$OS.type != "windows"
  res <- lapply(seq_len(nrow(sc)), function(j) {
    f <- function(i) one(sc$n_rater[j], sc$n_image[j])
    r <- if (use_mc) parallel::mclapply(seq_len(nsim), f, mc.cores = cores) else lapply(seq_len(nsim), f)
    r <- do.call(rbind, r)
    tibble::tibble(n_rater = sc$n_rater[j], n_image = sc$n_image[j],
                   power_two_stage = mean(r[, "p_two"] < 0.05, na.rm = TRUE),
                   power_one_stage = mean(r[, "p_one"] < 0.05, na.rm = TRUE),
                   n_ok = sum(!is.na(r[, "p_one"])))
  })
  dplyr::bind_rows(res)
}
