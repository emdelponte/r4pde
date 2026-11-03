#' Functional network meta-analysis of disease suppression (arm-based)
#'
#' @description
#' Fits an arm-based functional network model to disease progress data from
#' multiple trials (environments) in which treatments are unbalanced, and returns
#' the mean disease progress curve of each treatment (arm) and the disease
#' suppression curve of each treatment relative to a reference (e.g. an untreated
#' check), with pointwise and simultaneous bands from posterior simulation.
#'
#' @details
#' Every treatment, including the reference, is an arm. The model is a GAMM for
#' the plot-level response (proportion scale, logit link):
#'
#' \deqn{logit(y) = trt + s(time, by = trt) + s(time, env, bs = "fs") +
#'   s(time, env:trt, bs = "fs") + s(unit, bs = "re")}
#'
#' The environment factor-smooth absorbs the level, timing and shape of each
#' epidemic, which is shared by all arms of a trial; the treatment-by-environment
#' term (\code{ge_term}) is the heterogeneity of the network; the unit (plot)
#' random intercept accounts for repeated assessments of the same plot. Treatments
#' that never appear together in a trial are connected through the arms they share
#' with other trials (typically the reference).
#'
#' Arm curves are predicted for a typical environment: the linear predictor
#' averages the environment factor-smooth over all environments and the
#' treatment-by-environment deviations over the environments in which each
#' treatment was tested (as in \code{\link{functional_curves}}), and excludes the
#' unit effects. The suppression curve of treatment \eqn{j} is
#' \eqn{\mu_{ref}(t) - \mu_j(t)} (difference) or
#' \eqn{(\mu_{ref}(t) - \mu_j(t)) / \mu_{ref}(t)} (relative). Bands come from
#' \code{n_sim} draws of the coefficients from their approximate posterior
#' (\code{mgcv::rmvn(coef, Vp)}), so they are conditional on the smoothing
#' parameters and on the observed environments.
#'
#' @param data A data frame with one row per assessment unit (plot) and time.
#' @param time,response,treatment,environment Column names (character).
#' @param reference Name of the reference treatment (e.g. the untreated check).
#' @param block Optional block column; with \code{unit = NULL}, plots are
#'   identified as environment x treatment x block.
#' @param unit Optional column identifying the assessment unit (plot).
#' @param response_scale \code{"percent"} (default) or \code{"proportion"}.
#' @param ge_term Treatment-by-environment term: \code{"smooth"} (default),
#'   \code{"intercept"} or \code{"none"}.
#' @param k_trt,k_env,k_ge Basis dimensions of the treatment smooths, the
#'   environment factor-smooth and the treatment-by-environment factor-smooth.
#' @param unit_re Logical; include a unit (plot) random intercept. Default \code{TRUE}.
#' @param family \code{"quasibinomial"} (default) or \code{"betar"}.
#' @param contrast \code{"difference"} (default) or \code{"relative"}.
#' @param grid_n Number of time points in the prediction grid.
#' @param n_sim Number of posterior draws for the bands.
#' @param seed Seed for the posterior draws; the user's random stream is restored.
#' @param discrete,gamma Passed to \code{mgcv::bam()}.
#' @param keep_draws Logical; keep the posterior draws of the suppression curves.
#'
#' @return An object of class \code{"functional_network"}: a list with \code{gam},
#'   \code{arms} (arm curves with bands), \code{contrast} (suppression curves with
#'   bands; column \code{DSP}), \code{draws} (if \code{keep_draws}), \code{vars},
#'   \code{settings} (including the model formula) and \code{fit_time}.
#'
#' @seealso \code{\link{functional_curves}}, \code{\link{functional_suppression_profiles}}
#' @export
functional_network <- function(
    data,
    time,
    response,
    treatment,
    environment,
    reference,
    block = NULL,
    unit = NULL,
    response_scale = c("percent", "proportion"),
    ge_term = c("smooth", "intercept", "none"),
    k_trt = 8,
    k_env = 4,
    k_ge = 3,
    unit_re = TRUE,
    family = c("quasibinomial", "betar"),
    contrast = c("difference", "relative"),
    grid_n = 100,
    n_sim = 1000,
    seed = 1,
    discrete = TRUE,
    gamma = 1,
    keep_draws = TRUE
) {
  response_scale <- match.arg(response_scale)
  ge_term <- match.arg(ge_term)
  family <- match.arg(family)
  contrast <- match.arg(contrast)

  needed <- c(time, response, treatment, environment, block, unit)
  bad <- setdiff(needed, names(data))
  if (length(bad)) stop("Columns not found in `data`: ", paste(bad, collapse = ", "))

  df <- tibble::as_tibble(data)
  df <- df[stats::complete.cases(df[, needed]), , drop = FALSE]
  if (!reference %in% df[[treatment]]) {
    stop(sprintf("Reference treatment '%s' not found.", reference))
  }
  trt_levels <- c(reference, setdiff(sort(unique(as.character(df[[treatment]]))), reference))
  df$.trt <- factor(as.character(df[[treatment]]), levels = trt_levels)
  df$.env <- factor(as.character(df[[environment]]))
  df$.time <- as.numeric(df[[time]])
  y <- as.numeric(df[[response]])
  if (response_scale == "percent") y <- y / 100
  y <- pmin(pmax(y, 0), 1)
  if (family == "betar") {
    n <- length(y)
    y <- (y * (n - 1) + 0.5) / n
  }
  df$.y <- y
  df$.ge <- interaction(df$.env, df$.trt, drop = TRUE, sep = "|")
  if (!is.null(unit)) {
    df$.unit <- factor(as.character(df[[unit]]))
  } else if (!is.null(block)) {
    df$.unit <- interaction(df$.env, df$.trt, factor(df[[block]]), drop = TRUE, sep = "|")
  } else {
    df$.unit <- df$.ge
  }

  # Reference arm present in every environment?
  no_ref <- setdiff(levels(df$.env), unique(df$.env[df$.trt == reference]))
  if (length(no_ref)) {
    warning("Reference treatment missing in ", length(no_ref), " environment(s); ",
            "their arms are connected only through other shared treatments.", call. = FALSE)
  }

  terms <- c(".trt",
             sprintf("s(.time, by = .trt, k = %d)", k_trt),
             sprintf("s(.time, .env, bs = 'fs', k = %d, m = 1)", k_env))
  if (ge_term == "smooth") {
    terms <- c(terms, sprintf("s(.time, .ge, bs = 'fs', k = %d, m = 1)", k_ge))
  } else if (ge_term == "intercept") {
    terms <- c(terms, "s(.ge, bs = 're')")
  }
  unit_is_ge <- nlevels(df$.unit) == nlevels(df$.ge)
  if (unit_re && !unit_is_ge) terms <- c(terms, "s(.unit, bs = 're')")
  form <- stats::as.formula(paste(".y ~", paste(terms, collapse = " + ")))

  fam <- if (family == "betar") mgcv::betar(link = "logit") else stats::quasibinomial(link = "logit")
  t0 <- Sys.time()
  m <- mgcv::bam(form, family = fam, data = df, method = "fREML",
                 discrete = discrete, gamma = gamma)
  fit_time <- difftime(Sys.time(), t0, units = "mins")

  # ---- Mean linear-predictor matrix per arm (typical environment) ----
  t_grid <- seq(min(df$.time), max(df$.time), length.out = grid_n)
  labs <- vapply(m$smooth, function(s) s$label, character(1))
  cols_of <- function(lab) {
    sm <- m$smooth[[match(lab, labs)]]
    seq.int(sm$first.para, sm$last.para)
  }
  ge_lab <- labs[grepl(".ge", labs, fixed = TRUE)]
  unit_lab <- labs[grepl(".unit", labs, fixed = TRUE)]
  env_lev <- levels(df$.env)
  base_nd <- function(trt_i) {
    data.frame(.time = t_grid, .trt = factor(trt_i, levels = trt_levels),
               .env = factor(env_lev[1], levels = env_lev),
               .ge = factor(levels(df$.ge)[1], levels = levels(df$.ge)),
               .unit = factor(levels(df$.unit)[1], levels = levels(df$.unit)))
  }
  obs_ge <- unique(data.frame(trt = as.character(df$.trt), env = as.character(df$.env),
                              ge = as.character(df$.ge), stringsAsFactors = FALSE))
  arm_X <- lapply(trt_levels, function(trt_i) {
    nd <- base_nd(trt_i)
    X <- 0
    for (e in env_lev) {
      nd$.env <- factor(e, levels = env_lev)
      X <- X + mgcv::predict.bam(m, newdata = nd, type = "lpmatrix")
    }
    X <- X / length(env_lev)
    ge_cols <- unlist(lapply(ge_lab, cols_of))
    if (length(ge_cols)) {
      o <- obs_ge[obs_ge$trt == trt_i, , drop = FALSE]
      acc <- 0
      for (k in seq_len(nrow(o))) {
        nd$.env <- factor(o$env[k], levels = env_lev)
        nd$.ge <- factor(o$ge[k], levels = levels(df$.ge))
        acc <- acc + mgcv::predict.bam(m, newdata = nd, type = "lpmatrix")[, ge_cols, drop = FALSE]
      }
      X[, ge_cols] <- acc / max(1, nrow(o))
    }
    unit_cols <- unlist(lapply(unit_lab, cols_of))
    if (length(unit_cols)) X[, unit_cols] <- 0
    X
  })
  names(arm_X) <- trt_levels

  # ---- Posterior simulation ----
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  beta <- stats::coef(m)
  B <- mgcv::rmvn(n_sim, beta, m$Vp)
  linkinv <- m$family$linkinv
  scale_out <- if (response_scale == "percent") 100 else 1
  arm_mu <- lapply(arm_X, function(X) linkinv(X %*% t(B)) * scale_out)
  arm_hat <- lapply(arm_X, function(X) linkinv(as.numeric(X %*% beta)) * scale_out)

  band <- function(hat, draws) {
    se <- apply(draws, 1, stats::sd)
    z <- abs(sweep(draws, 1, hat)) / pmax(se, .Machine$double.eps)
    crit <- stats::quantile(apply(z, 2, max), 0.95, names = FALSE)
    tibble::tibble(lo = apply(draws, 1, stats::quantile, 0.025, names = FALSE),
                   hi = apply(draws, 1, stats::quantile, 0.975, names = FALSE),
                   sim_lo = hat - crit * se, sim_hi = hat + crit * se)
  }

  arms <- dplyr::bind_rows(lapply(trt_levels, function(trt_i) {
    dplyr::bind_cols(tibble::tibble(treatment = trt_i, time = t_grid, mu = arm_hat[[trt_i]]),
                     band(arm_hat[[trt_i]], arm_mu[[trt_i]]))
  }))

  ref_hat <- arm_hat[[reference]]
  ref_draws <- arm_mu[[reference]]
  contr <- function(ref, trt) {
    if (contrast == "difference") ref - trt else (ref - trt) / ref
  }
  others <- setdiff(trt_levels, reference)
  draws <- lapply(others, function(trt_i) contr(ref_draws, arm_mu[[trt_i]]))
  names(draws) <- others
  contrast_df <- dplyr::bind_rows(lapply(others, function(trt_i) {
    hat <- contr(ref_hat, arm_hat[[trt_i]])
    dplyr::bind_cols(tibble::tibble(treatment = trt_i, time = t_grid, DSP = hat),
                     band(hat, draws[[trt_i]]))
  }))
  names(contrast_df)[names(contrast_df) == "treatment"] <- treatment
  names(contrast_df)[names(contrast_df) == "time"] <- time
  names(arms)[names(arms) == "treatment"] <- treatment
  names(arms)[names(arms) == "time"] <- time

  structure(
    list(
      gam = m,
      arms = arms,
      contrast = contrast_df,
      draws = if (keep_draws) draws else NULL,
      grid = t_grid,
      observed_data = df,
      vars = list(time = time, response = response, treatment = treatment,
                  environment = environment, block = block, unit = unit),
      settings = list(reference = reference, response_scale = response_scale,
                      ge_term = ge_term, k_trt = k_trt, k_env = k_env, k_ge = k_ge,
                      unit_re = unit_re && !unit_is_ge, family = family,
                      contrast = contrast, n_sim = n_sim, seed = seed,
                      formula = paste(deparse(form, width.cutoff = 500L), collapse = " ")),
      fit_time = fit_time
    ),
    class = "functional_network"
  )
}

#' Print a functional network model
#'
#' @param x An object of class \code{"functional_network"}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.functional_network <- function(x, ...) {
  s <- summary(x$gam)
  cat("Functional network model of disease suppression (arm-based)\n")
  cat("Reference:", x$settings$reference, "| contrast:", x$settings$contrast, "\n")
  cat("Arms:", length(unique(x$arms[[x$vars$treatment]])),
      "| environments:", nlevels(x$observed_data$.env),
      "| units:", nlevels(x$observed_data$.unit), "| n:", s$n, "\n")
  cat("Formula:", x$settings$formula, "\n")
  cat(sprintf("Family: %s | R2 adj = %.3f | deviance explained = %.1f%% | fit time = %.1f min\n",
              x$settings$family, s$r.sq, 100 * s$dev.expl, as.numeric(x$fit_time)))
  invisible(x)
}
