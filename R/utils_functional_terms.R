# Internal helpers shared by functions that predict from a functional_curves() GAM.
# They keep the handling of environment, genotype-by-environment (G x E) and
# experimental-unit terms in one place, so every downstream function excludes
# and fills the same terms.

#' Labels of GAM smooths involving given variables
#' @noRd
.fc_smooth_labels <- function(m_gam, vars) {
  vars <- vars[!is.null(vars) & nzchar(vars)]
  if (is.null(m_gam$smooth) || length(vars) == 0) return(character(0))
  labs <- vapply(m_gam$smooth, function(s) {
    if (any(s$term %in% vars)) s$label else NA_character_
  }, character(1))
  unique(labs[!is.na(labs)])
}

#' Is the experimental unit the same partition as treatment x environment?
#' @noRd
.fc_unit_is_ge <- function(df, unit, ge) {
  if (is.null(unit) || is.null(ge) || !all(c(unit, ge) %in% names(df))) return(FALSE)
  u <- as.character(df[[unit]])
  g <- as.character(df[[ge]])
  n_u <- length(unique(u))
  n_u == length(unique(g)) && n_u == length(unique(paste(u, g, sep = "\r")))
}

#' Smooth labels for environment, G x E and unit terms of a functional_curves fit
#' @noRd
.fc_random_terms <- function(object) {
  m_gam <- object$gam
  v <- object$vars
  list(
    env  = .fc_smooth_labels(m_gam, v$environment),
    ge   = .fc_smooth_labels(m_gam, v$ge),
    unit = .fc_smooth_labels(m_gam, v$unit)
  )
}

#' Does a functional_curves fit contain a genotype-by-environment term?
#'
#' True when the model has an explicit G x E term (`ge_term`), or a unit
#' random effect whose units are exactly the treatment x environment curves.
#' @noRd
.fc_has_ge <- function(object) {
  rt <- .fc_random_terms(object)
  if (length(rt$ge) > 0) return(TRUE)
  length(rt$unit) > 0 && .fc_unit_matches_ge(object)
}

#' Are the experimental units exactly the treatment x environment curves?
#' Works for objects fitted before `unit_is_ge` was stored.
#' @noRd
.fc_unit_matches_ge <- function(object) {
  if (!is.null(object$settings$unit_is_ge)) return(isTRUE(object$settings$unit_is_ge))
  v <- object$vars
  df <- object$observed_data
  if (is.null(v$environment) || is.null(v$unit)) return(FALSE)
  df$.ge_tmp <- .fc_ge_id(df[[v$treatment]], df[[v$environment]])
  .fc_unit_is_ge(df, v$unit, ".ge_tmp")
}

#' Add factor columns a GAM needs for prediction but that newdata lacks
#'
#' Columns for environment, block, unit and the G x E identifier are filled
#' with the first observed level, so `predict.gam()` finds every variable.
#' Use together with `exclude =` for terms that should not contribute.
#' @noRd
.fc_fill_newdata <- function(newdata, object) {
  df <- object$observed_data
  v <- object$vars
  cols <- unique(c(v$environment, v$block, v$unit, v$ge))
  for (cc in cols) {
    if (is.null(cc) || !(cc %in% names(df)) || cc %in% names(newdata)) next
    lev <- levels(factor(df[[cc]]))
    newdata[[cc]] <- factor(lev[1], levels = lev)
  }
  newdata
}

#' Build the G x E identifier used by functional_curves()
#' @noRd
.fc_ge_id <- function(trt, env) {
  interaction(trt, env, drop = TRUE, sep = "|")
}

#' Linear predictor matrix for treatment mean curves
#'
#' Returns the `lpmatrix` rows that define the treatment mean curve at the
#' rows of `newdata` (treatment, time and any covariates), so that
#' `X %*% coef(gam)` is the mean curve on the link scale. Unit terms are
#' always excluded. Environment random intercepts are excluded (their
#' predicted values are centred). Environment factor-smooths are not centred
#' by construction, so they are averaged over all environments instead of
#' being dropped. Genotype-by-environment terms are averaged over the
#' environments in which each treatment was observed.
#' @noRd
.fc_mean_lpmatrix <- function(object, newdata) {
  m <- object$gam
  v <- object$vars
  df <- object$observed_data
  rt <- .fc_random_terms(object)

  center_env <- isTRUE(object$settings$env_smooth) && length(rt$env) > 0
  avg_ge <- length(rt$ge) > 0
  excl <- unique(c(rt$ge, rt$unit, if (!center_env) rt$env))

  nd <- .fc_fill_newdata(tibble::as_tibble(newdata), object)

  if (center_env) {
    env_lev <- levels(factor(df[[v$environment]]))
    X <- NULL
    for (e in env_lev) {
      nd_e <- nd
      nd_e[[v$environment]] <- factor(e, levels = env_lev)
      X_e <- mgcv::predict.gam(m, newdata = nd_e, type = "lpmatrix", exclude = excl)
      X <- if (is.null(X)) X_e else X + X_e
    }
    X <- X / length(env_lev)
  } else {
    X <- mgcv::predict.gam(m, newdata = nd, type = "lpmatrix", exclude = excl)
  }

  if (avg_ge) {
    ge_cols <- unlist(lapply(m$smooth, function(sm) {
      if (sm$label %in% rt$ge) seq.int(sm$first.para, sm$last.para) else NULL
    }))
    env_lev <- levels(factor(df[[v$environment]]))
    ge_lev <- levels(factor(df[[v$ge]]))
    obs <- unique(data.frame(
      trt = as.character(df[[v$treatment]]),
      env = as.character(df[[v$environment]]),
      ge  = as.character(df[[v$ge]]),
      stringsAsFactors = FALSE
    ))
    trt_new <- as.character(nd[[v$treatment]])
    for (g in unique(trt_new)) {
      rows <- which(trt_new == g)
      o <- obs[obs$trt == g, , drop = FALSE]
      if (nrow(o) == 0) next
      acc <- 0
      for (k in seq_len(nrow(o))) {
        nd_k <- nd[rows, , drop = FALSE]
        nd_k[[v$environment]] <- factor(o$env[k], levels = env_lev)
        nd_k[[v$ge]] <- factor(o$ge[k], levels = ge_lev)
        X_k <- mgcv::predict.gam(m, newdata = nd_k, type = "lpmatrix")
        acc <- acc + X_k[, ge_cols, drop = FALSE]
      }
      X[rows, ge_cols] <- acc / nrow(o)
    }
  }

  X
}

#' Evaluate code with a temporary RNG seed, restoring the caller's stream
#'
#' With `seed = NULL` the code uses the current random stream. Otherwise the
#' seed is set for the duration of `code` and the previous `.Random.seed` is
#' put back afterwards, so the user's own random sequence is not reset.
#' @noRd
.with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  set.seed(seed)
  code
}
