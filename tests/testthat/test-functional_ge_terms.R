# Environment factor-smooth and genotype-by-environment terms in functional_curves()

sim_ge_data <- function(seed = 1, n_env = 6, reps = 1) {
  set.seed(seed)
  geno <- paste0("G", 1:4)
  env <- paste0("L", rep(1:3, length.out = n_env), "_", rep(2021:2022, each = n_env / 2))
  d <- expand.grid(geno = geno, env = env, block = seq_len(reps),
                   time = seq(0, 60, by = 6), stringsAsFactors = FALSE)
  g_eff <- stats::setNames(c(-1.5, -0.5, 0.5, 1.5), geno)
  e_eff <- stats::setNames(stats::rnorm(n_env, 0, 0.6), env)
  key <- paste(d$geno, d$env)
  ge <- stats::setNames(stats::rnorm(length(unique(key)), 0, 0.5), unique(key))[key]
  eta <- -3.5 + g_eff[d$geno] + e_eff[d$env] + ge + 0.12 * (d$time - 30)
  mu <- stats::plogis(eta)
  d$y <- stats::rbeta(nrow(d), mu * 40, (1 - mu) * 40)
  d
}

fit_fc <- function(d, ...) {
  suppressWarnings(functional_curves(
    data = d, time = "time", response = "y", treatment = "geno",
    environment = "env", show_progress = FALSE, family_try = "quasibinomial", ...
  ))
}

test_that("default model is unchanged and has no G x E term", {
  d <- sim_ge_data()
  fc <- fit_fc(d)
  labs <- vapply(fc$gam$smooth, function(s) s$label, character(1))
  expect_true("s(env)" %in% labs)
  expect_false(any(grepl(".ge", labs, fixed = TRUE)))
  expect_false(fc$settings$env_smooth)
  expect_identical(fc$settings$ge_term, "none")
  expect_warning(functional_instability(fc), "no genotype-by-environment term")
})

test_that("env_smooth fits an environment factor-smooth", {
  d <- sim_ge_data()
  fc <- fit_fc(d, env_smooth = TRUE, k_env = 4)
  labs <- vapply(fc$gam$smooth, function(s) s$label, character(1))
  expect_true("s(time,env)" %in% labs)
  expect_false("s(env)" %in% labs)
  expect_match(fc$settings$formula, "bs = \"fs\"", fixed = TRUE)
  expect_false(anyNA(fc$curves$mu))
  expect_true(all(fc$curves$mu > 0 & fc$curves$mu < 1))
})

test_that("treatment mean curves stay close under the random intercept and the factor-smooth", {
  d <- sim_ge_data(seed = 3)
  fc_re <- fit_fc(d)
  fc_fs <- fit_fc(d, env_smooth = TRUE)
  # Factor-smooth deviations are averaged over environments rather than dropped,
  # so the mean curves must not collapse towards zero.
  expect_lt(max(abs(fc_re$curves$mu - fc_fs$curves$mu)), 0.1)
})

test_that("ge_term adds intercept or smooth G x E terms", {
  d <- sim_ge_data()
  fc_i <- fit_fc(d, ge_term = "intercept")
  fc_s <- fit_fc(d, env_smooth = TRUE, ge_term = "smooth", k_ge = 4)
  labs_i <- vapply(fc_i$gam$smooth, function(s) s$label, character(1))
  labs_s <- vapply(fc_s$gam$smooth, function(s) s$label, character(1))
  expect_true("s(.ge)" %in% labs_i)
  expect_true("s(time,.ge)" %in% labs_s)
  expect_identical(fc_s$vars$ge, ".ge")
  expect_false(anyNA(fc_s$curves$mu))
})

test_that("env_smooth and ge_term without environment warn and are ignored", {
  d <- sim_ge_data()
  expect_warning(
    fc <- functional_curves(d, time = "time", response = "y", treatment = "geno",
                            env_smooth = TRUE, ge_term = "smooth",
                            show_progress = FALSE, family_try = "quasibinomial"),
    "require `environment`"
  )
  expect_identical(fc$settings$ge_term, "none")
})

test_that("unit random intercept is dropped when units are the G x E curves", {
  d <- sim_ge_data()
  d$curve <- interaction(d$geno, d$env)
  expect_message(
    fc <- fit_fc(d, unit = "curve", ge_term = "smooth"),
    "unit random intercept is dropped"
  )
  labs <- vapply(fc$gam$smooth, function(s) s$label, character(1))
  expect_false("s(curve)" %in% labs)
  expect_true(fc$settings$unit_is_ge)
})

test_that("functional_instability works with G x E terms and with blocks", {
  d <- sim_ge_data(reps = 2)
  fc <- fit_fc(d, env_smooth = TRUE, ge_term = "smooth", block = "block")
  expect_true(fc$settings$unit_term)
  expect_silent(fi <- functional_instability(fc, env_sep = "_"))
  expect_equal(nrow(fi), 4)
  expect_true(all(is.finite(fi$nFI)))
  expect_true(all(c("nFI_space", "nFI_time") %in% names(fi)))
})

test_that("a curve-level unit intercept counts as a G x E term", {
  d <- sim_ge_data()
  d$curve <- interaction(d$geno, d$env)
  fc <- fit_fc(d, unit = "curve")
  expect_true(fc$settings$unit_is_ge)
  expect_silent(fi <- functional_instability(fc))
  expect_true(all(is.finite(fi$nFI)))
})

test_that("downstream functions use the same mean curves", {
  d <- sim_ge_data(reps = 2)
  d$cls <- ifelse(d$geno %in% c("G1", "G2"), "A", "B")
  fc <- fit_fc(d, env_smooth = TRUE, ge_term = "smooth", block = "block")

  fr <- suppressWarnings(functional_rate(fc, n_grid = length(fc$grid)))
  for (g in levels(fc$curves$geno)) {
    expect_equal(fr$data$fitted[fr$data$geno == g], fc$curves$mu[fc$curves$geno == g],
                 tolerance = 1e-6)
  }

  fd <- functional_distances(fc, cluster_k = 2, test_factor = "cls", n_perm = 19,
                             show_progress = FALSE, bootstrap = TRUE, boot_B = 5)
  expect_s3_class(fd, "functional_distances")
  expect_equal(nrow(fd$curve_distance), length(unique(fc$observed_data$.unit)))

  fres <- functional_resistance(fc, reference = "G4", group_method = "bootstrap",
                                n_boot = 5, n_groups = 2)
  expect_true("FRI" %in% names(fres$table))

  dg <- diagnose_curves(fc)
  expect_s3_class(dg, "r4pde_curve_diagnostics")
})

test_that("functional_resistance() accepts functional_instability() output", {
  d <- sim_ge_data()
  fc <- fit_fc(d, ge_term = "intercept")
  fi <- functional_instability(fc)
  fi_list <- functional_instability(fc, return_curves = TRUE)
  r1 <- functional_resistance(fc, instability = fi, reference = "G4")
  r2 <- functional_resistance(fc, instability = fi_list, reference = "G4")
  expect_true(all(is.finite(r1$table$nFI)))
  expect_equal(r1$table$SAFRI, r2$table$SAFRI)
})

test_that("permutation test keeps genotype curves together", {
  d <- sim_ge_data(n_env = 6)
  d$curve <- interaction(d$geno, d$env)
  d$cls <- ifelse(d$geno %in% c("G1", "G2"), "A", "B")
  fc <- fit_fc(d, unit = "curve")
  set.seed(5); before <- stats::runif(1)
  set.seed(5)
  fd <- functional_distances(fc, cluster_k = 2, test_factor = "cls", perm_unit = "geno",
                             n_perm = 49, test_mode = "global", show_progress = FALSE)
  expect_identical(stats::runif(1), before)
  # with 4 genotypes split 2/2 there are only 3 distinct relabellings, so the
  # unit-level null distribution takes at most 3 distinct F values
  expect_lte(length(unique(round(fd$test$perm_F, 8))), 3)
  fd2 <- functional_distances(fc, cluster_k = 2, test_factor = "cls", perm_unit = "geno",
                              n_perm = 49, test_mode = "global", show_progress = FALSE)
  expect_identical(fd$test$global, fd2$test$global)
})
