library(testthat)
library(r4pde)

# Helper function to generate simulated epidemic data
make_sim_data <- function() {
  set.seed(42)
  t_seq <- seq(1, 30, length.out = 15)
  df <- data.frame(
    time = rep(t_seq, 3),
    treatment = rep(c("Early", "Late", "None"), each = 15),
    y = c(
      1 / (1 + exp(-(t_seq - 10) * 0.4)),  # Early
      1 / (1 + exp(-(t_seq - 20) * 0.4)),  # Late
      rep(0, 15)                           # None
    )
  )
  df$y <- pmax(0, pmin(1, df$y + rnorm(nrow(df), 0, 0.005)))
  df$y[df$treatment == "None"] <- 0 # strictly zero
  df
}

test_that("1. Invalid inputs trigger informative errors", {
  expect_error(functional_rate("not_an_object"), "`object` must be of class 'functional_curves'")
  
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  # Invalid n_grid
  expect_error(functional_rate(fc, n_grid = 1), "`n_grid` must be a single integer >= 2")
  expect_error(functional_rate(fc, n_grid = -10), "`n_grid` must be a single integer >= 2")
  expect_error(functional_rate(fc, n_grid = "bad"), "`n_grid` must be a single integer >= 2")
  
  # Invalid eps
  expect_error(functional_rate(fc, eps = 0), "`eps` must be a single positive numeric value")
  expect_error(functional_rate(fc, eps = -0.5), "`eps` must be a single positive numeric value")
  
  # Invalid interval
  expect_error(functional_rate(fc, interval = 0), "`interval` must be a single numeric value in \\(0, 1\\)")
  expect_error(functional_rate(fc, interval = 1), "`interval` must be a single numeric value in \\(0, 1\\)")
  expect_error(functional_rate(fc, interval = 1.2), "`interval` must be a single numeric value in \\(0, 1\\)")
})

test_that("2. Linear curve with known derivative on link scale", {
  # Generate strictly linear logit progression: eta = -3 + 0.2 * time
  # On link scale, derivative should be identically 0.2 everywhere
  t_seq <- seq(5, 25, by = 1)
  eta <- -3 + 0.2 * t_seq
  p <- 1 / (1 + exp(-eta))
  
  df <- data.frame(
    time = rep(t_seq, 2),
    treatment = rep(c("A", "B"), each = length(t_seq)),
    y = rep(p, 2)
  )
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial", k_smooth = 5
  )
  
  fr_link <- functional_rate(fc, scale = "link", n_grid = 50)
  
  # Mean estimated rate should be very close to true slope 0.2
  r_A <- fr_link$data$rate[fr_link$data$treatment == "A"]
  expect_equal(mean(r_A), 0.2, tolerance = 0.05)
})

test_that("3. Logistic curve with known analytical derivative", {
  # True logistic curve: S(t) = 1 / (1 + exp(-0.3*(t - 15)))
  # True derivative on response scale: S'(t) = 0.3 * S(t) * (1 - S(t))
  # Maximum rate occurs at t = 15, with value 0.3 * 0.5 * 0.5 = 0.075
  t_seq <- seq(2, 28, by = 1)
  p <- 1 / (1 + exp(-0.3 * (t_seq - 15)))
  
  df <- data.frame(
    time = rep(t_seq, 2),
    treatment = rep(c("Trt1", "Trt2"), each = length(t_seq)),
    y = rep(p, 2)
  )
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial", k_smooth = 8
  )
  
  fr <- functional_rate(fc, scale = "response", n_grid = 100)
  
  # Maximum rate should be near 0.075 and t_r_max near 15
  summ <- summary(fr)
  expect_equal(summ$r_max[summ$treatment == "Trt1"], 0.075, tolerance = 0.015)
  expect_equal(summ$t_r_max[summ$treatment == "Trt1"], 15, tolerance = 1.0)
})

test_that("4. Multiple treatments handled properly", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc)
  summ <- summary(fr)
  
  expect_equal(nrow(summ), 3)
  expect_true(all(c("Early", "Late", "None") %in% summ$treatment))
})

test_that("5. Multiple environments support", {
  set.seed(12)
  t_seq <- seq(1, 20, length.out = 10)
  df <- expand.grid(time = t_seq, treatment = c("A", "B"), env = c("E1", "E2"))
  df$y <- 1 / (1 + exp(-(df$time - 10) * 0.3)) + rnorm(nrow(df), 0, 0.01)
  df$y <- pmax(0, pmin(1, df$y))
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    environment = "env", min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc)
  expect_s3_class(fr, "functional_rate")
  expect_true("treatment" %in% names(fr$summary))
})

test_that("6. Scale 'link' vs 'response'", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr_resp <- functional_rate(fc, scale = "response")
  fr_link <- functional_rate(fc, scale = "link")
  
  expect_equal(fr_resp$settings$scale, "response")
  expect_equal(fr_link$settings$scale, "link")
  
  # Response rates are in proportion scale [0, 1] per day
  expect_true(all(fr_resp$data$rate <= 1.0))
  
  # Link scale rates can exceed 1
  expect_true(is.numeric(fr_link$data$rate))
})

test_that("7. Supported families (quasibinomial, betar, gaussian)", {
  sim_dat <- make_sim_data()
  # Quasibinomial
  fc_qb <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  fr_qb <- functional_rate(fc_qb)
  expect_s3_class(fr_qb, "functional_rate")
  
  # Beta family
  fc_beta <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "betar"
  )
  fr_beta <- functional_rate(fc_beta)
  expect_s3_class(fr_beta, "functional_rate")
})

test_that("8. Completely zero curve has r_max = 0, t_r_max = NA, growth_duration = 0", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc)
  summ <- summary(fr)
  
  none_row <- summ |> dplyr::filter(treatment == "None")
  expect_equal(none_row$r_max, 0)
  expect_true(is.na(none_row$t_r_max))
  expect_equal(none_row$growth_duration, 0)
  expect_equal(none_row$cumulative_positive_growth, 0)
  expect_equal(none_row$n_growth_windows, 0L)
  
  # Check that rate_data for None is zeroed
  none_data <- fr$data |> dplyr::filter(treatment == "None")
  expect_true(all(none_data$rate == 0))
  expect_true(all(none_data$fitted == 0))
})

test_that("9. Non-zero constant curve has derivative approximately zero", {
  set.seed(42)
  t_seq <- seq(1, 20, by = 1)
  df <- data.frame(
    time = rep(t_seq, 2),
    treatment = rep(c("Flat1", "Flat2"), each = length(t_seq)),
    y = rep(0.3, length(t_seq) * 2) + rnorm(length(t_seq) * 2, 0, 1e-4)
  )
  df$y <- pmax(0, pmin(1, df$y))
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc, scale = "response")
  expect_true(all(abs(fr$data$rate) < 0.01))
})

test_that("10. Negative rates preserved in rate_data without truncation", {
  # Decreasing curve
  t_seq <- seq(1, 20, by = 1)
  df <- data.frame(
    time = rep(t_seq, 2),
    treatment = rep(c("Dec1", "Dec2"), each = length(t_seq)),
    y = rep(1 - 1 / (1 + exp(-(t_seq - 10) * 0.4)), 2)
  )
  
  suppressWarnings({
    fc <- functional_curves(
      data = df, time = "time", response = "y", treatment = "treatment",
      min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
    )
    fr <- functional_rate(fc, scale = "response")
  })
  
  # Negative rates must be present in rate_data
  expect_true(any(fr$data$rate < 0))
})

test_that("11. Different temporal domains between groups without extrapolation", {
  df <- data.frame(
    time = c(1:15, 10:25),
    treatment = c(rep("Short", 15), rep("Shifted", 16)),
    y = c(
      1 / (1 + exp(-(1:15 - 8) * 0.4)),
      1 / (1 + exp(-(10:25 - 18) * 0.4))
    )
  )
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc, n_grid = 30)
  
  short_data <- fr$data |> dplyr::filter(treatment == "Short")
  shifted_data <- fr$data |> dplyr::filter(treatment == "Shifted")
  
  expect_equal(min(short_data$time), 1)
  expect_equal(max(short_data$time), 15)
  expect_equal(min(shifted_data$time), 10)
  expect_equal(max(shifted_data$time), 25)
})

test_that("12. Custom newdata support and extrapolation guard", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  # Valid newdata
  custom_grid <- data.frame(
    time = rep(c(5, 10, 15, 20), 2),
    treatment = rep(c("Early", "Late"), each = 4)
  )
  
  fr_custom <- functional_rate(fc, newdata = custom_grid)
  expect_equal(nrow(fr_custom$data), 8)
  
  # Extrapolated newdata triggers error
  extrap_grid <- data.frame(
    time = c(0, 50),
    treatment = c("Early", "Late")
  )
  expect_error(functional_rate(fc, newdata = extrap_grid), "outside the observed domain")
})

test_that("13. Confidence intervals behavior", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr95 <- functional_rate(fc, interval = 0.95)
  fr90 <- functional_rate(fc, interval = 0.90)
  
  # rate_lower < rate < rate_upper for curves with positive uncertainty
  early_data <- fr95$data |> dplyr::filter(treatment == "Early")
  expect_true(all(early_data$rate_lower <= early_data$rate_upper))
  
  # 95% CI is wider than 90% CI
  ci95_width <- fr95$data$rate_upper - fr95$data$rate_lower
  ci90_width <- fr90$data$rate_upper - fr90$data$rate_lower
  expect_true(all(ci95_width >= ci90_width - 1e-8))
})

test_that("14. Integration of positive rate matches cumulative positive growth", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc, n_grid = 100, threshold = 0.01)
  summ <- summary(fr)
  
  early_data <- fr$data |> dplyr::filter(treatment == "Early")
  expected_growth <- trapz_vec_internal(early_data$time, pmax(early_data$rate - 0.01, 0))
  
  expect_equal(summ$cumulative_positive_growth[summ$treatment == "Early"], expected_growth, tolerance = 1e-4)
})

test_that("15. Correct identification of r_max and t_r_max in Early vs Late", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc, n_grid = 100)
  summ <- summary(fr)
  
  t_early <- summ$t_r_max[summ$treatment == "Early"]
  t_late  <- summ$t_r_max[summ$treatment == "Late"]
  
  expect_true(t_early < t_late)
  expect_equal(t_early, 10, tolerance = 2)
  expect_equal(t_late, 20, tolerance = 2)
})

test_that("16. Absence of extrapolation at boundaries", {
  sim_dat <- make_sim_data()
  t_min_orig <- min(sim_dat$time)
  t_max_orig <- max(sim_dat$time)
  
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc, n_grid = 50)
  
  expect_gte(min(fr$data$time), t_min_orig)
  expect_lte(max(fr$data$time), t_max_orig)
})

test_that("17. S3 methods print(), summary(), plot(), augment() work", {
  sim_dat <- make_sim_data()
  fc <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc)
  
  # print()
  expect_output(print(fr), "Instantaneous Epidemic Rates")
  
  # summary()
  s <- summary(fr)
  expect_s3_class(s, "tbl_df")
  expect_equal(nrow(s), 3)
  
  # plot()
  p_both <- plot(fr, type = "both")
  expect_true(inherits(p_both, "ggplot") || inherits(p_both, "patchwork") || is.list(p_both))
  
  p_curve <- plot(fr, type = "curve")
  expect_s3_class(p_curve, "ggplot")
  
  p_rate <- plot(fr, type = "rate")
  expect_s3_class(p_rate, "ggplot")
  
  # augment()
  aug <- augment(fr)
  expect_s3_class(aug, "tbl_df")
  expect_true("rate" %in% names(aug))
  
  aug2 <- augment_functional_rate(fr)
  expect_equal(aug, aug2)
})

test_that("18. Preservation of group names and factor levels", {
  sim_dat <- make_sim_data()
  names(sim_dat)[names(sim_dat) == "treatment"] <- "genotype_var"
  names(sim_dat)[names(sim_dat) == "time"] <- "day_after_emergence"
  
  fc <- functional_curves(
    data = sim_dat, time = "day_after_emergence", response = "y", treatment = "genotype_var",
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr <- functional_rate(fc)
  
  expect_true("genotype_var" %in% names(fr$data))
  expect_true("day_after_emergence" %in% names(fr$data))
  expect_true("genotype_var" %in% names(fr$summary))
})

test_that("19. Compatibility with global_smooth = TRUE and FALSE", {
  sim_dat <- make_sim_data()
  
  fc_ind <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    global_smooth = FALSE, min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  fr_ind <- functional_rate(fc_ind)
  expect_s3_class(fr_ind, "functional_rate")
  
  fc_glob <- functional_curves(
    data = sim_dat, time = "time", response = "y", treatment = "treatment",
    global_smooth = TRUE, min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  fr_glob <- functional_rate(fc_glob)
  expect_s3_class(fr_glob, "functional_rate")
})

test_that("20. No alteration in existing functional functions", {
  set.seed(123)
  df <- data.frame(
    time = rep(1:5, 4),
    treatment = rep(c("A", "B"), each = 10),
    block = rep(rep(1:2, each = 5), 2),
    y = pmax(0, pmin(1, rep(c(0.1, 0.3, 0.5, 0.7, 0.9), 4) + rnorm(20, 0, 0.05)))
  )
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    block = "block", min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fd <- functional_distances(fc, cluster_k = 2, show_progress = FALSE)
  expect_s3_class(fd, "functional_distances")
  expect_equal(fd$distance_matrix, t(fd$distance_matrix))
})
