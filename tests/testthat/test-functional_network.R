# Arm-based functional network model of disease suppression

sim_network_data <- function(seed = 1, n_env = 6, reps = 3) {
  set.seed(seed)
  trt <- c("Control", "A", "B", "C")
  env <- paste0("E", seq_len(n_env))
  d <- expand.grid(treatment = trt, env = env, block = seq_len(reps),
                   time = seq(0, 60, by = 10), stringsAsFactors = FALSE)
  # unbalanced network: B only in the first half of the environments,
  # C only in the second half
  d <- d[!(d$treatment == "B" & d$env %in% env[-(1:(n_env / 2))]), ]
  d <- d[!(d$treatment == "C" & d$env %in% env[1:(n_env / 2)]), ]
  eff <- c(Control = 0, A = -1.5, B = -0.8, C = -2)
  e_eff <- stats::setNames(stats::rnorm(n_env, 0, 0.5), env)
  eta <- -4 + 0.1 * d$time + eff[d$treatment] * (d$time / 60) + e_eff[d$env]
  mu <- stats::plogis(eta)
  d$sev <- 100 * stats::rbeta(nrow(d), mu * 50, (1 - mu) * 50)
  d
}

test_that("functional_network fits an arm-based model and returns suppression curves", {
  d <- sim_network_data()
  fn <- functional_network(
    data = d, time = "time", response = "sev", treatment = "treatment",
    environment = "env", reference = "Control", block = "block",
    k_trt = 5, k_env = 3, k_ge = 3, grid_n = 30, n_sim = 100
  )
  expect_s3_class(fn, "functional_network")
  expect_true(all(c("gam", "arms", "contrast", "draws", "settings") %in% names(fn)))

  # one arm curve per treatment (reference included), one suppression curve per non-reference
  expect_setequal(unique(fn$arms$treatment), c("Control", "A", "B", "C"))
  expect_setequal(unique(fn$contrast$treatment), c("A", "B", "C"))
  expect_equal(nrow(fn$contrast), 3 * 30)
  expect_true(all(c("DSP", "lo", "hi", "sim_lo", "sim_hi") %in% names(fn$contrast)))

  # bands contain the estimate and the simultaneous band is at least as wide
  expect_true(all(fn$contrast$lo <= fn$contrast$DSP + 1e-8 &
                  fn$contrast$hi >= fn$contrast$DSP - 1e-8))
  expect_true(mean(fn$contrast$sim_hi - fn$contrast$sim_lo) >=
                mean(fn$contrast$hi - fn$contrast$lo))

  # suppression = reference arm - treatment arm, on the response scale (percent)
  ref <- fn$arms$mu[fn$arms$treatment == "Control"]
  a <- fn$arms$mu[fn$arms$treatment == "A"]
  expect_equal(fn$contrast$DSP[fn$contrast$treatment == "A"], ref - a, tolerance = 1e-8)

  # treatment C (strongest effect, never together with B) suppresses most at the end
  end <- fn$contrast[fn$contrast$time == max(fn$contrast$time), ]
  expect_equal(end$treatment[which.max(end$DSP)], "C")

  # model terms
  labs <- vapply(fn$gam$smooth, function(s) s$label, character(1))
  expect_true(any(grepl(".env", labs, fixed = TRUE)))
  expect_true(any(grepl(".ge", labs, fixed = TRUE)))
  expect_true(any(grepl(".unit", labs, fixed = TRUE)))
  expect_output(print(fn), "Functional network model")
})

test_that("functional_network supports ge_term options and relative contrast", {
  d <- sim_network_data(seed = 2)
  fn_i <- functional_network(
    data = d, time = "time", response = "sev", treatment = "treatment",
    environment = "env", reference = "Control", block = "block",
    ge_term = "intercept", contrast = "relative",
    k_trt = 5, k_env = 3, grid_n = 20, n_sim = 50, keep_draws = FALSE
  )
  labs <- vapply(fn_i$gam$smooth, function(s) s$label, character(1))
  expect_true("s(.ge)" %in% labs)
  expect_null(fn_i$draws)
  expect_true(all(fn_i$contrast$DSP <= 1))
  expect_error(
    functional_network(d, "time", "sev", "treatment", "env", reference = "Nope"),
    "not found"
  )
})
