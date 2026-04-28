test_that("suggest_k works with explicit n_time and n_env", {
  out <- suggest_k(n_time = 5, n_env = 8, rule = "minimum", smoothness = "conservative")

  expect_s3_class(out, "suggest_k")
  expect_equal(out$effective_n_time, 5L)
  expect_equal(out$effective_n_env, 8L)
  expect_true(out$k_smooth >= 4)
  expect_true(out$k_trt < out$k_smooth)
  expect_equal(out$gamma, 1.8)
  expect_null(out$time_summary)
})

test_that("suggest_k infers time points per treatment-environment correctly", {
  df <- data.frame(
    time     = c(1:6, 1:5, 1:6, 1:5),   # trt A env E1: 6, trt A env E2: 5, trt B env E1: 6, trt B env E2: 5
    trt      = rep(c("A", "A", "B", "B"), times = c(6, 5, 6, 5)),
    env      = rep(c("E1", "E2", "E1", "E2"), times = c(6, 5, 6, 5)),
    severity = runif(22, 0, 0.5)
  )

  out_min <- suggest_k(data = df, time = time, treatment = trt,
                       environment = env, rule = "minimum")
  out_med <- suggest_k(data = df, time = time, treatment = trt,
                       environment = env, rule = "median")

  # minimum rule: smallest unique time count across trt x env is 5
  expect_equal(out_min$effective_n_time, 5L)
  # median rule: median of c(6,5,6,5) = 5.5 -> 5 as integer
  expect_equal(out_med$effective_n_time, 5L)

  # time_summary should have min = 5, max = 6
  expect_equal(out_min$time_summary[["minimum"]], 5)
  expect_equal(out_min$time_summary[["maximum"]], 6)
})

test_that("suggest_k infers environments per treatment correctly", {
  df <- data.frame(
    time     = rep(1:5, times = 9),
    trt      = rep(c("A", "B", "C"), each = 15),
    env      = rep(rep(c("E1", "E2", "E3"), each = 5), times = 3),
    severity = runif(45, 0, 0.5)
  )

  out <- suggest_k(data = df, time = time, treatment = trt,
                   environment = env, rule = "minimum")

  # each treatment appears in 3 environments
  expect_equal(out$environment_summary[["minimum"]], 3)
  expect_equal(out$environment_summary[["maximum"]], 3)
  expect_equal(out$effective_n_env, 3L)
})

test_that("rule = 'minimum' uses the smallest replication", {
  # Unbalanced: trt A has env E1=6 pts, E2=3 pts
  df <- data.frame(
    time = c(1:6, 1:3),
    trt  = c(rep("A", 6), rep("A", 3)),
    env  = c(rep("E1", 6), rep("E2", 3)),
    y    = runif(9)
  )
  out_min <- suggest_k(data = df, time = time, treatment = trt, environment = env, rule = "minimum")
  out_med <- suggest_k(data = df, time = time, treatment = trt, environment = env, rule = "median")

  expect_equal(out_min$effective_n_time, 3L)   # minimum is 3
  expect_equal(out_med$effective_n_time, 4L)   # median of c(6, 3) = 4.5 -> 4
})

test_that("rule = 'median' uses median replication", {
  df <- data.frame(
    time = c(1:3, 1:5, 1:7),
    trt  = rep(c("A", "B", "C"), times = c(3, 5, 7)),
    y    = runif(15)
  )
  out <- suggest_k(data = df, time = time, treatment = trt, rule = "median")
  # median of c(3, 5, 7) = 5
  expect_equal(out$effective_n_time, 5L)
})

test_that("sparse data returns low k values", {
  out <- suggest_k(n_time = 4, n_env = 3, smoothness = "conservative")
  expect_lte(out$k_smooth, 4)
  expect_lte(out$k_trt, out$k_smooth)
  expect_lte(out$k_env, 3)
  expect_equal(out$gamma, 1.8)  # most conservative gamma
})

test_that("smoothness levels produce ordered gamma values", {
  o_con <- suggest_k(n_time = 8, n_env = 10, smoothness = "conservative")
  o_mod <- suggest_k(n_time = 8, n_env = 10, smoothness = "moderate")
  o_fle <- suggest_k(n_time = 8, n_env = 10, smoothness = "flexible")

  expect_gt(o_con$gamma, o_mod$gamma)
  expect_gt(o_mod$gamma, o_fle$gamma)
})

test_that("suggest_k works without environment variable", {
  df <- data.frame(
    time = rep(1:6, times = 2),
    trt  = rep(c("A", "B"), each = 6),
    y    = runif(12)
  )
  out <- suggest_k(data = df, time = time, treatment = trt)
  expect_null(out$k_env)
  expect_null(out$environment_summary)
  expect_null(out$effective_n_env)
})

test_that("suggest_k errors informatively on bad inputs", {
  expect_error(suggest_k(), "n_time")
  expect_error(suggest_k(data = data.frame(a = 1), time = missing_col, treatment = a),
               "not found in")
})

test_that("print.suggest_k runs without error", {
  out <- suggest_k(n_time = 6, n_env = 4)
  expect_output(print(out), "functional_curves")
})
