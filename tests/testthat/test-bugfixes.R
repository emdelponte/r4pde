# Regression tests for the 0.2.3 bug fixes

test_that("AFSD summary mean SIF and CIF match the per-focus table", {
  set.seed(1)
  df <- data.frame(x = rep(1:20, each = 20), y = rep(1:20, 20),
                   i = sample(0:1, 400, replace = TRUE, prob = c(0.6, 0.4)))
  res <- AFSD(df)
  foci <- res[[2]]
  smry <- res[[1]]
  expect_equal(smry$value[smry$stats == "mean_SIF"], mean(foci$rows / foci$cols))
  expect_equal(smry$value[smry$stats == "mean_CIF"], mean(foci$CIF))
})

mock_fc_pca <- function() {
  t_grid <- seq(1, 10, length.out = 10)
  set.seed(2)
  curves <- expand.grid(time = t_grid, trt = paste0("T", 1:12))
  slope <- stats::setNames(rep(c(0.1, 0.5, 1), each = 4) + stats::rnorm(12, 0, 0.02),
                           paste0("T", 1:12))
  curves$mu <- slope[as.character(curves$trt)] * curves$time + 10 * (curves$time > 5) * (slope[as.character(curves$trt)] > 0.7)
  fc <- list(curves = curves, grid = t_grid, vars = list(time = "time", treatment = "trt"))
  class(fc) <- "functional_curves"
  fc
}

test_that("reconstruct_curves() undoes scaling", {
  fc <- mock_fc_pca()
  p <- functional_pca(fc, center = TRUE, scale = TRUE, var_explained = 1)
  rec <- reconstruct_curves(p, components = seq_len(ncol(p$pca$x)))
  expect_false(anyNA(rec$time))
  obs <- fc$curves
  m <- merge(rec, obs, by.x = c("curve_id", "time"), by.y = c("trt", "time"))
  expect_equal(nrow(m), nrow(obs))
  expect_equal(m$reconstructed, m$mu, tolerance = 1e-8)
})

test_that("functional_pca_clusters() chooses k by silhouette and elbow", {
  fc <- mock_fc_pca()
  p <- functional_pca(fc, center = TRUE, scale = FALSE, var_explained = 0.99)
  cl <- suppressMessages(functional_pca_clusters(p, choose_k = "silhouette", method = "hclust"))
  expect_equal(nlevels(cl$cluster), 3)
  expect_s3_class(attr(cl, "k_selection"), "tbl_df")
  cl2 <- suppressMessages(functional_pca_clusters(p, choose_k = "elbow", method = "hclust"))
  expect_true(nlevels(cl2$cluster) >= 2)
  expect_error(functional_pca_clusters(p), "Must specify")
})

test_that("persistence counts only time above the threshold", {
  # above at 1-3 and 6-8, gap at 4-5: 2 + 2, not the 7-unit span
  t <- 0:9
  y <- c(0, 2, 2, 2, 0, 0, 2, 2, 2, 0)
  expect_equal(r4pde:::time_above_internal(t, y, 0), 4)
  # no gap: same as the old first-to-last span
  expect_equal(r4pde:::time_above_internal(c(0, 5, 10, 15, 20), c(0, 10, 10, 10, 0), 5), 10)
})

test_that("functional_contrast() averages replicates explicitly", {
  d <- data.frame(trt = rep(c("Control", "F"), each = 6),
                  time = rep(c(0, 0, 10, 10, 20, 20), 2),
                  y = c(0, 2, 10, 12, 30, 34, 0, 0, 2, 4, 6, 8))
  expect_message(
    res <- functional_contrast(d, time = "time", response = "y", treatment = "trt",
                               reference = "Control"),
    "Averaging replicated"
  )
  expect_equal(res$DSP[res$time == 20], 32 - 7)
})

test_that("simulate_dsp_data() and the permutation test leave the RNG stream alone", {
  set.seed(42); a <- stats::runif(1)
  set.seed(42); invisible(simulate_dsp_data()); b <- stats::runif(1)
  expect_identical(a, b)
  expect_identical(simulate_dsp_data(), simulate_dsp_data())
})

test_that("fit_gradients() returns a numeric results table", {
  x <- c(0.8, 1.6, 2.4, 3.2, 4, 7.2, 12, 15.2, 21.6, 28.8)
  Y <- c(184.9, 113.3, 113.3, 64.1, 25, 8, 4.3, 2.5, 1, 0.8)
  mg <- fit_gradients(data.frame(x = x, Y = Y), C = 0.4)
  tab <- mg$results_table
  expect_s3_class(tab, "data.frame")
  expect_true(is.numeric(tab$b) && is.numeric(tab$R2) && is.numeric(tab$a_back))
  expect_equal(rownames(tab), c("Exponential", "Power", "Modified_Power"))
  expect_equal(tab["Power", "b"], round(unname(coef(lm(log(Y) ~ log(x)))[2]), 3))
})

test_that("theme_r4pde() returns a theme", {
  expect_s3_class(theme_r4pde(), "theme")
})
