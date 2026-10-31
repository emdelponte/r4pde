sad_x <- function(...) {
  sad_data(SADRustSoybean, wide = c(UN = "unaided", SAD = "aided"), group = "sad", ...)
}

test_that("sad_data builds conditions and detects a paired design", {
  x <- sad_x()
  expect_s3_class(x, "sad_data")
  expect_equal(nrow(x), 4000)
  expect_equal(levels(x$condition), c("UN_Old", "UN_New", "SAD_Old", "SAD_New"))
  expect_equal(attr(x, "design"), "paired")
  expect_equal(x$error, x$estimate - x$actual)
})

test_that("sad_data accepts long data and warns about problems", {
  long <- data.frame(rater = rep(1:3, each = 4), leaf = rep(1:2, 6),
                     actual = rep(c(10, 30), 6),
                     estimate = c(12, 28, 15, 33, 9, 120, 11, 31, 10, 29, 14, 35),
                     method = rep(rep(c("a", "b"), each = 2), 3))
  expect_warning(x <- sad_data(long), "outside 0-100")
  expect_equal(levels(x$method), c("a", "b"))
  bad <- long; bad$actual[1] <- 11
  expect_warning(sad_data(bad), "more than one actual")
  expect_error(sad_data(long, rater = "nope"), "not found")
})

test_that("Lin's statistics are correct", {
  set.seed(1)
  a <- runif(30, 0, 80); e <- a * 1.1 + 2 + rnorm(30)
  s <- r4pde:::.lin_stats(e, a)
  n <- length(a)
  s1 <- sqrt(var(e) * (n - 1) / n); s2 <- sqrt(var(a) * (n - 1) / n)
  rc <- 2 * cov(e, a) * (n - 1) / n / (s1^2 + s2^2 + (mean(e) - mean(a))^2)
  expect_equal(s$rho_c, rc)
  expect_equal(s$rho_c, s$r * s$Cb)
  expect_gt(s$u, 0)
  perfect <- r4pde:::.lin_stats(a, a)
  expect_equal(perfect$rho_c, 1)
})

test_that("sad_ccc and sad_compare reproduce Franceschi et al. (2020)", {
  skip_if_not_installed("emmeans"); skip_if_not_installed("lmerTest")
  skip_if_not_installed("multcompView")
  x <- sad_x()
  cc <- sad_ccc(x)
  expect_equal(nrow(cc), 80)
  cmp <- sad_compare(x, cc)
  expect_s3_class(cmp, "sad_compare")
  est <- cmp$estimates
  rho <- est[est$stat == "rho_c", ]
  expect_equal(round(rho$mean[rho$condition == "SAD_New"], 2), 0.96)
  expect_equal(rho$letter[rho$condition == "SAD_New"], "a")
  expect_true(all(c("condition", "N", "rho_c") %in% names(cmp$table)))
  expect_output(print(cmp), "same letter")
})

test_that("sad_compare works with a single SAD", {
  skip_if_not_installed("emmeans"); skip_if_not_installed("lmerTest")
  skip_if_not_installed("multcompView")
  x1 <- sad_data(subset(SADRustSoybean, sad == "New"), wide = c(UN = "unaided", SAD = "aided"))
  cmp <- sad_compare(x1, stats = "rho_c")
  expect_equal(nrow(cmp$table), 2)
  expect_equal(nrow(cmp$pairwise), 1)
})

test_that("sad_reliability matches known ICC and OCC", {
  x <- sad_x()
  rel <- sad_reliability(x)
  expect_equal(nrow(rel), 4)
  expect_equal(round(rel$ICC[rel$condition == "SAD_New"], 2), 0.94)
  expect_true(all(rel$ICC_lower < rel$ICC & rel$ICC < rel$ICC_upper))
  skip_if_not_installed("epiR")
  M <- as.matrix(tidyr::pivot_wider(as.data.frame(x[x$condition == "SAD_New", c("image", "rater", "estimate")]),
                                    names_from = "rater", values_from = "estimate")[, -1])
  expect_equal(rel$OCC[rel$condition == "SAD_New"], epiR::epi.occc(M, pairs = FALSE)$occc,
               tolerance = 1e-6)
})

test_that("sad_plots returns ggplot objects", {
  x <- sad_x()
  p <- sad_plots(x, baseline = "UN")
  expect_true(all(vapply(p, ggplot2::is.ggplot, logical(1))))
  expect_true("gain_vs_baseline" %in% names(p))
  expect_error(sad_plots(x, baseline = "nope"), "baseline")
})

test_that("one-stage model, curves, raters and images", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("MASS")
  x1 <- sad_data(subset(SADRustSoybean, sad == "New"), wide = c(UN = "unaided", SAD = "aided"))
  fit <- suppressMessages(sad_fit(x1, dispersion = "lin"))
  expect_s3_class(fit, "sad_fit")
  expect_output(print(fit), "One-stage model")
  cv <- sad_curve(fit, nsim = 500)
  expect_s3_class(cv, "sad_curve")
  expect_equal(unique(cv$differences$contrast), "SAD - UN")
  expect_true(all(cv$curves$expected_abs_error > 0))
  # the new SAD lowers the expected absolute error at 50% severity
  e50 <- sapply(split(cv$curves, cv$curves$condition), function(d)
    stats::approx(d$actual, d$expected_abs_error, 50)$y)
  expect_lt(e50["SAD"], e50["UN"])
  reg <- sad_regions(cv)
  expect_named(reg, c("contrast", "from", "to", "direction"))
  expect_s3_class(plot(cv, "difference"), "ggplot")
  expect_s3_class(plot(cv, "expected_abs_error"), "ggplot")
  rt <- sad_raters(fit)
  expect_equal(nrow(rt), 40)
  expect_true(all(abs(rt$slope - 1) < 1))
  expect_s3_class(plot_sad_raters(rt), "ggplot")
  im <- sad_images(fit)
  expect_equal(nrow(im), 50)
  seed_before <- .Random.seed
  cs <- sad_curve(fit, type = "smooth", nsim = 500)
  expect_identical(.Random.seed, seed_before)
  expect_s3_class(plot(cs), "ggplot")
  expect_error(plot(cs, "expected_abs_error"), "linear")
})

test_that("reference by method pairs conditions within groups", {
  x <- sad_x()
  d <- tibble::as_tibble(unclass(x))
  p <- r4pde:::.sad_pairs(d, "UN")
  expect_equal(p$a, c("SAD_Old", "SAD_New"))
  expect_equal(p$b, c("UN_Old", "UN_New"))
  expect_error(r4pde:::.sad_pairs(d, "nope"), "reference")
})

test_that("sad_simulate runs", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  s <- sad_simulate(n_rater = 6, n_image = 10, nsim = 2)
  expect_equal(nrow(s), 1)
  expect_true(all(c("power_two_stage", "power_one_stage") %in% names(s)))
})
