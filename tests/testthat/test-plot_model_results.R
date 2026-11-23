toy_run <- function(n = 21, H0 = 100, I0 = 1) {
  time <- seq(0, 100, length.out = n)
  frac <- 1 / (1 + exp(-0.1 * (time - 50)))
  data.frame(time = time,
             H = H0 * (1 - frac),
             L = H0 * frac * 0.2,
             I = H0 * frac * 0.3 + I0,
             R = H0 * frac * 0.5,
             P = 1)
}

test_that("HLIR-type plots have the expected compartments", {
  res <- toy_run()
  p <- plot_model_results(res, "HLIR", print = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$Compartment), c("H", "L", "I", "R", "D"))
  d <- p$data[p$data$Compartment == "D", ]
  expect_equal(d$Value, res$L + res$I + res$R)
  expect_equal(p$labels$y, "Number of hosts")
  expect_equal(p$labels$title, "HLIR model")

  g <- plot_model_results(res, "HLIRP (with host growth)", print = FALSE)
  expect_equal(levels(g$data$Compartment), c("H", "L", "I", "R", "D", "T_tissue"))
  tt <- g$data[g$data$Compartment == "T_tissue", ]
  expect_equal(tt$Value, res$H + res$L + res$I + res$R)
  expect_equal(g$labels$y, "Amount of host tissue")
})

test_that("printing returns the plot invisibly", {
  res <- toy_run()
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(p <- plot_model_results(res, "HLIRP"))
  expect_s3_class(p, "ggplot")
})

test_that("MULTI returns progress and inoculum plots", {
  multi <- list(n.yr = 3,
                results = lapply(1:3, function(k) toy_run(I0 = k)),
                inoc = c(1, 2, 4),
                n.max = 100)
  out <- plot_model_results(multi, "MULTI", print = FALSE)
  expect_named(out, c("progress", "inoculum"))
  expect_s3_class(out$progress, "ggplot")
  expect_equal(range(out$progress$data$time), c(0, (2 * 365 + 100) / 365))
  expect_equal(out$inoculum$data$log_Inoc, log(c(1, 2, 4)))
})

test_that("bad input gives informative errors", {
  expect_error(plot_model_results(toy_run()[, c("time", "H")], "HLIR"), "missing: L, I, R")
  expect_error(plot_model_results(toy_run(), "SIR"), "should be one of")
  expect_error(plot_model_results(toy_run(), "MULTI"), "must be a list")
  expect_error(plot_model_results(list(n.yr = 2, results = list(toy_run()), inoc = 1:2, n.max = 1),
                                  "MULTI"), "has 1 season")
})
