library(testthat)

test_that("functional_contrast errors correctly", {
  dat <- simulate_dsp_data()
  
  # Error if control doesn't exist
  expect_error(
    functional_contrast(dat, time = "time", response = "severity", treatment = "treatment", reference = "NonExistent"),
    "Reference level 'NonExistent' not found"
  )
  
  # Error if control missing in a group
  dat_missing_control <- dat |> dplyr::filter(!(treatment == "Control" & rep == 1))
  expect_error(
    functional_contrast(dat_missing_control, time = "time", response = "severity", treatment = "treatment", reference = "Control", group = "rep"),
    "Reference treatment is missing in some groups"
  )
})

test_that("functional_contrast calculates DSP correctly", {
  # Simple manual data
  dat <- data.frame(
    time = c(1, 2, 3, 1, 2, 3),
    trt = c("Control", "Control", "Control", "A", "A", "A"),
    y = c(10, 20, 30, 5, 5, 10)
  )
  
  res <- functional_contrast(dat, time = "time", response = "y", treatment = "trt", reference = "Control")
  
  expect_s3_class(res, "functional_dsp")
  expect_true("DSP" %in% names(res))
  
  # For treatment A: (10-5), (20-5), (30-10) -> 5, 15, 20
  expect_equal(res$DSP, c(5, 15, 20))
  
  # Check class is preserved and it's a tibble
  expect_true(tibble::is_tibble(res))
})

test_that("functional_summary calculations are correct", {
  dat <- data.frame(
    time = c(0, 10, 20),
    trt = rep("A", 3),
    DSP = c(0, 10, 0) # Triangle with base 20, height 10. Area = 100
  )
  
  # Make it look like a functional_dsp object
  attr(dat, "dsp_vars") <- list(time = "time", treatment = "trt", response = "DSP", group = NULL)
  class(dat) <- c("functional_dsp", class(dat))
  
  summ <- functional_summary(dat)
  
  expect_true(tibble::is_tibble(summ))
  expect_equal(summ$protected_area, 100) # (10*10)/2 + (10*10)/2
  expect_equal(summ$max_suppression, 10)
  expect_equal(summ$time_max_suppression, 10)
  
  # Persistence with threshold 0 should be 20 (from 0 to 20, since > 0 is not true at boundaries? 
  # Wait, DSP > 0 is only true at time 10. 
  # So persistence > 0 might be 0 since only 1 point is > 0.
  # Let's test with a different threshold.
  
  dat2 <- data.frame(
    time = c(0, 5, 10, 15, 20),
    trt = rep("B", 5),
    DSP = c(0, 10, 10, 10, 0) 
  )
  attr(dat2, "dsp_vars") <- list(time = "time", treatment = "trt", response = "DSP", group = NULL)
  class(dat2) <- c("functional_dsp", class(dat2))
  
  summ2 <- functional_summary(dat2, threshold = 5)
  # >5 at 5, 10, 15. Persistence = 15 - 5 = 10.
  expect_equal(summ2$persistence, 10)
})

test_that("integration with plot and ranks works", {
  dat <- simulate_dsp_data()
  dsp <- functional_contrast(dat, time = "time", response = "severity", treatment = "treatment", reference = "Control", group = "rep")
  
  # Just group by treatment for summary to average out reps, or summarize per group
  summ <- functional_summary(dsp, group = "rep")
  
  ranks <- rank_dsp(summ)
  
  expect_true("rank_protected_area" %in% names(ranks))
  expect_true("average_rank" %in% names(ranks))
  
  # plot
  p1 <- plot_dsp(dsp)
  p2 <- plot_dsp_rank_heatmap(ranks, group = "rep")
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("functional_summary works with custom column names", {
  dat <- data.frame(
    days = c(0, 10, 20),
    fungicide = rep("A", 3),
    contrast = c(0, 10, 0)
  )
  
  summ <- functional_summary(dat, time = "days", response = "contrast", treatment = "fungicide")
  
  expect_true(tibble::is_tibble(summ))
  expect_equal(summ$protected_area, 100)
  expect_equal(summ$max_suppression, 10)
  expect_equal(summ$time_max_suppression, 10)
})
