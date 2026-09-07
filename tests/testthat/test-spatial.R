test_that("AFSD detects foci structure and returns expected components", {
  df <- data.frame(
    x = c(1, 2, 3, 10, 11),
    y = c(1, 1, 1, 10, 10),
    i = c(1, 1, 1, 1, 1)
  )

  res <- AFSD(df)
  expect_type(res, "list")
  expect_length(res, 3)

  # Check cluster df (second element) has 2 distinct foci
  focus_df <- res[[2]]
  expect_equal(nrow(focus_df), 2)
  expect_equal(focus_df$size, c(3, 2))

  # Error on missing columns
  expect_error(AFSD(data.frame(x = 1, y = 1)), "must contain 'x', 'y', and 'i'")
  # Error on empty infected
  expect_error(AFSD(data.frame(x = 1, y = 1, i = 0)), "no rows where i == 1")
})

test_that("BPL fits binary power law and returns model parameters", {
  data("FHBWheat", package = "r4pde")
  res <- BPL(FHBWheat)

  expect_type(res, "list")
  expect_named(res, c("summary", "model_summary", "hypothesis_test", "ln_Ap", "slope"))

  # Slope should be close to 1 for FHBWheat
  expect_true(is.numeric(res$slope))
  expect_true(res$slope > 0.5 && res$slope < 1.5)
  expect_s3_class(res$model_summary, "summary.lm")
})

test_that("join_count calculates join count statistics correctly", {
  # 2x2 matrix with diagonal ones
  m <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
  jc <- join_count(m, verbose = FALSE)

  expect_type(jc, "list")
  expect_named(jc, c("observed", "expected", "sd", "zscore", "pattern"))
  expect_true(is.numeric(jc$observed$HD))
  expect_true(is.numeric(jc$observed$DD))
})

test_that("oruns_test calculates ordinary runs test and prints output", {
  # Sequence of runs: 1, 0, 1, 1, 0, 1, 0, 0, 1, 1
  # Runs: [1], [0], [1,1], [0], [1], [0,0], [1,1] -> 7 runs
  s <- c(1, 0, 1, 1, 0, 1, 0, 0, 1, 1)
  res <- oruns_test(s)

  expect_s3_class(res, "r4pde.oruns_test")
  expect_equal(res$U, 7)
  expect_true(res$EU > 0)
  expect_true(is.numeric(res$Z))
  expect_true(res$pvalue >= 0 && res$pvalue <= 1)
  expect_true(res$result %in% c("aggregation or clustering", "randomness"))

  # Print method output check
  expect_message(print(res), "Total Number of Runs")
})

test_that("windowpane aggregates variables over date windows", {
  df <- data.frame(
    date = as.Date("2021-01-01") + 0:19,
    end_date = as.Date("2021-01-20"),
    val = 1:20
  )

  wp <- windowpane(
    data = df,
    end_date_col = end_date,
    date_col = date,
    variable = val,
    summary_type = "mean",
    window_lengths = c(5, 10)
  )

  expect_s3_class(wp, "data.frame")
  expect_equal(nrow(wp), 1)
  expect_true(any(grepl("bw_length5", names(wp))))
  expect_true(any(grepl("bw_length10", names(wp))))
})
