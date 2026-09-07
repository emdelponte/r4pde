test_that("DSI calculates disease severity index correctly", {
  # Known example from documentation
  unit <- c(1, 2, 3, 4, 5, 6)
  class <- c(1, 2, 1, 2, 3, 1)
  max <- 3

  # Expected: sum of weights = (1*3 + 2*2 + 3*1) = 10
  # DSI = 10 / (6 * 3) * 100 = 55.55555...
  res <- DSI(unit, class, max)
  expect_equal(res, (10 / 18) * 100, tolerance = 1e-6)

  # Edge cases
  # All zero class
  expect_equal(DSI(1:4, rep(0, 4), max = 5), 0)
  # All max class
  expect_equal(DSI(1:4, rep(5, 4), max = 5), 100)
})

test_that("DSI2 calculates disease severity index from frequencies correctly", {
  # Example from documentation
  res <- DSI2(c(0, 1, 2, 3, 4), c(2, 0, 5, 0, 5), 4)
  # Weighted sum: 0*2 + 1*0 + 2*5 + 3*0 + 4*5 = 30
  # Total freq: 2 + 0 + 5 + 0 + 5 = 12
  # Max: 4 -> 30 / (12 * 4) * 100 = 62.5
  expect_equal(res, 62.5)

  # Edge case: All zeros
  expect_equal(DSI2(c(0, 1, 2), c(10, 0, 0), max = 2), 0)
  # Edge case: All max
  expect_equal(DSI2(c(0, 1, 2), c(0, 0, 10), max = 2), 100)
})

test_that("CompMuCens executes when interval package is available", {
  skip_if_not_installed("interval")

  inputData <- data.frame(
    treatment = c(rep("A", 10), rep("B", 10)),
    x = c(rep(2, 5), rep(3, 5), rep(4, 5), rep(5, 5))
  )
  res <- CompMuCens(
    dat = inputData,
    scale = c(0, 3, 6, 12, 25, 50, 75, 88, 94, 97, 100, 100),
    ckData = TRUE
  )
  expect_type(res, "list")
  expect_true(!is.null(res$Test_result))
})
