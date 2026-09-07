test_that("get_era5 validates input coordinates and column names without network call", {
  df_no_coords <- data.frame(
    date = as.Date("2021-01-01"),
    study = 1
  )

  expect_error(
    get_era5(df_no_coords, days_around = 1, date_col = "date"),
    "Longitude column not found"
  )

  df_no_lat <- data.frame(
    longitude = -45.0,
    date = as.Date("2021-01-01"),
    study = 1
  )

  expect_error(
    get_era5(df_no_lat, days_around = 1, date_col = "date"),
    "Latitude column not found"
  )
})

test_that("get_brdwgd validates coordinate columns and file path", {
  df_no_coords <- data.frame(
    date = as.Date("2021-01-01"),
    study = 1
  )

  expect_error(
    get_brdwgd(df_no_coords, days_around = 1, date_col = "date", path = "nonexistent_dir"),
    "Longitude column not found"
  )

  df_with_coords <- data.frame(
    longitude = -45.0,
    latitude = -20.0,
    date = as.Date("2021-01-01"),
    study = 1
  )

  expect_warning(
    res <- get_brdwgd(df_with_coords, days_around = 1, date_col = "date", path = "nonexistent_dir"),
    "No files for"
  )
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0)
})

test_that("get_nasapower validates coordinate columns without network call", {
  df_no_coords <- data.frame(
    date = as.Date("2021-01-01"),
    study = 1
  )

  expect_error(
    get_nasapower(df_no_coords, days_around = 1, date_col = "date"),
    "Longitude column not found"
  )

  df_no_lat <- data.frame(
    longitude = -45.0,
    date = as.Date("2021-01-01"),
    study = 1
  )

  expect_error(
    get_nasapower(df_no_lat, days_around = 1, date_col = "date"),
    "Latitude column not found"
  )
})
