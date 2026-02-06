# Fetch NASA POWER Data for Multiple Locations with a Progress Bar

This function downloads daily NASA POWER data for specified weather
variables over a specified number of days around a given date column for
multiple locations. It includes a progress bar to show the download
progress.

## Usage

``` r
get_nasapower(
  data,
  days_around,
  date_col,
  study_col = "study",
  pars = c("T2M", "RH2M", "PRECTOTCORR", "T2M_MAX", "T2M_MIN", "T2MDEW"),
  direction = "both"
)
```

## Arguments

- data:

  A data frame containing the input data, including columns for
  latitude, longitude, and the date column.

- days_around:

  An integer specifying the number of days before and after the date in
  the date column to download data.

- date_col:

  A character string specifying the name of the date column in the data
  frame.

- study_col:

  A character string specifying the name of the column containing the
  study identifier in the input data frame (default: "study").

- pars:

  A character vector specifying the weather variables to fetch from NASA
  POWER (default: c("T2M", "RH2M", "PRECTOTCORR", "T2M_MAX", "T2M_MIN",
  "T2MDEW")).

- direction:

  Character string specifying the direction of the date range relative
  to the reference date. Options are "both" (default), "back", or
  "forth". "back" retrieves data from `date - days_around` to `date`.
  "forth" retrieves data from `date` to `date + days_around`. "both"
  retrieves data from `date - days_around` to `date + days_around`.

## Value

A data frame with the downloaded weather data from NASA POWER, combined
for all specified locations. Includes a variable `study` indicating the
study identifier from the input data. Returns an empty data frame if no
data is retrieved.

## Details

The function uses the `get_power` function from the `nasapower` package
to fetch weather data for a range of dates around the specified date
column for each location. A progress bar is shown during the data
download process, and the results are combined into a single data frame.

## See also

Other Disease modeling:
[`get_brdwgd()`](https://emdelponte.github.io/r4pde/reference/get_brdwgd.md),
[`get_era5()`](https://emdelponte.github.io/r4pde/reference/get_era5.md),
[`windowpane()`](https://emdelponte.github.io/r4pde/reference/windowpane.md)
