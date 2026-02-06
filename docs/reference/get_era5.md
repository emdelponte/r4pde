# Fetch ERA5 Data from Open-Meteo for Multiple Locations

This function downloads ERA5 historical weather data from the Open-Meteo
API for specified locations and dates. It mimics the functionality of
`get_nasapower` by fetching weather variables and returning a
daily-aggregated data frame. It uses the Open-Meteo Archive API and
aggregates hourly data to daily values to ensure all variables
(including relative humidity and dew point) are available.

## Usage

``` r
get_era5(
  data,
  days_around,
  date_col,
  study_col = "study",
  pars = c("temperature_2m", "relative_humidity_2m", "precipitation", "dewpoint_2m"),
  models = NULL
)
```

## Arguments

- data:

  A data frame containing the input data, including columns for
  `latitude`, `longitude`, and the date column specified by `date_col`.

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

  A character vector specifying the weather variables to fetch from
  Open-Meteo. These are mapped to `nasapower` equivalents:
  "temperature_2m" (T2M), "relative_humidity_2m" (RH2M), "precipitation"
  (PRECTOTCORR), etc.

- models:

  A character string specifying the reanalysis model(s) to use. If
  `NULL` (default), Open-Meteo uses its "best_match" logic (ERA5,
  ERA5-Land, etc.).

## Value

A data frame (tibble) with daily weather data. Columns are named to
match `nasapower` conventions: `date`, `T2M`, `T2M_MAX`, `T2M_MIN`,
`RH2M`, `PRECTOTCORR`, `T2MDEW`, `longitude`, `latitude`, and `study`.

## Details

Since the Open-Meteo Archive API `daily` endpoint does not provide all
variables (like mean relative humidity) directly, this function fetches
`hourly` data and performs the daily aggregation manually:

- `T2M`: Mean of hourly `temperature_2m`

- `T2M_MAX`: Max of hourly `temperature_2m`

- `T2M_MIN`: Min of hourly `temperature_2m`

- `RH2M`: Mean of hourly `relative_humidity_2m`

- `PRECTOTCORR`: Sum of hourly `precipitation`

- `T2MDEW`: Mean of hourly `dewpoint_2m`

## See also

Other Disease modeling:
[`get_brdwgd()`](https://emdelponte.github.io/r4pde/reference/get_brdwgd.md),
[`get_nasapower()`](https://emdelponte.github.io/r4pde/reference/get_nasapower.md),
[`windowpane()`](https://emdelponte.github.io/r4pde/reference/windowpane.md)
