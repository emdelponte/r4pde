# Fetch BR-DWGD Data for Multiple Locations

This function extracts daily weather data from the Brazilian Daily
Weather Gridded Data (BR-DWGD) NetCDF files for specified locations and
dates. It mimics the functionality of `get_nasapower` and `get_era5` by
returning a daily-aggregated data frame with columns named according to
`nasapower` conventions.

## Usage

``` r
get_brdwgd(
  data,
  days_around,
  date_col,
  study_col = "study",
  path = "netcdf_files",
  vars = c("pr", "Tmax", "Tmin", "Rs", "RH")
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

- path:

  A character string specifying the directory where the BR-DWGD NetCDF
  files are stored.

- vars:

  A character vector specifying the weather variables to fetch. These
  are expected to be the prefix of the NetCDF files (e.g., "pr",
  "Tmax"). Default: `c("pr", "Tmax", "Tmin", "Rs", "RH")`.

## Value

A data frame (tibble) with daily weather data. Columns are named to
match `nasapower` conventions where possible: `date`, `PRECTOTCORR`
(from `pr`), `T2M_MAX` (from `Tmax`), `T2M_MIN` (from `Tmin`),
`ALLSKY_SFC_SW_DWN` (from `Rs`), `RH2M` (from `RH`), `longitude`,
`latitude`, and `study`.

## Details

The function requires the `terra` package to read NetCDF files and the
`tidyr` package for data reshaping. It expects NetCDF files to be named
starting with the variable name (e.g., `pr_*.nc`) and to contain the
variable with the same name.

## See also

Other Disease modeling:
[`get_era5()`](https://emdelponte.github.io/r4pde/reference/get_era5.md),
[`get_nasapower()`](https://emdelponte.github.io/r4pde/reference/get_nasapower.md),
[`windowpane()`](https://emdelponte.github.io/r4pde/reference/windowpane.md)
