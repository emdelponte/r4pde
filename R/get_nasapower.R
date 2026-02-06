#' Fetch NASA POWER Data for Multiple Locations with a Progress Bar
#'
#' This function downloads daily NASA POWER data for specified weather variables over a specified number of days
#' around a given date column for multiple locations. It includes a progress bar to show the download progress.
#'
#' @param data A data frame containing the input data, including columns for latitude, longitude, and 
#' the date column.
#' @param days_around An integer specifying the number of days before and after the date in the date column
#' to download data.
#' @param date_col A character string specifying the name of the date column in the data frame.
#' @param study_col A character string specifying the name of the column containing the study identifier 
#' in the input data frame (default: "study").
#' @param pars A character vector specifying the weather variables to fetch from NASA POWER
#' (default: c("T2M", "RH2M", "PRECTOTCORR", "T2M_MAX", "T2M_MIN", "T2MDEW")).
#'
#' @return A data frame with the downloaded weather data from NASA POWER, combined for all specified locations.
#' Includes a variable `study` indicating the study identifier from the input data.
#' Returns an empty data frame if no data is retrieved.
#'
#' @details
#' The function uses the `get_power` function from the `nasapower` package to fetch weather data for a range of
#' dates around the specified date column for each location. A progress bar is shown during the data download
#' process, and the results are combined into a single data frame.
#'
#' @importFrom dplyr mutate bind_rows
#' @importFrom lubridate days
#' @importFrom nasapower get_power
#' @importFrom purrr map_dfr
#' @importFrom progress progress_bar
#' @export
#' @family Disease modeling
get_nasapower <- function(data, days_around, date_col, study_col = "study", 
                          pars = c("T2M", "RH2M", "PRECTOTCORR", "T2M_MAX", "T2M_MIN", "T2MDEW")) {
  # Helper function to fetch data for a single location
  get_nasa_data <- function(data, index, days_around, date_col, study_col, pars) {
    # Dynamically extract the date column
    given_date <- as.Date(data[[date_col]][index])

    # Calculate start and end dates based on the specified given date and days_around
    start_date <- given_date - days(days_around)
    end_date <- given_date + days(days_around)

    # Identify coordinate columns (flexible matching)
    lon_col <- if ("longitude" %in% names(data)) "longitude" else if ("lon" %in% names(data)) "lon" else stop("Longitude column not found.")
    lat_col <- if ("latitude" %in% names(data)) "latitude" else if ("lat" %in% names(data)) "lat" else stop("Latitude column not found.")

    # Fetch data using get_power
    response <- get_power(
      community = "ag",
      temporal_api = "daily",
      dates = c(start_date, end_date),
      lonlat = c(data[[lon_col]][index], data[[lat_col]][index]),
      pars = pars
    )

    if (nrow(response) > 0) {
      # Use the specified study_col to assign the study identifier
      response <- response %>% mutate(study = data[[study_col]][index])
    }

    return(response)
  }

  # Initialize a progress bar
  pb <- progress_bar$new(
    format = "  Downloading NASA POWER [:bar] :percent in :elapsed",
    total = nrow(data), clear = FALSE, width = 60
  )

  # Use map_dfr with progress bar
  results <- map_dfr(1:nrow(data), function(index) {
    pb$tick()  # Update the progress bar
    get_nasa_data(data, index, days_around, date_col, study_col, pars)
  })

  return(results)
}
