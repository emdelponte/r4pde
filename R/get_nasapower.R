#' Fetch NASA POWER Data for Multiple Locations with a Progress Bar
#'
#' This function downloads daily NASA POWER data for specified weather variables over a specified number of days
#' around a given date column for multiple locations. It includes a progress bar to show the download progress.
#'
#' @param data A data frame containing the input data, including columns for latitude, longitude, study identifier,
#' and the date column.
#' @param days_around An integer specifying the number of days before and after the date in the date column
#' to download data.
#' @param date_col A character string specifying the name of the date column in the data frame.
#' @param pars A character vector specifying the weather variables to fetch from NASA POWER
#' (default: c("T2M", "RH2M", "PRECTOTCORR", "T2M_MAX", "T2M_MIN", "T2MDEW")).
#'
#' @return A data frame with the downloaded weather data from NASA POWER, combined for all specified locations.
#' Includes a new variable `study` indicating the study identifier from the input data.
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
#'
#' @examples
#' \dontrun{
#' # Sample data
#' data <- data.frame(
#'   study = c("Study1", "Study2"),
#'   latitude = c(-15.78, -20.45),
#'   longitude = c(-47.93, -54.82),
#'   heading = c("2022-05-10", "2022-05-15")
#' )
#'
#' # Fetch weather data with a progress bar
#' results <- get_nasapower(data, days_around = 28, date_col = "heading")
#' }
#'
#' @export
#' @family Disease modeling
get_nasapower <- function(data, days_around, date_col, pars = c("T2M", "RH2M", "PRECTOTCORR", "T2M_MAX", "T2M_MIN", "T2MDEW")) {
  # Helper function to fetch data for a single location
  get_nasa_data <- function(data, index, days_around, date_col, pars) {
    # Dynamically extract the date column
    given_date <- as.Date(data[[date_col]][index])

    # Calculate start and end dates based on the specified given date and days_around
    start_date <- given_date - days(days_around)
    end_date <- given_date + days(days_around)

    # Fetch data using get_power
    response <- get_power(
      community = "ag",
      temporal_api = "daily",
      dates = c(start_date, end_date),
      lonlat = c(data$longitude[index], data$latitude[index]),
      pars = pars
    )

    if (nrow(response) > 0) {
      response <- response %>% mutate(study = data$study[index])
    }

    return(response)
  }

  # Initialize a progress bar
  pb <- progress_bar$new(
    format = "  Downloading [:bar] :percent in :elapsed",
    total = nrow(data), clear = FALSE, width = 60
  )

  # Use map_dfr with progress bar
  results <- map_dfr(1:nrow(data), function(index) {
    pb$tick()  # Update the progress bar
    get_nasa_data(data, index, days_around, date_col, pars)
  })

  return(results)
}
