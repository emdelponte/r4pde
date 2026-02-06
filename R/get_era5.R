#' Fetch ERA5 Data from Open-Meteo for Multiple Locations
#'
#' This function downloads ERA5 historical weather data from the Open-Meteo API
#' for specified locations and dates. It mimics the functionality of `get_nasapower`
#' by fetching weather variables and returning a daily-aggregated data frame.
#' It uses the Open-Meteo Archive API and aggregates hourly data to daily values
#' to ensure all variables (including relative humidity and dew point) are available.
#'
#' @param data A data frame containing the input data, including columns for `latitude`, `longitude`,
#' and the date column specified by `date_col`.
#' @param days_around An integer specifying the number of days before and after the date in the date column
#' to download data.
#' @param date_col A character string specifying the name of the date column in the data frame.
#' @param study_col A character string specifying the name of the column containing the study identifier
#' in the input data frame (default: "study").
#' @param pars A character vector specifying the weather variables to fetch from Open-Meteo.
#' These are mapped to `nasapower` equivalents:
#' "temperature_2m" (T2M), "relative_humidity_2m" (RH2M), "precipitation" (PRECTOTCORR), etc.
#' @param models A character string specifying the reanalysis model(s) to use. If `NULL` (default),
#' Open-Meteo uses its "best_match" logic (ERA5, ERA5-Land, etc.).
#' @param direction Character string specifying the direction of the date range relative to the reference date.
#' Options are "both" (default), "back", or "forth".
#' "back" retrieves data from `date - days_around` to `date`.
#' "forth" retrieves data from `date` to `date + days_around`.
#' "both" retrieves data from `date - days_around` to `date + days_around`.
#'
#' @return A data frame (tibble) with daily weather data. Columns are named to match
#' `nasapower` conventions: `date`, `T2M`, `T2M_MAX`, `T2M_MIN`, `RH2M`, `PRECTOTCORR`,
#' `T2MDEW`, `longitude`, `latitude`, and `study`.
#'
#' @details
#' Since the Open-Meteo Archive API `daily` endpoint does not provide all variables
#' (like mean relative humidity) directly, this function fetches `hourly` data and
#' performs the daily aggregation manually:
#' - `T2M`: Mean of hourly `temperature_2m`
#' - `T2M_MAX`: Max of hourly `temperature_2m`
#' - `T2M_MIN`: Min of hourly `temperature_2m`
#' - `RH2M`: Mean of hourly `relative_humidity_2m`
#' - `PRECTOTCORR`: Sum of hourly `precipitation`
#' - `T2MDEW`: Mean of hourly `dewpoint_2m`
#'
#' @importFrom dplyr mutate group_by summarize ungroup as_tibble select everything
#' @importFrom lubridate days as_date
#' @importFrom purrr map_dfr
#' @importFrom progress progress_bar
#' @importFrom httr GET content
#' @importFrom jsonlite fromJSON
#' @export
#' @family Disease modeling
get_era5 <- function(data, days_around, date_col, study_col = "study",
                     pars = c(
                       "temperature_2m", "relative_humidity_2m",
                       "precipitation", "dewpoint_2m"
                     ),
                     models = NULL, direction = "both") {
  # Helper function to fetch and aggregate data for a single location
  get_single_era5 <- function(data_row, days_around, date_col, study_col, pars, models, direction) {
    # Dynamically extract the date column
    given_date <- as.Date(data_row[[date_col]])

    # Calculate start and end dates
    # Calculate start and end dates
    if (direction == "back") {
      start_date <- given_date - lubridate::days(days_around)
      end_date <- given_date
    } else if (direction == "forth") {
      start_date <- given_date
      end_date <- given_date + lubridate::days(days_around)
    } else {
      # "both" or default
      start_date <- given_date - lubridate::days(days_around)
      end_date <- given_date + lubridate::days(days_around)
    }

    # Open-Meteo Archive API URL
    url <- "https://archive-api.open-meteo.com/v1/archive"

    # Identify coordinate columns (flexible matching)
    lon_col <- if ("longitude" %in% names(data_row)) "longitude" else if ("lon" %in% names(data_row)) "lon" else stop("Longitude column not found.")
    lat_col <- if ("latitude" %in% names(data_row)) "latitude" else if ("lat" %in% names(data_row)) "lat" else stop("Latitude column not found.")

    # Query parameters - requesting HOURLY to ensure we can calculate all means
    query <- list(
      latitude = data_row[[lat_col]],
      longitude = data_row[[lon_col]],
      start_date = as.character(start_date),
      end_date = as.character(end_date),
      hourly = paste(pars, collapse = ","),
      timezone = "auto"
    )

    if (!is.null(models)) query$models <- models

    # Fetch data
    res <- tryCatch(
      {
        httr::GET(url, query = query)
      },
      error = function(e) {
        warning(paste("Failed to fetch data for location:", data_row[[lat_col]], data_row[[lon_col]], "-", e$message))
        return(NULL)
      }
    )

    if (is.null(res)) {
      return(data.frame())
    }

    if (httr::status_code(res) != 200) {
      error_msg <- tryCatch(
        {
          jsonlite::fromJSON(httr::content(res, as = "text"))$reason
        },
        error = function(e) "Unknown error"
      )

      warning(paste("Open-Meteo API error", httr::status_code(res), " (", error_msg, ") for location:", data_row[[lat_col]], data_row[[lon_col]]))
      return(data.frame())
    }

    content_json <- httr::content(res, as = "text", encoding = "UTF-8")
    json_data <- jsonlite::fromJSON(content_json)

    if (is.null(json_data$hourly)) {
      return(data.frame())
    }

    # Convert hourly results to a data frame
    hourly_df <- as.data.frame(json_data$hourly)
    hourly_df$date <- as.Date(hourly_df$time)

    # Aggregate to daily
    # We check available columns first to avoid issues with names(.)
    avail_cols <- names(hourly_df)

    daily_df <- hourly_df %>%
      dplyr::group_by(.data$date) %>%
      dplyr::summarize(
        T2M = if ("temperature_2m" %in% avail_cols) mean(.data$temperature_2m, na.rm = TRUE) else NA_real_,
        T2M_MAX = if ("temperature_2m" %in% avail_cols) max(.data$temperature_2m, na.rm = TRUE) else NA_real_,
        T2M_MIN = if ("temperature_2m" %in% avail_cols) min(.data$temperature_2m, na.rm = TRUE) else NA_real_,
        RH2M = if ("relative_humidity_2m" %in% avail_cols) mean(.data$relative_humidity_2m, na.rm = TRUE) else NA_real_,
        PRECTOTCORR = if ("precipitation" %in% avail_cols) sum(.data$precipitation, na.rm = TRUE) else NA_real_,
        T2MDEW = if ("dewpoint_2m" %in% avail_cols) mean(.data$dewpoint_2m, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )

    # Add metadata
    daily_df$longitude <- json_data$longitude
    daily_df$latitude <- json_data$latitude
    daily_df$study <- data_row[[study_col]]

    return(daily_df)
  }

  # Initialize a progress bar
  pb <- progress_bar$new(
    format = "  Downloading ERA5 [:bar] :percent in :elapsed",
    total = nrow(data), clear = FALSE, width = 60
  )

  # Use map_dfr with progress bar
  results <- purrr::map_dfr(1:nrow(data), function(i) {
    pb$tick()
    get_single_era5(data[i, , drop = FALSE], days_around, date_col, study_col, pars, models, direction)
  })

  if (nrow(results) == 0) {
    return(dplyr::as_tibble(data.frame()))
  }

  return(dplyr::as_tibble(results))
}
