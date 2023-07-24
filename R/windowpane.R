#' Sliding Window Analysis for Time-Series Data
#'
#' This function performs a sliding window analysis over a given data frame for the specified variable.
#' It calculates the mean, sum, or count above a given threshold for the variable values in each window.
#' The window size is defined by the user-specified window lengths.
#'
#' @param data A data frame that contains the data for analysis. It must include a `YYYYMMDD` date column.
#' @param end_date A character string that specifies the end date for the analysis in the format 'YYYY-MM-DD'.
#' @param variable The name of the variable in the data frame on which to perform the analysis.
#' @param summary_type The type of summary statistic to compute for each window. Can be 'mean', 'sum', or 'above_threshold'.
#' @param threshold An optional numerical value. If `summary_type` is 'above_threshold', the function counts how many variable values in each window are above this threshold.
#' @param window_lengths A numerical vector that specifies the lengths of the windows for the analysis.
#'
#' @return A data frame in wide format where each column corresponds to a window of a certain length and start day, and contains the calculated summary statistic for that window.
#' The column names are in the format 'variable_startday_endday'.
#'
#' @importFrom dplyr filter summarize
#' @importFrom lubridate days
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' data <- tibble(YYYYMMDD = seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "day"),
#'                variable = rnorm(365))
#' windowpane(data, "2023-12-31", variable, "mean", NULL, c(7, 14, 21))
#' }

windowpane <- function(data,
                       end_date,
                       variable,
                       summary_type,
                       threshold = NULL,
                       window_lengths) {
  variable_name <- deparse(substitute(variable))
  data$date <- as.Date(data$YYYYMMDD)
  end_date <- as.Date(end_date)
  start_date <- end_date - days(max(window_lengths))

  results <- tibble(
    start_day = numeric(),
    end_day = numeric(),
    value = numeric()
  )
  for (window_size in window_lengths) {
    for (i in 0:(max(window_lengths) - window_size)) {
      window_end_date <- end_date - days(i)
      window_start_date <- window_end_date - days(window_size - 1)
      window_data <- data %>%
        filter(date >= window_start_date & date <= window_end_date)
      window_summary <- case_when(
        summary_type == "mean" ~ window_data %>%
          summarize(value = mean({{variable}}, na.rm = TRUE)),
        summary_type == "sum" ~ window_data %>%
          summarize(value = sum({{variable}}, na.rm = TRUE)),
        summary_type == "above_threshold" ~ window_data %>%
          summarize(value = sum({{variable}} > threshold, na.rm = TRUE))
      )

      # Save results
      results <- rbind(results, tibble(
        start_day = -i, # Start day is always 0 in this case
        end_day = -(i + window_size - 1), # End day is the window size minus 1
        value = window_summary$value
      ))
    }
  }

  # Pivot to wide format
  results_wide <- results %>%
    unite("column_name", start_day, end_day, sep = "_") %>%
    pivot_wider(names_from = column_name, values_from = value)

  # Adding the variable name to the column names
  colnames(results_wide) <- paste0(variable_name, "_", colnames(results_wide))

  return(results_wide)
}

