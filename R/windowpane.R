#' Windowpane Analysis of Epidemiological Data
#'
#' This function calculates rolling statistics (e.g., mean, sum) within defined time windows
#' for a specified variable in epidemiological data.
#'
#' @param data A dataframe containing the data.
#' @param end_date_col A character string specifying the column name for the end date.
#' @param date_col A character string specifying the column name for the date.
#' @param variable A character string specifying the name of the variable to analyze.
#' @param study_col A character string specifying the column name for study IDs.
#' @param summary_type A character string specifying the type of summary ("mean", "sum", "above_threshold").
#' @param threshold A numeric value for the threshold (used only if summary_type is "above_threshold").
#' @param window_lengths A numeric vector specifying the window lengths to use (default is c(2, 4, 8, 16)).
#' @param direction A character string specifying the direction of the rolling window ("backward" or "forward").
#'
#' @return A dataframe with the results of the windowpane analysis in wide format.
#' @import dplyr tidyr lubridate progress
#' @examples
#' # Simulated data example
#' set.seed(123)
#' n <- 100
#' sim_data <- data.frame(
#'   study = sample(1:5, n, replace = TRUE),
#'   heading = as.Date('2024-01-01') + sample(0:30, n, replace = TRUE),
#'   YYYYMMDD = as.Date('2024-01-01') + sample(0:30, n, replace = TRUE),
#'   T2M = runif(n, 15, 30)
#' )
#' result <- windowpane(
#'   data = sim_data,
#'   end_date_col = "heading",
#'   date_col = "YYYYMMDD",
#'   variable = "T2M",
#'   study_col = "study",
#'   summary_type = "mean",
#'   window_lengths = c(2, 4, 8, 16),
#'   direction = "backward"
#' )
#' print(result)
#'
#' @export
windowpane <- function(data,
                       end_date_col,
                       date_col,
                       variable,
                       study_col,
                       summary_type,
                       threshold = NULL,
                       window_lengths = c(2, 4, 8, 16),
                       direction = "backward") {
  variable_name <- deparse(substitute(variable))

  results_list <- list()  # Initialize list to store results

  # Unique combinations of study and end_date
  unique_combinations <- data %>%
    distinct(!!sym(study_col), !!sym(end_date_col))

  # Calculate total steps accurately based on window_lengths
  total_steps <- 0
  for (window_size in window_lengths) {
    total_steps <- total_steps + (max(window_lengths) - window_size + 1)
  }
  total_steps <- total_steps * nrow(unique_combinations)

  # Initialize the progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent in :elapsed seconds",
    total = total_steps,
    clear = FALSE, width = 60
  )

  # Iterate through each combination of study and end_date
  for (i in 1:nrow(unique_combinations)) {
    study_id <- unique_combinations[[study_col]][i]
    end_date <- unique_combinations[[end_date_col]][i]

    subset_data <- data %>%
      filter(!!sym(study_col) == study_id, !!sym(end_date_col) == end_date)

    start_date <- end_date - days(max(window_lengths))

    results <- tibble(
      start_day = numeric(),
      end_day = numeric(),
      value = numeric()
    )

    # Window calculations
    for (window_size in window_lengths) {
      for (j in 0:(max(window_lengths) - window_size)) {

        if (direction == "backward") {
          window_end_date <- end_date - days(j)
          window_start_date <- window_end_date - days(window_size - 1)
        } else {
          window_start_date <- start_date + days(j)
          window_end_date <- window_start_date + days(window_size - 1)
        }

        # Filter data for the current window
        window_data <- subset_data %>%
          filter(!!sym(date_col) >= window_start_date & !!sym(date_col) <= window_end_date)

        # Calculate summary value based on the summary_type
        if (summary_type == "mean") {
          value <- mean(window_data[[variable]], na.rm = TRUE)
        } else if (summary_type == "sum") {
          value <- sum(window_data[[variable]], na.rm = TRUE)
        } else if (summary_type == "above_threshold" & !is.null(threshold)) {
          value <- sum(window_data[[variable]] > threshold, na.rm = TRUE)
        } else {
          value <- NA
        }

        # Save results
        results <- bind_rows(results, tibble(
          start_day = ifelse(direction == "backward", -j, j),
          end_day = ifelse(direction == "backward", -(j + window_size - 1), j + window_size - 1),
          value = value
        ))

        # Update progress bar after each iteration
        pb$tick()
      }
    }

    # Pivot to wide format
    results_wide <- results %>%
      unite("column_name", start_day, end_day, sep = "_") %>%
      pivot_wider(names_from = column_name, values_from = value)

    # Append variable name to column names
    colnames(results_wide) <- paste0(variable_name, "_", colnames(results_wide))

    # Convert 'end_date' to Date format if it is numeric
    if (is.numeric(end_date)) {
      end_date <- as.Date(end_date, origin = "1970-01-01")
    }

    # Add end_date and study information
    results_wide <- results_wide %>%
      mutate(!!end_date_col := end_date,
             !!study_col := study_id) %>%
      relocate(!!sym(end_date_col), .before = everything()) %>%
      relocate(!!sym(study_col), .before = everything())

    results_list[[paste0(study_id, "_", end_date)]] <- results_wide
  }

  # Combine all results into a single dataframe
  final_results <- bind_rows(results_list)

  return(final_results)
}






