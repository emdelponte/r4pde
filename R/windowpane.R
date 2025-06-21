#' Window Pane for Epidemiological Analysis
#'
#' This function calculates summary statistics within specified windows around a given end date
#' in a dataset, facilitating epidemiological analysis. It allows backward, forward, or both
#' directions of window calculations based on a user-defined variable and window lengths.
#'
#' @param data A data frame containing the input data.
#' @param end_date_col A string specifying the name of the column representing the end date.
#' @param date_col A string specifying the name of the column representing the date variable.
#' @param variable A string specifying the name of the column for which summary statistics are calculated.
#' @param summary_type A string specifying the type of summary to calculate. Options are "mean", "sum",
#'   "above_threshold", or "below_threshold".
#' @param threshold Optional numeric value used when \code{summary_type} is "above_threshold" or "below_threshold".
#' @param window_lengths A numeric vector specifying the window lengths (in days) for the calculations.
#' @param direction A string specifying the direction of the window. Options are "backward" (default),
#'   "forward", or "both".
#' @param group_by_cols Optional vector of strings specifying column names for grouping the data.
#' @param date_format A string specifying the format of the date columns. Default is "%Y-%m-%d".
#'
#' @return A data frame with the calculated summary values for each window.
#'
#' @import tidyr
#' @export
windowpane <- function(data,
                       end_date_col,
                       date_col,
                       variable,
                       summary_type,
                       threshold = NULL,
                       window_lengths,
                       direction = "backward",
                       group_by_cols = NULL,
                       date_format = "%Y-%m-%d") {

  # Convert columns to symbols for tidy evaluation
  end_date_col <- enquo(end_date_col)
  date_col <- enquo(date_col)
  variable <- enquo(variable)

  # Prepare data and convert date columns
  data <- data %>%
    mutate(
      {{ date_col }} := as.Date(as.character({{ date_col }}), format = date_format),
      {{ end_date_col }} := as.Date(as.character({{ end_date_col }}), format = date_format)
    )

  # Handle grouping
  if (!is.null(group_by_cols)) {
    grouped_data <- data %>% group_by(across(all_of(group_by_cols)))
  } else {
    grouped_data <- data
  }

  # Get unique combinations of grouping variables and end_date_col
  unique_combinations <- grouped_data %>%
    distinct(across(c(group_by_cols, !!end_date_col))) %>%
    ungroup()

  results_list <- list()

  # Iterate over each unique combination
  for (i in seq_len(nrow(unique_combinations))) {
    current_group <- unique_combinations[i, ]
    end_date <- pull(current_group, !!end_date_col)

    # Filter data for the current group
    subset_data <- data
    if (!is.null(group_by_cols)) {
      for (col in group_by_cols) {
        subset_data <- subset_data %>% filter(.data[[col]] == current_group[[col]])
      }
    }
    subset_data <- subset_data %>% filter((!!end_date_col) == end_date)

    start_date <- end_date - max(window_lengths)

    # Temporary list to store results for each window size
    temp_results <- list()

    # Backward window calculations
    if (direction == "backward" || direction == "both") {
      for (window_size in window_lengths) {
        window_results <- tibble(
          start_day = numeric(),
          end_day = numeric(),
          value = numeric()
        )

        for (j in 0:(max(window_lengths) - window_size)) {
          window_end_date <- end_date - j
          window_start_date <- window_end_date - (window_size - 1)

          # Filter data for the current window
          window_data <- subset_data %>%
            filter((!!date_col) >= window_start_date & (!!date_col) <= window_end_date)

          # Calculate summary value
          value <- switch(
            summary_type,
            "mean" = mean(pull(window_data, !!variable), na.rm = TRUE),
            "sum" = sum(pull(window_data, !!variable), na.rm = TRUE),
            "above_threshold" = if (!is.null(threshold)) sum(pull(window_data, !!variable) > threshold, na.rm = TRUE),
            "below_threshold" = if (!is.null(threshold)) sum(pull(window_data, !!variable) < threshold, na.rm = TRUE),
            NA
          )

          # Save results
          window_results <- bind_rows(window_results, tibble(
            start_day = -j,
            end_day = -(j + window_size - 1),
            value = value
          ))
        }

        # Pivot to wide format for current window size
        window_results_wide <- window_results %>%
          unite("column_name", start_day, end_day, sep = "_") %>%
          pivot_wider(names_from = column_name, values_from = value)

        # Add prefix to column names for backward direction
        window_results_wide <- window_results_wide %>%
          rename_with(~ paste0("bw_length", window_size, "_", as_label(variable), "_", summary_type, "_", .))

        # Store results in temporary list
        temp_results[[paste0("bw_", window_size)]] <- window_results_wide
      }
    }

    # Forward window calculations
    if (direction == "forward" || direction == "both") {
      for (window_size in window_lengths) {
        window_results <- tibble(
          start_day = numeric(),
          end_day = numeric(),
          value = numeric()
        )

        for (j in 0:(max(window_lengths) - window_size)) {
          window_start_date <- start_date + j
          window_end_date <- window_start_date + (window_size - 1)

          # Filter data for the current window
          window_data <- subset_data %>%
            filter((!!date_col) >= window_start_date & (!!date_col) <= window_end_date)

          # Calculate summary value
          value <- switch(
            summary_type,
            "mean" = mean(pull(window_data, !!variable), na.rm = TRUE),
            "sum" = sum(pull(window_data, !!variable), na.rm = TRUE),
            "above_threshold" = if (!is.null(threshold)) sum(pull(window_data, !!variable) > threshold, na.rm = TRUE),
            "below_threshold" = if (!is.null(threshold)) sum(pull(window_data, !!variable) < threshold, na.rm = TRUE),
            NA
          )

          # Save results
          window_results <- bind_rows(window_results, tibble(
            start_day = j,
            end_day = j + window_size - 1,
            value = value
          ))
        }

        # Pivot to wide format for current window size
        window_results_wide <- window_results %>%
          unite("column_name", start_day, end_day, sep = "_") %>%
          pivot_wider(names_from = column_name, values_from = value)

        # Add prefix to column names for forward direction
        window_results_wide <- window_results_wide %>%
          rename_with(~ paste0("fw_length", window_size, "_", as_label(variable), "_", summary_type, "_", .))

        # Store results in temporary list
        temp_results[[paste0("fw_", window_size)]] <- window_results_wide
      }
    }

    # Combine results for all window sizes for the current group
    combined_results <- bind_cols(temp_results)

    # Add grouping columns and end_date column
    if (!is.null(group_by_cols)) {
      for (col in group_by_cols) {
        combined_results[[col]] <- current_group[[col]]
      }
    }
    combined_results[[as_label(end_date_col)]] <- as.character(end_date)

    # Store the results for the current group
    results_list[[i]] <- combined_results
  }

  # Combine all results into a single dataframe
  final_results <- bind_rows(results_list)

  return(final_results)
}
