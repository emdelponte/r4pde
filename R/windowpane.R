#' Window Pane for Epidemiological Analysis
#'
#' This function calculates summary statistics within specified windows around a given end date
#' in a dataset, facilitating epidemiological analysis. It allows backward or forward window
#' calculations based on a user-defined variable and window lengths.
#'
#' @param data A dataframe containing the input data.
#' @param end_date_col The name of the column representing the end date.
#' @param date_col The name of the column representing the date variable.
#' @param variable The name of the column for which summary statistics are calculated.
#' @param summary_type The type of summary to calculate. Options are "mean", "sum",
#'   or "above_threshold".
#' @param threshold Optional numeric value used when summary_type is "above_threshold".
#' @param window_lengths A vector of window lengths (in days) for the calculations.
#' @param direction The direction of the window. Options are "backward" (default) or "forward".
#' @param group_by_cols Optional vector of column names for grouping the data.
#' @param date_format The format of the date columns. Default is "%Y-%m-%d".
#'
#' @return A dataframe with the calculated summary values for each window.
#'
#' @examples
#' # Example usage
#' data <- data.frame(
#'   study = rep(1, 10),
#'   date = as.Date("2023-01-01") + 0:9,
#'   value = rnorm(10),
#'   end_date = rep(as.Date("2023-01-10"), 10)
#' )
#' windowpane(data, end_date_col = end_date, date_col = date,
#'            variable = value, summary_type = "mean",
#'            window_lengths = c(3, 5), direction = "backward")
#'
#' @import dplyr tidyr
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

    results <- tibble(
      start_day = numeric(),
      end_day = numeric(),
      value = numeric()
    )

    # Window calculations
    for (window_size in window_lengths) {
      for (j in 0:(max(window_lengths) - window_size)) {

        if (direction == "backward") {
          window_end_date <- end_date - j
          window_start_date <- window_end_date - (window_size - 1)
        } else {
          window_start_date <- start_date + j
          window_end_date <- window_start_date + (window_size - 1)
        }

        # Filter data for the current window
        window_data <- subset_data %>%
          filter((!!date_col) >= window_start_date & (!!date_col) <= window_end_date)

        # Calculate summary value
        if (nrow(window_data) == 0) {
          value <- 0
        } else {
          if (summary_type == "mean") {
            value <- mean(pull(window_data, !!variable), na.rm = TRUE)
          } else if (summary_type == "sum") {
            value <- sum(pull(window_data, !!variable), na.rm = TRUE)
          } else if (summary_type == "above_threshold" & !is.null(threshold)) {
            value <- sum(pull(window_data, !!variable) > threshold, na.rm = TRUE)
          } else {
            value <- NA
          }
        }

        # Save results
        results <- bind_rows(results, tibble(
          start_day = ifelse(direction == "backward", -j, j),
          end_day = ifelse(direction == "backward", -(j + window_size - 1), j + window_size - 1),
          value = value
        ))
      }
    }

    # Pivot to wide format
    results_wide <- results %>%
      unite("column_name", start_day, end_day, sep = "_") %>%
      pivot_wider(names_from = column_name, values_from = value)

    # Identify grouping columns present in results_wide
    existing_group_cols <- intersect(names(results_wide), c(group_by_cols, as_label(end_date_col)))

    # Add the variable name as a prefix to all non-grouping columns
    results_wide <- results_wide %>%
      rename_with(~ paste0(as_label(variable), "_", .), .cols = -all_of(existing_group_cols))

    # Add grouping columns and end_date column
    if (!is.null(group_by_cols)) {
      for (col in group_by_cols) {
        results_wide[[col]] <- current_group[[col]]
      }
    }
    results_wide[[as_label(end_date_col)]] <- as.character(end_date)

    # Store the results for the current group
    results_list[[i]] <- results_wide
  }

  # Combine all results into a single dataframe
  final_results <- bind_rows(results_list)

  return(final_results)
}
