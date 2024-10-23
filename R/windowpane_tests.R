#' Windowpane Tests for Correlation Analysis
#'
#' This function performs bootstrapped correlation analysis for multiple predictors against a response variable.
#' It applies the Simes method for global significance testing and calculates individual correlations,
#' p-values, and bootstrap statistics.
#'
#' @param data A data frame containing the predictors and the response variable.
#' @param response_var A string representing the name of the response variable in the data frame.
#' @param corr_type A string specifying the correlation method to use; options are "spearman" (default),
#'        "pearson", or "kendall".
#' @param R An integer indicating the number of bootstrap replications. Default is 1000.
#' @param global_alpha A numeric value representing the global alpha level for the Simes correction.
#'        Default is 0.05.
#' @param individual_alpha A numeric value for the individual alpha threshold for testing individual predictors.
#'        Default is 0.005.
#'
#' @details The function calculates correlations between the response variable and each predictor in the
#' data frame, using bootstrapping to generate mean, standard deviation, and median estimates of the
#' correlation. The Simes correction is applied to control for multiple testing, providing a global
#' p-value (\emph{Pg}). The function also returns the maximum observed correlation.
#'
#' @return A list containing the following elements:
#' \item{results}{A data frame with columns: \code{variable}, \code{correlation}, \code{p_value},
#' \code{mean_corr}, \code{sd_corr}, \code{median_corr}, \code{rank}, \code{simes_threshold},
#' \code{significant_simes}, and \code{individual_significant}.}
#' \item{summary_table}{A data frame summarizing the global p-value (\emph{Pg}) and maximum correlation.}
#' \item{global_significant}{A logical value indicating whether the global test is significant.}
#'
#' @examples
#' # Example usage
#' data <- data.frame(x1 = rnorm(100), x2 = rnorm(100), y = rnorm(100))
#' result <- windowpane_tests(data, "y", corr_type = "spearman", R = 500)
#' result$results
#' result$summary_table
#'
#' @import boot
#' @import dplyr
#' @export
windowpane_tests <- function(data, response_var, corr_type = "spearman", R = 1000, global_alpha = 0.05, individual_alpha = 0.005) {
  # Define predictors and response
  predictors <- setdiff(names(data), response_var)
  response <- data[[response_var]]

  # Ensure predictors are numeric
  data[predictors] <- lapply(data[predictors], as.numeric)

  # Initialize the results data frame
  results <- data.frame(variable = predictors,
                        correlation = NA,
                        p_value = NA,
                        mean_corr = NA,
                        sd_corr = NA,
                        median_corr = NA)

  # Define the internal function for bootstrapping
  calc_correlation <- function(data, indices, var, response_var) {
    # Subset data for bootstrap sample
    sample_data <- data[indices, ]
    var_data <- sample_data[[var]]
    response_data <- sample_data[[response_var]]

    # Calculate correlation for the given sample
    corr_result <- cor.test(var_data, response_data,
                            method = corr_type, exact = FALSE)

    # Return the correlation estimate and p-value
    return(c(as.numeric(corr_result$estimate), corr_result$p.value))
  }

  # Loop through each predictor
  for (var in predictors) {
    var_data <- data[[var]]

    # Ensure the variable is numeric and not constant
    if (is.numeric(var_data) && length(unique(var_data)) > 1) {
      # Check if enough complete cases exist
      complete_cases <- complete.cases(var_data, response)
      if (sum(complete_cases) > 5) {  # Minimum sample size of 5
        # Prepare data for bootstrapping
        data_boot <- data.frame(var_data = var_data[complete_cases], response = response[complete_cases])

        # Run the bootstrap (using boot package)
        boot_result <- boot(data = data_boot, statistic = function(data, indices) {
          calc_correlation(data, indices, "var_data", "response")
        }, R = R)

        # Prepare bootstrap summary
        bootstrap_df <- as.data.frame(boot_result$t)
        colnames(bootstrap_df) <- c("correlation", "p_value")

        # Calculate mean, standard deviation, and median of the correlation estimates
        mean_corr <- mean(bootstrap_df$correlation, na.rm = TRUE)
        sd_corr <- sd(bootstrap_df$correlation, na.rm = TRUE)
        median_corr <- median(bootstrap_df$correlation, na.rm = TRUE)

        # Extract the initial correlation and p-value from the bootstrap
        corr <- boot_result$t0[1]
        p_value <- boot_result$t0[2]

        # Store results
        results[results$variable == var, c('correlation', 'p_value', 'mean_corr', 'sd_corr', 'median_corr')] <-
          c(corr, p_value, mean_corr, sd_corr, median_corr)
      }
    }
  }

  # Apply Simes method to adjust for multiple testing
  results <- results %>%
    arrange(p_value) %>%
    mutate(rank = row_number(),
           m = n(),
           simes_threshold = global_alpha * rank / m,
           significant_simes = p_value <= simes_threshold,
           individual_significant = p_value <= individual_alpha)  # Use individual_alpha = 0.005

  # Calculate the global p-value (Pg) as the minimum of the Simes-adjusted p-values
  Pg <- min(results$p_value / (results$rank / results$m), na.rm = TRUE)

  # Determine global significance
  global_significant <- Pg < global_alpha

  # Find the maximum correlation
  max_correlation <- max(results$correlation, na.rm = TRUE)

  # Add global Pg and max correlation as a separate row
  summary_table <- data.frame(
    Metric = c("Global P-value (Pg)", "Max Correlation"),
    Value = c(Pg, max_correlation)
  )

  # Return the results data frame and summary table
  return(list(
    results = results,
    summary_table = summary_table,
    global_significant = global_significant
  ))
}
