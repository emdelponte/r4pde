#' Binary Power Law Analysis for Spatial Disease Patterns
#'
#' This function calculates the Binary Power Law (BPL) parameters for spatial disease patterns,
#' fits a linear model, and performs a hypothesis test for the slope.
#'
#' @param data A data frame containing the following columns:
#'   - \code{field}: The field identifier.
#'   - \code{n}: The number of observations in each quadrat.
#'   - \code{i}: The incidence count in each quadrat.
#'
#' @return A list containing the following elements:
#'   - \code{summary}: A data frame summarizing the input data by field, including total observations (\code{n_total}),
#'     mean incidence (\code{incidence_mean}), observed variance (\code{V}), and binomial variance (\code{Vbin}).
#'   - \code{model_summary}: A summary of the linear model fitted to the log-transformed variances.
#'   - \code{hypothesis_test}: The result of the hypothesis test for the slope being equal to 1.
#'   - \code{ln_Ap}: The intercept of the linear model, representing the natural logarithm of the parameter \( A_p \).
#'   - \code{slope}: The slope of the linear model.
#'
#' @details The function performs the following steps:
#'   1. Summarizes the data by field to calculate the total number of observations (\code{n_total}),
#'      mean incidence (\code{incidence_mean}), observed variance (\code{V}), and binomial variance (\code{Vbin}).
#'   2. Log-transforms the variances.
#'   3. Fits a linear model to the log-transformed variances.
#'   4. Tests the hypothesis that the slope of the linear model is equal to 1.
#' @examples
#' \donttest{
#' # Example usage with a sample data frame
#' result <- BPL(FHBWheat)
#' print(result$summary)
#' print(result$model_summary)
#' print(result$hypothesis_test)
#' print(paste("ln(Ap):", result$ln_Ap))
#' print(paste("Slope (b):", result$slope))
#' }
#' @export
#' @family Spatial analysis

BPL <- function(data) {
  # Summarize the data by field
  df_summary <- data %>%
    group_by(field) %>%
    mutate(p = i / n) %>%
    summarise(
      n_total = sum(n),
      incidence_mean = mean(p, na.rm = TRUE),
      V = var(p, na.rm = TRUE),
      Vbin = (mean(p, na.rm = TRUE) * (1 - mean(p, na.rm = TRUE))) / (mean(n, na.rm = TRUE))
    )

  # Log-transform the variances
  df_summary <- df_summary %>%
    mutate(
      log_V = log(V),
      log_Vbin = log(Vbin)
    )

  # Fit the linear model
  model <- lm(log_V ~ log_Vbin, data = df_summary)

  # Summary of the model
  model_summary <- summary(model)

  # Hypothesis test for slope = 1
  hypothesis_test <- linearHypothesis(model, "log_Vbin = 1")

  # Extract coefficients
  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  # Calculate ln(Ap)
  ln_Ap <- intercept

  return(list(
    summary = df_summary,
    model_summary = model_summary,
    hypothesis_test = hypothesis_test,
    ln_Ap = ln_Ap,
    slope = slope
  ))
}
