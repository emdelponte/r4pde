
#' Fit Gradient Models to Data
#'
#' This function fits three gradient models (exponential, power, and modified power) to given data.
#' It then ranks the models based on their R-squared values and returns diagnostic plots for each model.
#'
#' @param data A dataframe containing the data, with columns "x" representing distances and
#' "Y" representing the corresponding measurements or counts.
#' @param C A constant to be used in the modified power model. Defaults to 1.
#'
#' @return A list containing:
#' \item{data}{The input data, which will include an additional column 'mod_x'.}
#' \item{results_table}{A table of the model parameters and R-squared values.}
#' \item{plot_exponential}{Diagnostic plot for the exponential model.}
#' \item{plot_power}{Diagnostic plot for the power model.}
#' \item{plot_modified_power}{Diagnostic plot for the modified power model.}
#' \item{plot_exponential_original}{Plot of the original data with the exponential model fit.}
#' \item{plot_power_original}{Plot of the original data with the power model fit.}
#' \item{plot_modified_power_original}{Plot of the original data with the modified power model fit.}
#'
#' @examples
#' x <- c(0.8, 1.6, 2.4, 3.2, 4, 7.2, 12, 15.2, 21.6, 28.8)
#' Y <- c(184.9, 113.3, 113.3, 64.1, 25, 8, 4.3, 2.5, 1, 0.8)
#' grad1 <- data.frame(x = x, Y = Y)
#' library(ggplot2)
#' mg <- fit_gradients(grad1, C = 0.4)
#' mg$plot_power_original +
#'   labs(title = "", x = "Distance from focus (m)", y = "Count of lesions")
#'
#' @export
#' @family Spatial analysis

fit_gradients <- function(data, C = 1) {
  # Ensure column names are as expected
  if (!all(c("x", "Y") %in% colnames(data))) {
    stop("The dataframe should have columns named 'x' and 'Y'")
  }

  # Ensure ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required. Please install it with install.packages('ggplot2').")
  }

  # Fit the models
  exponential <- lm(log(Y) ~ x, data = data)
  power <- lm(log(Y) ~ log(x), data = data)
  modified_power <- lm(log(Y) ~ log(x + C), data = data)

  # Extract parameters and R-squared
  get_params <- function(model) {
    s <- summary(model)
    a <- round(coef(model)[1], 3)
    b <- round(coef(model)[2], 3)
    R2 <- round(s$r.squared, 3)
    a_back_transformed <- round(exp(a), 3)

    se_a <- round(s$coefficients["(Intercept)", "Std. Error"], 3)
    predictor_name <- names(coef(model))[2]
    se_b <- round(s$coefficients[predictor_name, "Std. Error"], 3)

    sig_symbol <- function(p_val) {
      if (p_val < 0.01) return("**")
      else if (p_val < 0.05) return("*")
      else return("")
    }

    sig_a <- sig_symbol(s$coefficients["(Intercept)", "Pr(>|t|)"])
    sig_b <- sig_symbol(s$coefficients[predictor_name, "Pr(>|t|)"])

    return(list(
      a = a, se_a = se_a, sig_a = sig_a,
      b = b, se_b = se_b, sig_b = sig_b,
      a_back = a_back_transformed, R2 = R2
    ))
  }

  exp_params <- get_params(exponential)
  power_params <- get_params(power)
  mod_power_params <- get_params(modified_power)

  # Combine results into a table
  results <- rbind(
    Exponential = unlist(exp_params),
    Power = unlist(power_params),
    Modified_Power = unlist(mod_power_params)
  )

  # Plotting
  plot_exponential <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = log(Y))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black") +
    ggplot2::ggtitle("Exponential Model")

  plot_power <- ggplot2::ggplot(data, ggplot2::aes(x = log(x), y = log(Y))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black") +
    ggplot2::ggtitle("Power Model")

  data$mod_x <- log(data$x + C)
  plot_modified_power <- ggplot2::ggplot(data, ggplot2::aes(x = mod_x, y = log(Y))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black") +
    ggplot2::ggtitle("Modified Power Model")

  plot_exponential_original <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = Y)) +
    ggplot2::geom_point() +
    ggplot2::stat_function(fun = function(x) exp_params$a_back * exp(exp_params$b * x), color = "black") +
    ggplot2::ggtitle("Exponential Model Fit")

  plot_power_original <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = Y)) +
    ggplot2::geom_point() +
    ggplot2::stat_function(fun = function(x) power_params$a_back * (x^power_params$b), color = "black") +
    ggplot2::ggtitle("Power Model Fit")

  plot_modified_power_original <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = Y)) +
    ggplot2::geom_point() +
    ggplot2::stat_function(fun = function(x) mod_power_params$a_back * ((x + C)^mod_power_params$b), color = "black") +
    ggplot2::ggtitle("Modified Power Model Fit")

  return(list(
    data = data,
    results_table = results,
    plot_exponential = plot_exponential,
    plot_power = plot_power,
    plot_modified_power = plot_modified_power,
    plot_exponential_original = plot_exponential_original,
    plot_power_original = plot_power_original,
    plot_modified_power_original = plot_modified_power_original
  ))
}
