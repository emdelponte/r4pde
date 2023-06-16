#' Calculate the Disease severity Index (DSI) (frequency of each class)
#'
#' This function calculates the Disease Severity Index (DSI) given a vector of classes,
#' a vector of frequencies, and a maximum possible class value. The DSI is calculated
#' as a weighted sum of class values, where each class is multiplied by its corresponding
#' frequency, then divided by the product of the total frequency and maximum class value,
#' and finally multiplied by 100 to get a percentage.
#'
#' @param class A numeric vector representing the classes.
#' @param freq A numeric vector representing the frequency of each class.
#'     Must be the same length as 'class'.
#' @param max A numeric value representing the maximum possible class value.
#'
#' @return Returns a single numeric value representing the DSI.
#'
#' @examples
#' DSI2(c(0, 1, 2, 3, 4), c(2, 0, 5, 0, 5), 4)
#'
#' @export
DSI2 <- function(class, freq, max) {
  df <- data.frame(class, freq)

  # Calculate weighted sum of class * frequency
  weighted_sum <- sum(df$class * df$freq)

  # Calculate total frequency
  total_freq <- sum(df$freq)

  # Calculate DSI
  DSI_value <- (weighted_sum / (total_freq * max)) * 100

  return(DSI_value)
}
