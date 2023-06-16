##' Calculate the Disease Severity Index (DSI) (class for each unit)
#'
#' This function calculates the Disease Severity Index (DSI) based on the provided unit,
#' class, and maximum class value. The DSI is computed by aggregating the classes,
#' calculating weights by multiplying the frequency of each class by the class itself,
#' and then dividing the sum of these weights by the product of the total number of
#' entries and the maximum class value, then multiplying by 100.
#'
#' @param unit A vector representing the units.
#' @param class A vector representing the classes corresponding to the units.
#' @param max A numeric value representing the maximum possible class value.
#'
#' @return Returns a single numeric value representing the DSI.
#'
#' @examples
#' # Example usage:
#' unit <- c(1, 2, 3, 4, 5, 6)
#' class <- c(1, 2, 1, 2, 3, 1)
#' max <- 3
#' DSI(unit, class, max)
#'
#' @export
DSI <- function(unit, class, max) {
  df <- data.frame(unit, class)
  tab <- aggregate(df$class, list(num=class), length)
  tab$weight = tab$num * tab$x
  (sum(as.numeric(tab$weight)))/(nrow(df)*max)*100
}



