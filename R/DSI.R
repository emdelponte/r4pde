#' DSI
#' @description
#' Calculates the Disease Severity Index
#' @param unit A numeric vector representing the unit of measure
#' @param class A factor vector representing the class labels
#' @param max A numeric scalar representing the maximum possible value of the unit of measure
#' @return A numeric scalar representing the DSI
#' @importFrom stats aggregate
#' @export
#' @examples
#' unit <- c(1:12)
#' class <- c(2,3,1,1,3,4,5,0,2,5,2,1)
#' ratings <- data.frame(unit, class)
#' DSI(unit = ratings$unit, class = ratings$class, max = 6)
DSI <- function(unit, class, max) {
  df <- data.frame(unit, class)
  tab <- aggregate(df$class, list(num=class), length)
  tab$weight = tab$num * tab$x
  (sum(as.numeric(tab$weight)))/(nrow(df)*max)*100
}


