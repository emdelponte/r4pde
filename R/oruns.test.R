#' Runs Test
#'
#' Perform a runs test on the input data to test for clustering or randomness.
#'
#' @param x A numeric vector representing the input data
#' @return A list with elements U (number of runs), EU (expected number of runs), pvalue, and result (either "clustering" or "randomness").
#' @examples
#' oruns.test(c(1, 0, 1, 1, 0, 1, 0, 0, 1, 1))
#'
#' @importFrom stats pnorm
#' @export
oruns.test <- function(x) {
  # S is the input data
  S <- x

  # U is the number of runs in the data
  U = max(cumsum(c(1, diff(S) != 0)))

  # m is the sum of the input data
  m = sum(S)

  # N is the length of the input data
  N = length(S)

  # EU is the expected number of runs
  EU = 1 + (2 * m * (N - m) / N)

  # sU is the standard deviation of the expected number of runs
  sU = sqrt(2 * m * (N - m) * (2 * m * (N - m) - N) / (N ^ 2 * (N - 1)))

  # Z is the Z-score of the observed number of runs
  Z = (U - EU) / sU

  # pvalue is the p-value of the Z-score
  pvalue <- (2 * pnorm(abs(Z), lower.tail = FALSE))

  # result is "clustering" if the Z-score is less than 1.64, otherwise it is "randomness"
  result <- ifelse(Z < 1.64,
                   c("clustering"),
                   c("randomness"))

  # Return the results as a list
  return(list(U = U, EU = EU, pvalue = pvalue, result = result))
}
