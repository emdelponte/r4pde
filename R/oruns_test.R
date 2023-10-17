#' Runs Test
#'
#' Perform a runs test on the input data to test for clustering or randomness.
#'
#' @param x A numeric vector representing the input data
#' @return A list with elements U (number of runs), EU (expected number of runs), pvalue, and result (either "clustering" or "randomness").
#' @examples
#' oruns_test(c(1, 0, 1, 1, 0, 1, 0, 0, 1, 1))
#'
#' @importFrom stats pnorm
#' @export
oruns_test <- function(x) {
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
  pvalue <- (1 - pnorm(abs(Z), lower.tail = TRUE))

  result <- ifelse(Z < -1.64,
                   c("aggregation or clustering"),
                   c("randomness"))
  message <- paste(
    "Ordinary Runs Test of Data Sequence:\n",
    "-------------------------------------\n",
    sprintf("Total Number of Runs (U): %d\n", U),
    sprintf("Expected Number of Runs (EU): %.2f\n", EU),
    sprintf("Standard Deviation of Runs (sU): %.2f\n", sU),
    sprintf("Z-score: %.2f\n", Z),
    sprintf("P-value: %.4f\n\n", pvalue),
    "Interpretation:\n",
    sprintf("Based on the Z-score, the sequence exhibits '%s'.\n", result)

  )

  cat(message)
}
