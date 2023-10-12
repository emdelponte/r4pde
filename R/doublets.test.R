#' Doublets Test
#'
#' Perform a doublets test on the input data to test for aggregation or clustering or randomness.
#'
#' @param x A numeric vector representing the input data
#' @return A list with elements Db (number of doublets), EDb (expected number of doublets), pvalue, and result (either "aggregation or clustering" or "randomness").
#' @examples
#' doublets.test(c(1, 0, 1, 1, 0, 1, 0, 0, 1, 1))
#'
#' @importFrom stats pnorm
#' @export
doublets.test <- function(x) {
  # S is the input data
  S <- x

  # Create a matrix by combining each element in S with the element that follows it
  matrix <- cbind(S[-length(S)], S[-1])

  # Count the number of times each pair appears in the matrix
  pairs <- table(data.frame(matrix))

  # Db is the number of pairs that appear twice
  Db <- pairs[2, 2]

  # N is the length of the input data
  N <- length(S)

  # m is the sum of the input data
  m = sum(S)

  # EDb is the expected number of doublets
  EDb = m * ((m - 1) / N)

  # SDb is the standard deviation of the expected number of doublets
  SDb = sqrt(EDb * (1 - (2 / N)))

  # ZDb is the Z-score of the observed number of doublets
  ZDb = (Db - EDb) / SDb

  # pvalue is the p-value of the Z-score
  pvalue <- (pnorm(abs(ZDb), lower.tail = FALSE))

  # result is "aggregation or clustering" if the Z-score is greater than 1.64, otherwise it is "randomness"
  result <- ifelse(abs(ZDb) >= 1.64,
                   c("aggregation or clustering"),
                   c("randomness"))

  # Return the results as a list
  return(list(Db = Db, EDb = EDb, ZDb = ZDb, pvalue = pvalue, result = result))
}
