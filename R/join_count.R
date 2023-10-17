
#' Test for Spatial Join Count Statistics
#'
#' @description
#' The function `join_count` calculates spatial statistics for a matrix.
#' It identifies patterns of aggregation for values in a binary matrix based on
#' join count statistics. The results determine whether the observed spatial
#' arrangement is aggregated, random, or uniform regarding the presence of 1s
#' ("D") and 0s ("H").
#'
#' @param matrix_data A binary matrix (with elements 0 and 1) representing the
#'   spatial distribution of two types of points. Rows and columns of the matrix
#'   indicate coordinates while the 0s and 1s represent different categories of
#'   observations.
#'
#' @return A list containing several statistical measures:
#'   \itemize{
#'     \item Observed_01_10: The number of observed 01 or 10 sequences.
#'     \item Observed_11: The number of observed 11 sequences.
#'     \item Expected_ER_01: The expected number of 01 or 10 sequences under
#'           random distribution.
#'     \item Expected_ER_11: The expected number of 11 sequences under random
#'           distribution.
#'     \item sR_01: Standard deviation for 01 or 10 sequences.
#'     \item sR_11: Standard deviation for 11 sequences.
#'     \item ZHD: The Z-score for 01 or 10 sequences.
#'     \item ZDD: The Z-score for 11 sequences.
#'     \item PatternHD: Descriptive result for 01 or 10 sequences - "Aggregated"
#'           or "Not Aggregated".
#'     \item PatternDD: Descriptive result for 11 sequences - "Aggregated" or
#'           "Not Aggregated".
#'   }
#'
#' @details
#' The function first calculates the count of specific binary sequences in the
#' matrix both horizontally and vertically. It then computes expected values
#' and standard deviations based on the spatial arrangement dimensions and uses
#' these for Z-score calculations to identify the pattern of the matrix.
#'
#' @examples
#' matrix_data <- matrix(c(1,1,1,0,0,
#'                         1,1,1,0,0,
#'                         1,1,1,0,0,
#'                         1,1,1,0,0,
#'                         0,0,0,0,0), ncol = 5, byrow = TRUE)
#' join_count(matrix_data)
#'
#' @export
join_count <- function(matrix_data) {

  count_sequences <- function(matrix_data) {
    count_01_10 <- 0  # previously count_HD_DH
    count_11 <- 0     # previously count_DD

    for (i in 1:nrow(matrix_data)) {
      for (j in 1:(ncol(matrix_data)-1)) {
        if ((matrix_data[i,j] == 0 && matrix_data[i,j+1] == 1) ||
            (matrix_data[i,j] == 1 && matrix_data[i,j+1] == 0)) {
          count_01_10 <- count_01_10 + 1
        } else if (matrix_data[i,j] == 1 && matrix_data[i,j+1] == 1) {
          count_11 <- count_11 + 1
        }
      }
    }

    for (j in 1:ncol(matrix_data)) {
      for (i in 1:(nrow(matrix_data)-1)) {
        if ((matrix_data[i,j] == 0 && matrix_data[i+1,j] == 1) ||
            (matrix_data[i,j] == 1 && matrix_data[i+1,j] == 0)) {
          count_01_10 <- count_01_10 + 1
        } else if (matrix_data[i,j] == 1 && matrix_data[i+1,j] == 1) {
          count_11 <- count_11 + 1
        }
      }
    }

    return(list(One_Zero_or_Zero_One = count_01_10, One_One = count_11))
  }

  results <- count_sequences(matrix_data)

  # Calculate matrix dimensions and S-values
  num_rows <- nrow(matrix_data)
  num_cols <- ncol(matrix_data)

  S0 <- 2 * ((2*num_rows * num_cols) - num_rows - num_cols)
  S1 <- 2 * S0
  S2 <- 8 * ((8*num_rows * num_cols) - (7 * num_rows) - (7 * num_cols) + 4)

  # Calculate y as the mean of the numeric matrix
  y <- mean(matrix_data)

  # Calculate expected values and standard deviations
  ER_01 <- S0 * y * (1 - y)  # Expected random 01 or 10 sequences
  sR_01 <- sqrt(S2 * y * (1 - y) + 4 * (S1 - S2) * y^2 * (1 - y)^2) / 2  # Std dev

  ER_11 <- S0 * y^2 / 2  # Expected random 11 sequences
  sR_11 <- sqrt(S1 * y^2 + (S2 - 2 * S1) * y^3 + (S1 - S2) * y^4) / 2  # Std dev

  # Calculate Z value
  ZHD <- (results$One_Zero_or_Zero_One - ER_01 - 0.5) / sR_01

  pattern1 <- ifelse(ZHD < -1.64, "Aggregated", "Not Aggregated")

  ZDD <- (results$One_One - ER_11 + 0.5) / sR_11

  pattern2 <- ifelse(ZDD > 1.64, "Aggregated", "Not Aggregated")

  # Return results
  return(list(Observed_01_10 = results$One_Zero_or_Zero_One, Observed_11 = results$One_One,
              Expected_ER_01 = ER_01, Expected_ER_11 = ER_11,
              sR_01 = sR_01, sR_11 = sR_11,
              ZHD = ZHD, ZDD = ZDD, PatternHD = pattern1, PatternDD = pattern2))
}



