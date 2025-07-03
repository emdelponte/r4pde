#' Test for Spatial Join Count Statistics
#'
#' @description
#' The function `join_count` calculates spatial join count statistics for a binary matrix,
#' identifying patterns of aggregation or randomness.
#'
#' @param matrix_data A binary matrix (with elements 0 and 1) representing the spatial distribution
#'   of two types of points: 0 for healthy plants (H) and 1 for diseased plants (D). This matrix reflects
#'   the geographical distribution or layout of plants in the studied area.
#' @param verbose Logical. If TRUE (default), prints a formatted message to the console.
#' @return A comprehensive, rich-text formatted string of results that includes:
#'   \itemize{
#'     \item Statistical counts of specific binary sequences (e.g., "01 or 10", "11")
#'     \item Expected counts under the assumption of Complete Spatial Randomness (CSR)
#'     \item Standard deviations and Z-scores (ZHD for "01 or 10" sequences, ZDD for "11" sequences)
#'     \item Interpretation of whether the spatial distribution for each sequence type is
#'           "Aggregated" or "Not Aggregated" based on Z-scores
#'     \item A summary explaining the implications of these statistics and patterns
#'   }
#'   The return value aims to provide a clear understanding of the spatial arrangement's
#'   characteristics, aiding in further spatial analysis or research.
#'
#' @details
#' The function conducts an analysis by first counting the occurrence of specific sequences
#' ("01 or 10" and "11" - equivalent to HD and DD) in the binary matrix.
#'  It then calculates expected values, standard
#' deviations, and Z-scores to determine the spatial randomness or aggregation. The analysis
#' considers both horizontal and vertical adjacency (rook case) in the matrix.
#' @references
#' Madden, L. V., Hughes, G., & van den Bosch, F. (2007). The Study of Plant Disease Epidemics.
#'   The American Phytopathological Society.
#'
#' @export
#' @family Spatial analysis
join_count <- function(matrix_data, verbose = TRUE) {

  count_sequences <- function(matrix_data) {
    count_01_10 <- 0
    count_11 <- 0

    for (i in 1:nrow(matrix_data)) {
      for (j in 1:(ncol(matrix_data) - 1)) {
        if ((matrix_data[i, j] == 0 && matrix_data[i, j + 1] == 1) ||
            (matrix_data[i, j] == 1 && matrix_data[i, j + 1] == 0)) {
          count_01_10 <- count_01_10 + 1
        } else if (matrix_data[i, j] == 1 && matrix_data[i, j + 1] == 1) {
          count_11 <- count_11 + 1
        }
      }
    }

    for (j in 1:ncol(matrix_data)) {
      for (i in 1:(nrow(matrix_data) - 1)) {
        if ((matrix_data[i, j] == 0 && matrix_data[i + 1, j] == 1) ||
            (matrix_data[i, j] == 1 && matrix_data[i + 1, j] == 0)) {
          count_01_10 <- count_01_10 + 1
        } else if (matrix_data[i, j] == 1 && matrix_data[i + 1, j] == 1) {
          count_11 <- count_11 + 1
        }
      }
    }

    list(One_Zero_or_Zero_One = count_01_10, One_One = count_11)
  }

  results <- count_sequences(matrix_data)
  num_rows <- nrow(matrix_data)
  num_cols <- ncol(matrix_data)
  S0 <- 2 * ((2 * num_rows * num_cols) - num_rows - num_cols)
  S1 <- 2 * S0
  S2 <- 8 * ((8 * num_rows * num_cols) - (7 * num_rows) - (7 * num_cols) + 4)
  y <- mean(matrix_data)

  ER_01 <- S0 * y * (1 - y)
  sR_01 <- sqrt(S2 * y * (1 - y) + 4 * (S1 - S2) * y^2 * (1 - y)^2) / 2
  ER_11 <- S0 * y^2 / 2
  sR_11 <- sqrt(S1 * y^2 + (S2 - 2 * S1) * y^3 + (S1 - S2) * y^4) / 2

  ZHD <- (results$One_Zero_or_Zero_One - ER_01 - 0.5) / sR_01
  ZDD <- (results$One_One - ER_11 + 0.5) / sR_11

  pattern_HD <- ifelse(ZHD < -1.64, "aggregated", "not aggregated")
  pattern_DD <- ifelse(ZDD > 1.64, "aggregated", "not aggregated")

  output <- list(
    observed = list(HD = results$One_Zero_or_Zero_One, DD = results$One_One),
    expected = list(HD = ER_01, DD = ER_11),
    sd = list(HD = sR_01, DD = sR_11),
    zscore = list(HD = ZHD, DD = ZDD),
    pattern = list(HD = pattern_HD, DD = pattern_DD)
  )

  if (verbose) {
    message(
      "Join Count Analysis of Spatial Patterns of Plant Diseases:\n",
      "----------------------------------------------------------\n",
      "1) 'HD' Sequences:\n",
      sprintf("   - Observed Count: %d\n", output$observed$HD),
      sprintf("   - Expected Count: %.2f\n", output$expected$HD),
      sprintf("   - Standard Deviation: %.2f\n", output$sd$HD),
      sprintf("   - Z-score: %.2f\n", output$zscore$HD),
      sprintf("   - Pattern: %s\n\n", output$pattern$HD),
      "2) 'DD' Sequences:\n",
      sprintf("   - Observed Count: %d\n", output$observed$DD),
      sprintf("   - Expected Count: %.2f\n", output$expected$DD),
      sprintf("   - Standard Deviation: %.2f\n", output$sd$DD),
      sprintf("   - Z-score: %.2f\n", output$zscore$DD),
      sprintf("   - Pattern: %s\n", output$pattern$DD),
      "----------------------------------------------------------\n"
    )
  }

  invisible(output)
}
