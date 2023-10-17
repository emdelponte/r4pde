#' Test for Spatial Join Count Statistics
#'
#' @description
#' The function `join_count` calculates spatial join count statistics for a binary matrix,
#' identifying patterns of aggregation or randomness.
#'
#' @param matrix_data A binary matrix (with elements 0 and 1) representing the spatial distribution
#'   of two types of points: 0 for healthy plants (H) and 1 for diseased plants (D). This matrix reflects
#'   the geographical distribution or layout of plants in the studied area.
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
#'
#' @examples
#' \dontrun{
#' matrix_data <- matrix(c(1,1,1,0,0,
#'                         1,1,1,0,0,
#'                         1,1,1,0,0,
#'                         1,1,1,0,0,
#'                         0,0,0,0,0),
#'                         ncol = 5, byrow = TRUE)
#' result_text <- join_count(matrix_data)
#' cat(result_text)  # Print the rich text output
#' }
#'
#' @references
#' Madden, L. V., Hughes, G., & van den Bosch, F. (2007). The Study of Plant Disease Epidemics.
#'   The American Phytopathological Society.
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

  # Construct the result sentence based on Z-scores.
  spatial_pattern_HD <- ifelse(ZHD < -1.64, "aggregated", "not aggregated")
  spatial_pattern_DD <- ifelse(ZDD > 1.64, "aggregated", "not aggregated")

  # Prepare a rich, informative text output.
  message <- paste(
    "Join Count Analysis of Spatial Patterns of Plant Diseases:\n",
    "----------------------------------------------------------\n",
    "This analysis is based on the assessment of spatial patterns within the provided matrix data, focusing on the occurrences of 'HD' (pair of healthy and diseased plant) and 'DD' (pair of diseased plants) sequences. Specifications and formulations as shown in Madden et al. (2007).\n\n",
    "1) 'HD' Sequences:\n",
    sprintf("   - Observed Count: %d\n", results$One_Zero_or_Zero_One),
    sprintf("   - Expected Count : %.2f\n", ER_01),
    sprintf("   - Standard Deviation: %.2f\n", sR_01),
    sprintf("   - Z-score: %.2f\n", ZHD),
    sprintf("Based on the Z score, the pattern for 'HD' sequences is '%s'.\n\n", spatial_pattern_HD),
    "2) 'DD' Sequences:\n",
    sprintf("   - Observed Count: %d\n", results$One_One),
    sprintf("   - Expected Count: %.2f\n", ER_11),
    sprintf("   - Standard Deviation: %.2f\n", sR_11),
    sprintf("   - Z-score: %.2f\n", ZDD),
    sprintf("Based on the Z score, the pattern for 'DD' sequences is '%s'.\n", spatial_pattern_DD),
    "-----------------------------------------------------\n"

  )

  # Print the message for immediate visibility in the console.
  cat(message)

  # You can also return this message as a value, in addition to other statistical details, if needed.
  #return(message)
}






