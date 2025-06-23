#' Count the Number of Ones in Subareas of a Matrix
#'
#' This function takes a binary matrix (0s and 1s) and divides it into rectangular subareas,
#' counting the number of ones in each. Subareas are defined by the number of rows and columns
#' specified by the user. If the matrix dimensions are not perfectly divisible by the subarea size,
#' edge subareas may be smaller.
#'
#' @param matrix_data A matrix of 0s and 1s to analyze.
#' @param sub_rows Number of rows in each subarea.
#' @param sub_cols Number of columns in each subarea.
#'
#' @return A matrix where each cell corresponds to a subarea and contains the count of ones.
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(sample(c(0, 1), 12 * 16, replace = TRUE), nrow = 16, ncol = 12)
#' count_matrix <- count_subareas(mat, sub_rows = 3, sub_cols = 3)
#' print(count_matrix)
#'
#' @export
#' @family Spatial analysis


# Function to count 1s in subareas
count_subareas <- function(matrix_data, sub_rows, sub_cols) {
  # Get matrix dimensions
  total_rows <- nrow(matrix_data)
  total_cols <- ncol(matrix_data)

  # Calculate number of subareas in each dimension
  n_sub_rows <- ceiling(total_rows / sub_rows)
  n_sub_cols <- ceiling(total_cols / sub_cols)

  # Initialize result matrix
  result <- matrix(NA, nrow = n_sub_rows, ncol = n_sub_cols)

  # Iterate through each subarea
  for (i in 1:n_sub_rows) {
    for (j in 1:n_sub_cols) {
      # Calculate row and column indices for current subarea
      row_start <- (i-1)*sub_rows + 1
      row_end <- min(i*sub_rows, total_rows)
      col_start <- (j-1)*sub_cols + 1
      col_end <- min(j*sub_cols, total_cols)

      # Extract submatrix
      submatrix <- matrix_data[row_start:row_end, col_start:col_end]

      # Count 1s in submatrix
      result[i,j] <- sum(submatrix == 1, na.rm = TRUE)
    }
  }

  return(result)
}

