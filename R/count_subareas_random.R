#' @title Random Subgrid Sampling of a Binary Matrix
#' @description Randomly samples submatrices (quadrats) of specified size from a binary matrix,
#' and returns the positions, submatrices, and count of 1s in each sampled quadrat.
#'
#' @param matrix_data A binary matrix of 0s and 1s.
#' @param sub_rows Number of rows in each subgrid sample.
#' @param sub_cols Number of columns in each subgrid sample.
#' @param n_samples Number of subgrid samples to draw.
#'
#' @return A list of sampled subgrids. Each element is a list with:
#' \item{position}{Row and column start position of the sample.}
#' \item{submatrix}{The sampled subgrid matrix.}
#' \item{count}{Number of 1s in the sampled submatrix.}
#' @export
#' @family Spatial analysis

count_subareas_random <- function(matrix_data, sub_rows = 3, sub_cols = 3, n_samples = 100) {
  total_rows <- nrow(matrix_data)
  total_cols <- ncol(matrix_data)
  max_row_start <- total_rows - sub_rows + 1
  max_col_start <- total_cols - sub_cols + 1

  samples <- lapply(1:n_samples, function(i) {
    row_start <- sample(1:max_row_start, 1)
    col_start <- sample(1:max_col_start, 1)

    submatrix <- matrix_data[
      row_start:(row_start + sub_rows - 1),
      col_start:(col_start + sub_cols - 1)
    ]

    list(
      position = c(row_start, col_start),
      submatrix = submatrix,
      count = sum(submatrix == 1, na.rm = TRUE)
    )
  })

  list(
    samples = samples,
    counts = sapply(samples, function(x) x$count)
  )
}
