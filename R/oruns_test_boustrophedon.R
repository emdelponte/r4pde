#' @title Boustrophedon Run Test for Binary Matrix
#' @description Applies the ordinary runs test to a binary matrix using boustrophedon-style traversal.
#' The function supports two modes: row-wise and column-wise boustrophedon. Each traversal flattens the matrix into a 1D sequence
#' which is then tested using `oruns_test`.
#'
#' @param mat A binary matrix (containing 0s and 1s, and possibly NAs).
#' @return A list with two elements:
#' 	\item{rowwise_boustrophedon}{List containing the sequence and result of `oruns_test` for row-wise traversal.}
#' 	\item{colwise_boustrophedon}{List containing the sequence and result of `oruns_test` for column-wise traversal.}
#' @export
#' @seealso \code{\link{oruns_test}}
#' @importFrom stats na.omit

oruns_test_boustrophedon <- function(mat) {
  if (!is.matrix(mat) || !all(mat %in% c(0, 1, NA))) {
    stop("Input must be a binary matrix (0/1), possibly with NAs.")
  }

  # Row-wise boustrophedon
  row_seq <- unlist(lapply(seq_len(nrow(mat)), function(i) {
    row_i <- mat[i, ]
    if (i %% 2 == 0) rev(row_i) else row_i
  }))
  row_seq <- na.omit(row_seq)

  # Column-wise boustrophedon
  col_seq <- unlist(lapply(seq_len(ncol(mat)), function(j) {
    col_j <- mat[, j]
    if (j %% 2 == 0) rev(col_j) else col_j
  }))
  col_seq <- na.omit(col_seq)

  row_run <- oruns_test(row_seq)
  col_run <- oruns_test(col_seq)

  list(
    rowwise_boustrophedon = list(
      sequence = row_seq,
      test = row_run
    ),
    colwise_boustrophedon = list(
      sequence = col_seq,
      test = col_run
    )
  )
}
