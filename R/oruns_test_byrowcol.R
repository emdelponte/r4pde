#' @title Runs Test for Each Row and Column of a Binary Matrix
#' @description Applies the ordinary runs test to each row and column of a binary matrix individually.
#'
#' @param mat A binary matrix (containing 0s and 1s, and possibly NAs).
#' @return A list with four elements:
#' 	\item{row_results}{Data frame with test results for each row.}
#' 	\item{col_results}{Data frame with test results for each column.}
#' 	\item{row_summary}{Percentage summary of interpretation for rows.}
#' 	\item{col_summary}{Percentage summary of interpretation for columns.}
#' @export
#' @seealso \code{\link{oruns_test}}
#' @importFrom stats na.omit
#' @family Spatial analysis

oruns_test_byrowcol <- function(mat) {
  if (!is.matrix(mat) || !all(mat %in% c(0, 1, NA)))
    stop("Input must be a binary matrix (0/1), possibly with NAs.")

  apply_oruns <- function(x) {
    x <- na.omit(x)
    if (length(x) < 2 || all(x == x[1])) return(NULL)
    oruns_test(x)
  }

  build_df <- function(res) {
    if (is.null(res)) {
      return(data.frame(U = NA, EU = NA, sU = NA, Z = NA, pvalue = NA, result = NA))
    } else {
      return(data.frame(U = res$U, EU = res$EU, sU = res$sU, Z = res$Z, pvalue = res$pvalue, result = res$result))
    }
  }

  row_tests <- lapply(seq_len(nrow(mat)), function(i) apply_oruns(mat[i, ]))
  col_tests <- lapply(seq_len(ncol(mat)), function(j) apply_oruns(mat[, j]))

  row_df <- do.call(rbind, lapply(row_tests, build_df))
  col_df <- do.call(rbind, lapply(col_tests, build_df))

  row_summary <- round(prop.table(table(na.omit(row_df$result))) * 100, 2)
  col_summary <- round(prop.table(table(na.omit(col_df$result))) * 100, 2)

  list(
    row_results = row_df,
    col_results = col_df,
    row_summary = as.data.frame(row_summary),
    col_summary = as.data.frame(col_summary)
  )
}
