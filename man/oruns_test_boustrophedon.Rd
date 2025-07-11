% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oruns_test_boustrophedon.R
\name{oruns_test_boustrophedon}
\alias{oruns_test_boustrophedon}
\title{Boustrophedon Run Test for Binary Matrix}
\usage{
oruns_test_boustrophedon(mat)
}
\arguments{
\item{mat}{A binary matrix (containing 0s and 1s, and possibly NAs).}
}
\value{
A list with two elements:
\item{rowwise_boustrophedon}{List containing the sequence and result of \code{oruns_test} for row-wise traversal.}
\item{colwise_boustrophedon}{List containing the sequence and result of \code{oruns_test} for column-wise traversal.}
}
\description{
Applies the ordinary runs test to a binary matrix using boustrophedon-style traversal.
The function supports two modes: row-wise and column-wise boustrophedon. Each traversal flattens the matrix into a 1D sequence
which is then tested using \code{oruns_test}.
}
\seealso{
\code{\link{oruns_test}}

Other Spatial analysis: 
\code{\link{AFSD}()},
\code{\link{BPL}()},
\code{\link{count_subareas}()},
\code{\link{count_subareas_random}()},
\code{\link{fit_gradients}()},
\code{\link{join_count}()},
\code{\link{oruns_test}()},
\code{\link{oruns_test_byrowcol}()},
\code{\link{plot_AFSD}()}
}
\concept{Spatial analysis}
