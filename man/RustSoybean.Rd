\name{RustSoybean}
\alias{RustSoybean}
\docType{data}
\title{
Soybean rust severity in Brazil
}
\description{
Soybean rust severity
}
\usage{data("RustSoybean")}
\format{
  A data frame with 34 observations on the following 7 variables.
  \describe{
    \item{\code{epidemia}}{a numeric vector}
    \item{\code{latitude}}{a numeric vector}
    \item{\code{longitude}}{a numeric vector}
    \item{\code{local}}{a character vector}
    \item{\code{planting}}{a Date}
    \item{\code{detection}}{a Date}
    \item{\code{severity}}{a numeric vector}
  }
}
\details{
Text
}
\references{
Del Ponte et al. (2006)
}
\examples{
data(RustSoybean)
## maybe str(RustSoybean) ; plot(RustSoybean) ...
}
\keyword{datasets}
