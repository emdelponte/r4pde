\name{WhiteMoldSoybean}
\alias{WhiteMoldSoybean}
\docType{data}
\title{
White mold on soybean in Brazil
}
\description{
The raw data was organized in the long or tidy format where each treatment (observation) in a fungicide trial (hereafter study) is placed in its row and each variable in its column.
}
\usage{data("WhiteMoldSoybean")}
\format{
  A data frame with 382 observations on the following 17 variables.
  \describe{
    \item{\code{study}}{a numeric vector}
    \item{\code{treat}}{a numeric vector}
    \item{\code{season}}{a character vector}
    \item{\code{harvest_year}}{a numeric vector}
    \item{\code{location}}{a character vector}
    \item{\code{state}}{a character vector}
    \item{\code{country}}{a character vector}
    \item{\code{elevation}}{a numeric vector}
    \item{\code{region}}{a character vector}
    \item{\code{elevation_class}}{a character vector}
    \item{\code{inc_check}}{a numeric vector}
    \item{\code{inc_class}}{a character vector}
    \item{\code{yld_check}}{a numeric vector}
    \item{\code{yld_class}}{a character vector}
    \item{\code{inc}}{a numeric vector}
    \item{\code{scl}}{a numeric vector}
    \item{\code{yld}}{a numeric vector}
  }
}

\source{
https://github.com/emdelponte/paper-white-mold-meta-analysis
}
\references{
Lehner, M. S., Pethybridge, S. J., Meyer, M. C., and Del Ponte, E. M. 2016. Meta-analytic modelling of the incidenceyield and incidencesclerotial production relationships in soybean white mould epidemics. Plant Pathology. 66:460–468 Available at: http://dx.doi.org/10.1111/ppa.12590.
}
\examples{
data(WhiteMoldSoybean)
## maybe str(WhiteMoldSoybean) ; plot(WhiteMoldSoybean) ...
}
\keyword{datasets}
