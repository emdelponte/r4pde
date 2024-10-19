\name{BlastWheat}
\alias{BlastWheat}
\docType{data}
\title{
Wheat blast data in 143 field trials in Brazil
}
\description{
Data on wheat blast incidence and severity and yield collected in 143 wheat field trials in Brazil
}
\usage{data("BlastWheat")}
\format{
  A data frame with 143 observations on the following 10 variables.
  \describe{
    \item{\code{study}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{location}}{a character vector}
    \item{\code{state}}{a character vector}
    \item{\code{latitude}}{a numeric vector}
    \item{\code{longitude}}{a numeric vector}
    \item{\code{heading}}{a character vector}
    \item{\code{inc_mean}}{a numeric vector}
    \item{\code{index_mean}}{a numeric vector}
    \item{\code{yld_mean}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
De Cól, M., Coelho, M., & Del Ponte, E. M. (2024). Weather-based logistic regression models for predicting wheat head blast epidemics. Plant Disease, 108(7), 2206-2213.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(BlastWheat)
## maybe str(BlastWheat) ; plot(BlastWheat) ...
}
\keyword{datasets}
