% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters.R
\name{getAllMz}
\alias{getAllMz}
\title{Get all mz value of an MALDIassay-object}
\usage{
getAllMz(object, excludeNormMz = FALSE)
}
\arguments{
\item{object}{Object of class MALDIassay}

\item{excludeNormMz}{Logical, remove normMz from list of mz values.}
}
\value{
numeric vector of mz values
}
\description{
Get all mz value of an MALDIassay-object
}
\examples{
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getAllMz(Blank2022res))
}
