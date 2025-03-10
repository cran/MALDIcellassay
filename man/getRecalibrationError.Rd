% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters.R
\name{getRecalibrationError}
\alias{getRecalibrationError}
\title{Calculate remaining calibration error of a MALDIassay object}
\usage{
getRecalibrationError(object)
}
\arguments{
\item{object}{Object of class MALDIassay}
}
\value{
A tibble containing statistics about remaining calibration error
}
\description{
Calculate remaining calibration error of a MALDIassay object
}
\examples{
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getRecalibrationError(Blank2022res)
}
