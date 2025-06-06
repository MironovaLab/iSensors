% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iSensor.R
\name{add_iSensor_labels}
\alias{add_iSensor_labels}
\title{Add Sample Labels to iSensor Object}
\usage{
add_iSensor_labels(iSensor_obj, labels)
}
\arguments{
\item{iSensor_obj}{An iSensor object created by \code{create_iSensor}}

\item{labels}{Character vector of sample labels (length must match number of samples)}
}
\value{
The modified iSensor object with updated sample labels
}
\description{
Adds or updates sample labels in an iSensor object. The labels must match
the number of samples (columns) in the expression data matrix.
}
\examples{
# Create iSensor object
iSensor_obj <- create_iSensor(data = matrix(rnorm(100), nrow = 10))

# Add sample labels
sample_labels <- paste0("Sample", 1:10)
iSensor_obj <- add_iSensor_labels(iSensor_obj, sample_labels)

}
