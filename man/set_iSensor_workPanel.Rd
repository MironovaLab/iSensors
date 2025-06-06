% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iSensor.R
\name{set_iSensor_workPanel}
\alias{set_iSensor_workPanel}
\title{Set Active Gene Panels for iSensor Analysis}
\usage{
set_iSensor_workPanel(iSensor_obj, panelSet = "allSpec_allHorm_allType")
}
\arguments{
\item{iSensor_obj}{An iSensor analysis object created by \code{create_iSensor}}

\item{panelSet}{Name of the panel set to activate (default: 'allSpec_allHorm_allType')}
}
\value{
The modified iSensor object with updated active gene panels
}
\description{
Selects a predefined set of gene panels to use for subsequent analysis steps
in the iSensor pipeline. The specified panel set must exist in the iSensor object.
}
\examples{
\dontrun{
# After creating an iSensor object:
iSensor_obj <- set_iSensor_workPanel(iSensor_obj, panelSet = "my_panels")
}

}
