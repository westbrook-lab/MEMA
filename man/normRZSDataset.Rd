% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_functions.R
\name{normRZSDataset}
\alias{normRZSDataset}
\title{Normalize selected values in a dataset on a plate basis}
\usage{
normRZSDataset(dt)
}
\arguments{
\item{dt}{A data.table with a \code{Barcode} column numeric values to be RZS normalized 
using all ECM proteins in the FBS well}
}
\value{
A datatable with the normalized values
}
\description{
A wrapper function for \code{normRZSWellsWithinPlate} that selects the
\code{_CP_|_QI_|_PA_|SpotCellCount|Lineage} columns of dt if they exist and 
normalizes them on a plate basis
}

