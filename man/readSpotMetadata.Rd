% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gal_parse_functions.R
\name{readSpotMetadata}
\alias{readSpotMetadata}
\title{Read the spot metadata from a gal file}
\usage{
readSpotMetadata(galFile)
}
\arguments{
\item{galFile}{the name of a gal file to be read}
}
\value{
A dataframe with the contents of each spot in the gal file. The spots
 are converted from print block, row and column space to array row and column
 space. The array space is aligned with the print space. That is,
 the 1,1,1 position in the print space is the first row and first column of the
 array space. The Name and ID fields of spots that were not printed have their
 dashes replaced with the word blank.
}
\description{
\code{readSpotMetadata} returns a dataframe of the contents in a gal file and
adds array indices.
}

