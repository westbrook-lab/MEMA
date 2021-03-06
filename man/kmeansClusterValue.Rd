% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derivedParameter_functions.R
\name{kmeansClusterValue}
\alias{kmeansClusterValue}
\title{Cluster using kmeans}
\usage{
kmeansClusterValue(x, centers = 2)
}
\arguments{
\item{x}{A numeric vector to be clustered}

\item{centers}{The number of centers or clusters to find.}
}
\value{
The cluster assignments for x using the base kmeans command.
}
\description{
\code{kmeansClusterValue} is a wrapper function for perfoming kmeans clustering
}

