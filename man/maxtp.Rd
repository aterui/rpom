% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{maxtp}
\alias{maxtp}
\title{Utility: get maximum trophic position}
\usage{
maxtp(foodweb, occupancy, weight = TRUE)
}
\arguments{
\item{foodweb}{Binary food web matrix.}

\item{occupancy}{Equilibrium occupancies of constituent species}

\item{weight}{Logical.
If \code{TRUE}, maximum trophic position is weighted by relative occupancies.}
}
\description{
Utility: get maximum trophic position
}
\author{
Akira Terui, \email{hanabi0111@gmail.com}
}
