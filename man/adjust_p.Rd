% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics_functions.R
\name{adjust_p}
\alias{adjust_p}
\title{Adjust p-values stored in a vector or matrix using the Holm-Bonferroni correction}
\usage{
adjust_p(mat_with_p, method = "holm")
}
\arguments{
\item{mat_with_p}{A double. Vector or matrix containing a column of p-values to be
adjusted.}

\item{method}{A string defining the name of the p-value adjustment method to
be used. Defautls to Holm's method.}
}
\value{
A matrix identical to the original, with an added column
  containing the adjusted p-values.
}
\description{
Adjust p-values stored in a vector or matrxi using the Holm-Bonferroni
correction. The cell name containing the p.values needs to be named
'p.value'.
}
