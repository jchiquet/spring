% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{diffOp}
\alias{diffOp}
\title{Taken from genlasso package. Usefull to construct the structring
matrix}
\usage{
diffOp(p, K = 1)
}
\arguments{
\item{p}{an integer for the number of features}

\item{K}{an integer for the order (default is 1) for the order of
the derivative operator}
}
\value{
an sparse matrix with class \code{dgCMatrix}.
}
\description{
Taken from genlasso package. Usefull to construct the structring
matrix
}
