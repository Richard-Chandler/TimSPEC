\name{dmvnorm.SpD}
\alias{dmvnorm.SpD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density of multivariate normal distribution}
\description{To calculate the density of a multivariate normal distribution parameterised in terms of its mean vector and the spectral decomposition (eigendecomposition) of its covariance matrix.}
\usage{
dmvnorm.SpD(x, mu, Cov.SpD, logged = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \code{n*d} matrix, where each row represents a point in \eqn{R^d} for which the density is required.
%%     ~~Describe \code{x} here~~
}
  \item{mu}{Vector of length \code{d}: the mean vector of the distribution}
  \item{Cov.SpD}{A list, with components \code{values} and \code{vectors}, defining the spectral decomposition of the covariance matrix of the distribution. Typically the result of a call to \code{\link{eigen}}. The eigenvalues must all be strictly positive.}
  \item{logged}{Return the density (\code{FALSE}, the default) or its logarithm?}
}
\value{a vector of length \code{n}, containing the (log) densities.}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\note{Output from this routine has been cross-checked against that from 
\code{dmvnorm} in the \pkg{mvtnorm} library.}
