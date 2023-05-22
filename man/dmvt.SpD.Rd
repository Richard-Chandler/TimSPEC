\name{dmvt.SpD}
\alias{dmvt.SpD}
\title{Density of multivariate t distribution}
\description{To calculate the density of a multivariate tdistribution parameterised in terms of its mode, spectral decomposition (eigendecomposition) of its dispersion matrix and degrees of freedom.}
\usage{
dmvt.SpD(x, mu, Cov.SpD, df, logged = FALSE)
}
\arguments{
  \item{x}{An \code{n*d} matrix, where each row represents a point in \eqn{R^d} for which the density is required.
%%     ~~Describe \code{x} here~~
}
  \item{mu}{Vector of length \code{d}: the mode of the distribution}
  \item{Cov.SpD}{A list, with components \code{values} and \code{vectors}, defining the spectral decomposition of the dispersion matrix of the distribution. Typically the result of a call to \code{\link{eigen}}. The eigenvalues must all be strictly positive.}
  \item{df}{Numeric scalare: the degrees of freedom of the distribution}
  \item{logged}{Return the density (\code{FALSE}, the default) or its logarithm?}
}
\value{a vector of length \code{n}, containing the (log) densities.}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\note{Output from this routine has been cross-checked against that from 
\code{dmvt} in the \pkg{mvtnorm} library.}
