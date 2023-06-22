\name{which.diffuse}
\alias{which.diffuse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find diffuse elements of state vector in a dynamic linear model}
\description{Takes the transition matrix of a state space model and identifies which elements of the state vector are nonstationary. No checking is done; this is designed to be called from other functions.}
\usage{
which.diffuse(GG)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GG}{The transition matrix in the state equation of the model. }
}
\value{A logical vector that is \code{TRUE} for the nonstationary components of the state vector and \code{FALSE} for the stationary components.}
\references{
Chandler, R. E. and Scott, E. M. (2011). \emph{Statistical methods for trend detection and analysis in the environmental sciences}. John Wiley & Sons, Chichester.
}
\note{Users will rarely need to call this routine directly: it is used by other fitting and smoothing routines within the \pkg{TimSPEC} package.

The function is based on the software supplied with Chandler and Scott (2011).}
