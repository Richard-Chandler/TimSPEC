\name{make.SSinits}
\alias{make.SSinits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initialises a state space model}
\description{Computes initial values for the state vector and associated error covariances, for a time-invariant state space model.}
\usage{
make.SSinits(GG, W, kappa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GG}{Time-invariant transition matrix for the state process in the model}
  \item{W}{Innovation covariance matrix for the state process in the model}
  \item{kappa}{Variance used to initialise diffuse elements of the state vector.}
}
\details{The covariance matrix for stationary components of the state vector is calculated according to equation (3.3.21) of Harvey (1989).}
\value{A list containing components: 
\item{m0}{estimate of the state vector prior to seeing any data}
\item{P0}{associated error covariance matrix.}
}
\note{Users will rarely need to call this routine directly: it is called by the wrapper \code{...Smooth} functions.

The function is based on the software supplied with Chandler and Scott (2011).}
\references{
Chandler, R. E. and Scott, E. M. (2011). \emph{Statistical methods for trend detection and analysis in the environmental sciences}. John Wiley & Sons, Chichester.

Harvey, A.C. (1989). \emph{Forecasting, structural time series models and the Kalman filter}. Cambridge University Press, Cambridge
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
