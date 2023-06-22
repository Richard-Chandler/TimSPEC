\name{SLLT.IniPar}
\alias{SLLT.IniPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Preliminary parameter estimates for smooth local linear trend model}
\description{To find initial values for the two variance parameters in a smooth local linear trend (SLLT) model, as described in \code{\link{SLLT.modeldef}}. }
\usage{
SLLT.IniPar(Y, method = "arima", collapse = TRUE, Tiny = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Either a vector or a matrix containing observed series. If a matrix, the columns are considered to contain the series. }
  \item{method}{Either \code{arima} (the default) or \code{moments}. \code{arima} exploits the fact that the SLLT model is stochastically equivalent to an ARIMA(0,2,2) process and obtains the required values using the \code{CSS-ML} algorithm inplemented by \code{\link{arima}}. \code{moments} uses a method of moments based on the variance and lag-1 covariance of the second differences of the series. }
  \item{collapse}{If \code{Y} is a matrix then, if \code{collapse} is \code{TRUE}, a single pair of initial values is produced using the data from all columns; if \code{collapse} is \code{FALSE} then values are produced for each column separately. }
  \item{Tiny}{A small value used to prevent negative variance estimates. }
}
\value{if \code{Y} is a vector or if \code{collapse} is \code{TRUE}, a vector of length 2 containing initial values for the logarithms of the measurement and innovation variances in the SLLT model. If \code{Y} is a matrix and \code{collapse} is \code{FALSE}, a 2-row matrix containing separate pairs for each column of \code{Y}. }
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\note{Although the function provides the option of using a method of moments to find initial values by setting \code{method="moments"}, this is very inefficient for ARIMA  models containing a moving average component: its only advantage over \code{method="arima"} is that it's a bit quicker, although experience suggests that it can indeed produce wildly inaccurate initial values. }

\note{Users will rarely need to call this routine directly: it is called by other fitting and smoothing functions.}
\seealso{\code{\link{SLLT.modeldef}}, \code{\link{SLLTSmooth}}. }
