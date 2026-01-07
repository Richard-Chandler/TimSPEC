\name{SampleStates}
\alias{SampleStates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sampling from posterior distribution of state vector}
\description{To sample from the posterior distribution of the state vector for a dynamic linear model. This is essentially a wrapper for \code{\link[dlm]{dlmBSample}}. }
\usage{
SampleStates(Thetas, build, Y, NonNeg = NULL, debug = FALSE, 
             messages = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Thetas}{A matrix containing sampled parameter values, as produced by (e.g.) \code{\link{dlm.ThetaSample}} }
  \item{build}{A function that will construct the \code{dlm} object representing the model being fitted. Typically one of the \code{.modeldef} functions in the \code{TimSPEC} package. }
    \item{Y}{The data: a matrix containing a column for each time series, in the form anticipated by the model constructed by \code{build}. This is used to calculate the exact log-posterior for each row of \code{Thetas}. }
  \item{NonNeg}{An optional vector containing indices of state elements that should not be negative (e.g. when analysing untransformed precipitation series). }
  \item{debug}{As in \code{\link[dlm]{dlmLL}}. }
  \item{messages}{Controls whether to print progress to screen. }
  \item{\dots}{Additional arguments to \code{build}. Note, however, that if \code{build} has a \code{RunsPerTS} argument then the routine will take this from the \code{RunsPerTS} attribute of \code{Y} if present e.g. arising from a call to \code{\link{CompactEns}}: in this case therefore, the user should not supply this argument separately.}
}
\details{Most users won't need to call this function: use \code{\link{PostPredSample}} instead, which is a wrapper to all of the function calls needed to produce a posterior predictive sample. }

\value{The function returns an array of dimension \code{c(N, T, p)} where \code{N} is the number of samples, \code{T} the number of time points in \code{Y}, and \code{p} the dimension of the state vector. Samples for which the Kalman Filter failed are left as \code{NA}, with a warning, as are samples containing negative values in the columns selected by \code{NonNeg} if this is non-\code{NULL}. In this case, the result also has a \code{FailedIDs} attribute which is a numeric vector containing the indices of the failed samples. 

\strong{NOTE} that the routine removes the initial "time zero" values inserted by the \pkg{dlm} package. }
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{PostPredSample}}}
