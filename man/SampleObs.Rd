\name{SampleObs}
\alias{SampleObs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sampling from posterior predictive distribution.}
\description{To sample from the posterior predictive distribution of (elements of) the observation vector for a dynamic linear model, conditioned on a set of sampled state vectors.}
\usage{
SampleObs(Thetas, States, build, Y, WhichEls = 1:ncol(Y), NonNeg = FALSE, ReplaceAll = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Thetas}{A matrix containing sampled parameter values, as produced by (e.g.) \code{\link{dlm.ThetaSample}}. }
  \item{States}{An array containing sampled time series of the state vector for each row of \code{Thetas} produced (e.g.) by \code{\link{SampleStates}}. }
  \item{build}{Function to construct the \code{dlm} object being used. Typically one of the \code{.modeldef} functions in the \pkg{PECTSP} package. }
  \item{Y}{The data: a matrix containing a column for each time series, in the form anticipated by the model constructed by \code{build}.}
  \item{WhichEls}{Numeric vector indicating which column(s) of \code{Y} to simulate. }
  \item{NonNeg}{Indicates whether the values to be returned should be non-negative. If not (the default), sampling is done from a normal distribution. Otherwise the samples are then transformed to the corresponding quantiles of a gamma distribution with the same mean and variance. }
  \item{ReplaceAll}{If \code{TRUE}, the returned object will contain random samples for all time points. If \code{FALSE} (the default), only the values missing from \code{Y} are simulated (and the non-missing values are included in the output). }
  \item{\dots}{Additional arguments to \code{build}. }
}
\details{Most users won't need to call this function: use \code{\link{PostPredSample}} instead, which is a wrapper to all of the function calls needed to produce a posterior predictive sample. }
\value{The function returns an array of dimension \code{c(T, p, N)} where \code{N} is the number of samples, \code{T} the number of time points in \code{Y}, and \code{p} the number of series being simulated (i.e. the length of \code{WhichEls}). The ordering of dimensions ensures that if a single series is being simulated then the result will be a matrix with the series in columns. }
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{PostPredSample}}}
