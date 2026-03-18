\name{summary.dlmSample}
\alias{summary.dlmSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate posterior sampling for dynamic linear models}
\description{Generates samples from the approximate posterior distribution for the parameter vector of a dynamic linear model}
\usage{
\method{summary}{dlmSample}(object, weights=NULL, digits=4, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The result of a call to \code{\link{dlm.ThetaSample}} or \code{\link{PostPredSample}}}.
  \item{weights}{Either a logical scalar indicating whether or not to calculate a weighted summary using information in the \code{Weights} component of \code{object}, or a vector containing importance sampling weights (e.g. from a separate call to \code{\link{dlm.ImportanceWts}}. The default behaviour is to calculate the weighted summary if \code{object} contains a \code{Weights} component.}
  \item{digits}{The minimum number of significant digits to be used for screen output.}
  \item{...}{For consistency with the \code{summary()} generic.}
}

\value{A data frame containing a column for each of the parameters, and rows giving their posterior means and standard deviations. If \code{object} contains a non-\code{NULL} \code{Weights} component and \code{weights} is \code{TRUE}, or if \code{weights} is a vector, then the data frame contains rows for both weighted and unweighted means and standard deviations (so four rows in total). Otherwise it contains just the unweighted versions.

The function also displays the summary on the standard output device.}

\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

\seealso{\code{\link{PostPredSample}} for an example; \code{\link{dlm.ThetaSample}}.}
