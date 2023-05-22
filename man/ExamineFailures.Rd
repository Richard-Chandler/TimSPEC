\name{ExamineFailures}
\alias{ExamineFailures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphics to examine failures of the Kalman Filter algorithm}
\description{To examine the parameter values for which \code{\link{SampleStates}} failed due to problems with the Kalman Filter, typically associated with near-singular matrix calculations.}
\usage{
ExamineFailures(Thetas, FailedIDs, alpha = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Thetas}{A matrix containing sampled parameter values, as produced by (e.g.) \code{\link{dlm.ThetaSample}} }
  \item{FailedIDs}{A vector containing the row numbers of \code{ThetaSamples} for which \code{\link{SampleStates}} failed. }
  \item{alpha}{Transparency to use for plotting non-failed samples. Values close to 1 are opaque, those close to 0 are more transparent and create the impression of a hazy background data cloud }
}
\details{The function produces a pairs plot of sampled parameter vectors, with the majority of the data points plotted as a background "data cloud" in grey; and with the points for failed samples plotted in the foreground in red. This is primarily intended for diagnostic purposes, hence few user-customisable options.}
\value{The function produces a plot on the current graphics device.}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{SampleStates}}}