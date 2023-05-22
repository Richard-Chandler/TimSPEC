\name{plot.ImportanceWts}
\alias{plot.ImportanceWts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualise importance sampling weights}
\description{To produce a pairs plot of sampled parameter values from a dynamic linear model, with points shaded according to the value of their importance sampling weights.}
\usage{
plot.ImportanceWts(samples, weights, colours = NULL, oma = c(7, 3, 5, 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{A matrix containing sampled parameter values, with attributes as set by (e.g.) \code{\link{dlm.ThetaSample}} }
  \item{weights}{A 3-column matrix in which the third column contains importance weights for each row of samples, as produced by (e.g.) \code{\link{dlm.ImportanceWts}} }
  \item{colours}{A vector of colours to use when shading the plots. If \code{NULL} (the default), a blue-red scale is used: dark blue denotes points with very low importance weights and dark red denotes those with very high weights. }
  \item{oma}{This is a fudge, used to define the "outer margins" of the plot area in order to put a legend at the bottom (the underlying \code{\link{pairs}} command doesn't allow for this). The default value may need to be tweaked, depending on the dimensions of the output device. }
  \item{\dots}{Other arguments to \code{\link{pairs}}. }
}
\value{The function produces a plot on the current graphics device}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{dlm.ThetaSample}}, \code{\link{dlm.ImportanceWts}}, \code{\link{CumulativeWeightPlot}} }
