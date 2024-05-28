\name{SmoothPlot}
\alias{SmoothPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots of postprocessing results}
\description{To plot the smoothing results from a dynamic linear model (DLM), together with the data used to produce them.}
\usage{
SmoothPlot(Data, Smooth = NULL, Samples = NULL, AntiLog = FALSE, 
           NonNeg = FALSE, PlotMu0 = TRUE, PlotConsensus = FALSE, 
           Weights = NULL, DatColours, DatTypes = c(1, 1), PredColours, 
           alpha = c(0.6, 0.1), SmoothPIs = FALSE, 
           Units = expression(degree * "C"), plot.title = "", 
		   Legend=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{A data frame or matrix containing time (years) in column 1, real-world series in column 2 and ensemble members in the remaining columns. }
  \item{Smooth}{A list containing at least the named components \code{Model} (an object of class \code{dlm}) and \code{Smooth} (the result of Kalman Smoothing a dataset using the DLM stored in \code{Model}, via  \code{\link{dlmSmooth}}. Typically produced via one of the \code{*Smooth} functions in the \pkg{TimSPEC} package. Ignored (with a warning) if \code{Samples} is non-\code{NULL}. }
  \item{Samples}{List containing at least the named elements \code{States} (an array containing posterior samples of the state vector) \code{Obs} (a matrix containing samples from the posterior predictive distribution of the observed series). Such a list will usually be produced using \code{\link{PostPredSample}}. If non-\code{NULL}, it is used to compute prediction intervals and posterior means of relevant quantities. }
  \item{AntiLog}{If \code{TRUE}, all results will be exponentiated(with appropriate adjustments for the posterior means) prior to plotting. This allows for situations where the DLM has been fitted to logged data (e.g. to ensure that postprocessed precipitation samples are guaranteed to be non-negative). }
  \item{NonNeg}{If \code{FALSE} (the default), prediction intervals will be computed in the usual way as mean +/- 1.96 standard deviations --- these being the 2.5th and 97.5th percentiles of a normal distribution. If \code{TRUE}, the quantiles of a gamma distribution will be used instead to ensure that the intervals are non-negative. }
  \item{PlotMu0}{Logical scalar controlling whether or not to plot the estimated trend in the real-world series. }
  \item{PlotConsensus}{Controls whether to add the ensemble consensus to the plot. }
  \item{Weights}{Optional vector of importance weights associated with the samples in \code{Samples}. Ignored if \code{Samples} itself contains an element named \code{Weights} (e.g. if it was produced via a call to \code{PostPredSample}). }
  \item{DatColours}{Vector of colours used to plot the ensemble members: for details, see header \code{\link{PlotEnsTS}}. }
  \item{DatTypes}{Vector of line types, similarly. }
  \item{PredColours}{Vector of two colours: the first is the base colour for plotting estimates of the observed trend \eqn{\mu_{0}(\cdot)}{mu[0](.)}, and the second for plotting the ensemble consensus. }
  \item{alpha}{Vector of two transparency values used for plotting estimates and prediction intervals respectively. Values close to 1 are opaque, those close to 0 are more transparent. }
  \item{SmoothPIs}{If \code{TRUE} and if \code{Samples} is provided, the prediction intervals on the plot will be smoothed using \code{\link{loess}}. Otherwise (the default) the raw quantiles from the sampled serieswill be plotted. }
  \item{Units}{Length-1 character vector, used to label the vertical axis of the plot. Default value is appropriate for plots of temperature in degrees Celsius. }
  \item{plot.title}{Self-explanatory, hopefully. }
  \item{Legend}{Logical scalar controlling whether or not to add a legend to the plot. Default is to add it.}
  \item{\dots}{Other arguments to \code{\link{PlotEnsTS}}. }
}
\value{The function produces a plot on the current graphics device. It also returns, invisibly, a data frame containing the values of the estimated "real" and "ensemble consensus" trends, together with the upper and lower 95\% limits for the real trend. The interpretation of the limits depends on the inputs: for a \code{Smooth} object, they're pointwise confidence bands while for a \code{Samples} object they are pointwise posterior predictive intervals. }
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

\seealso{\code{\link{SLLTSmooth}}, \code{\link{EnsSLLTSmooth}}, \code{\link{EBMtrendSmooth}} for examples. Also \code{\link{PlotEnsTS}}, \code{\link{PostPredSample}}. }
