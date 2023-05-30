\name{PlotEnsTS}
\alias{PlotEnsTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualising ensembles of climate projection time series}
\description{Plots observed time series (with a single value per year), individual ensemble members and ensemble mean. }
\usage{
PlotEnsTS(Data, Colours, EnsTransp = 1, Types = c(1, 1), 
          Units = expression(degree * "C"), plot = TRUE, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{Matrix or data frame containing time (i.e. year) in column 1, real-world series in column 2 and ensemble members in the remaining columns.}
  \item{Colours}{Vector containing the colour specifications for the observations and ensemble members respectively. If of length 2 then the two colours will be used to represent the observations and ensemble members, respectively. If of length \code{ncol(Data)} then the first element will be used for the observations, the second for the ensemble mean and the remainder for the individual ensemble members.}
  \item{EnsTransp}{Transparency (alpha) value to use when plotting the ensemble members. Values close to 1 are opaque, those close to 0 are more transparent.}
  \item{Types}{Vector containing the line types to use for observations and ensemble members respectively. Specification as for \code{Cols}, defaults to \code{c(1,1)} i.e. solid lines throughout. }
  \item{Units}{Length-1 character vector, used to label the vertical axis of the plot. Default value is appropriate for plots of temperature in degrees Celsius.}
  \item{plot}{If \code{TRUE} (the default) then plot the results; otherwise just set up the plot region but don't actually plot anything.}
  \item{add}{If \code{TRUE}, add lines to an existing plot; if \code{FALSE} (the default) then create a new plot. }
  \item{\dots}{Other arguments to \code{\link{matplot}}. }
}
\value{The function returns a \code{NULL} value: its purpose is to produce a plot, not to perform calculations for subsequent use. }
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{SmoothPlot}}}
