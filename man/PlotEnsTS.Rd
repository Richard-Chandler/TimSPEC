\name{PlotEnsTS}
\alias{PlotEnsTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualising ensembles of climate projection time series}
\description{Plots observed time series (with a single value per year), individual ensemble members and ensemble mean. }
\usage{
PlotEnsTS(Data, Colours, EnsStruct=NULL, na.rm=TRUE, EnsTransp = 1,  
          Types = c(1, 1), Units = expression(degree * "C"), 
          plot = TRUE, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{Matrix or data frame containing time (i.e. year) in column 1, real-world series in column 2 and ensemble members in the remaining columns.}
  \item{Colours}{Vector containing the colour specifications for the observations and ensemble members respectively. If of length 2 then the two colours will be used to represent the observations and ensemble members, respectively. If of length \code{ncol(Data)} then the first element will be used for the observations, the second for the ensemble mean and the remainder for the individual ensemble members.}
  \item{EnsStruct}{Vector or matrix representing ensemble structure, used to control calculation of the ensemble mean series. See \code{\link{EnsMean}} for details.}
  \item{na.rm}{Controls handling of missing values when calculating the ensemble mean series. See \code{\link{EnsMean}} for details.}
  \item{EnsTransp}{Transparency (alpha) value to use when plotting the ensemble members. Values close to 1 are opaque, those close to 0 are more transparent.}
  \item{Types}{Vector containing the line types to use for observations and ensemble members respectively. Specification as for \code{Cols}, defaults to \code{c(1,1)} i.e. solid lines throughout. }
  \item{Units}{Length-1 character vector, used to label the vertical axis of the plot. Default value is appropriate for plots of temperature in degrees Celsius.}
  \item{plot}{If \code{TRUE} (the default) then plot the results; otherwise just set up the plot region but don't actually plot anything.}
  \item{add}{If \code{TRUE}, add lines to an existing plot; if \code{FALSE} (the default) then create a new plot. }
  \item{\dots}{Other arguments to \code{\link{matplot}}. }
}
\value{The function returns a \code{NULL} value: its purpose is to produce a plot, not to perform calculations for subsequent use. }
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

\seealso{\code{\link{SmoothPlot}}, \code{\link{SSP585data}}}
\examples{
## Load package data 
data(SSP585data)

## Plot with simple colour scales
PlotEnsTS(GlobalTemps585, Colours=c("black", "darkred"), EnsTransp=0.2, 
  main="Observed and simulated global mean temperatures, 1950-2100")
legend("topleft", lwd=c(3,3,1), col=c("black", "darkred", "darkred"),
  legend=c("Observations", "Ensemble mean", "Ensemble members"))

## And with each climate model identified (kudos to Clair Barnes
## for the colour choices)
GCMCols <- c("CanESM5" = "darkgoldenrod1",
              "UKESM1-0-LL" = "chartreuse4",
              "MIROC-ES2L-f2" = "brown4",
              "NorESM2-MM" = "blue",
              "FGOALS-g3" = "violetred3")
GCMIDs <- 
  unlist(lapply(strsplit(names(GlobalTemps585)[-(1:2)], "\\\\."), 
                FUN=function(x) x[1]))

PlotEnsTS(GlobalTemps585, 
  Colours=c("black", "darkred", GCMCols[GCMIDs]),  
  main="Observed and simulated global mean temperatures, 1950-2100")
legend("topleft", lwd=c(3,3,rep(1, 5)), col=c("black", "darkred", GCMCols),
  legend=c("Observations", "Ensemble mean", names(GCMCols)), ncol=2)

## And using the EnsStruct argument to represent the mean of the GCMs
## rather than that of the ensemble members

PlotEnsTS(GlobalTemps585, EnsStruct=GCMIDs,
  Colours=c("black", "darkred", GCMCols[GCMIDs]),  
  main="Observed and simulated global mean temperatures, 1950-2100")
legend("topleft", lwd=c(3,3,rep(1, 5)), col=c("black", "darkred", GCMCols),
  legend=c("Observations", "Ensemble mean", names(GCMCols)), ncol=2)

}
