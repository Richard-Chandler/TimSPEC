\name{CompactEns}
\alias{CompactEns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compacting an ensemble of time series}
\description{To produce a compact version of an ensemble of time series, by averaging members from the same simulator or simulator combination.
}
\usage{
CompactEns(Data, Groups, impute=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{Data frame containing an observed series in column 1 and ensemble members in the remaining columns. \strong{NB} column 1 is ignored: the structure of this argument is consistent with that of the \code{Data} argument for other routines in the \pkg{TimSPEC} package.
}
  \item{Groups}{A vector (numeric, character or factor) specifying the grouping structure in the ensemble, of length \code{ncol(Data)-1} corresponding to the ensemble-derived series in \code{Data}. Series with the same value of \code{Groups} are considered to come from the same source (e.g. the same GCM or GCM:RCM combination but under different initial conditions).
}
  \item{impute}{Controls whether to attempt to impute missing values in the ensemble first using \code{\link{EnsImpute}}}.
}

\details{The calculations assume that all series within a group contain non-missing values at the same set of time points. If this isn't the case, the routine will terminate with an error that some series have non-empty but incomplete groups. For an ensemble with such structure, either remove incomplete series or set \code{impute} to \code{TRUE} (the default) to solve the problem.
}

\value{The function returns an object of the same class as \code{Data}, but containing the observed series in column 1 and the group mean  series in subsequent columns. The result has an attribute \code{RunsPerTS}, which is a vector containing the numbers of series contributing to each of the group means. If \code{Data} itself has such an attribute, then the \code{RunsPerTS} attribute of the  result is obtained by summing the corresponding elements and the "group mean" series are weighted accordingly; otherwise it's just a count of the numbers of series being averaged.
}
\author{
Richard E. Chandler <r.chandler@ucl.ac.uk>
}
\seealso{\code{\link{EnsMean}}, \code{\link{PlotEnsTS}}, \code{\link{SSP585data}}}
\examples{
##  Load the example for the package data
example("SSP585data") # Note the creation of "GCMIDs" which defines groups

##
##  Produce mean series from each simulator (NB the first column
##  of the supplied data is Year which must be discarded)
##
CompactData <- CompactEns(GlobalTemps585[,-1], Groups=GCMIDs)

##
##  Numbers of series contributing to each of the group means
##  (compare with the table in the output from example("SSP585data") above)
##
print(attr(CompactData, "RunsPerTS"))

##
##  And plot (need to reinstate the "Year" column here). Compare 
##  with the example for PlotEnsTS
##
GCMCols <- c("CanESM5" = "darkgoldenrod1",
              "UKESM1-0-LL" = "chartreuse4",
              "MIROC-ES2L-f2" = "brown4",
              "NorESM2-MM" = "blue",
              "FGOALS-g3" = "violetred3")
SeriesIDs <- names(attr(CompactData, "RunsPerTS"))

PlotEnsTS(cbind(GlobalTemps585$Year, CompactData), 
          Colours=c("black", "darkred", GCMCols[SeriesIDs]), 
          main="Observed and simulated global mean temperatures, 1950-2100")
legend("topleft", lwd=c(3,3,rep(1, 5)), col=c("black", "darkred", GCMCols),
  legend=c("Observations", "Ensemble mean", names(GCMCols)), ncol=2)
}