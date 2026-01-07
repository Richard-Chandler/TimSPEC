\name{EnsMean}
\alias{EnsMean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating ensemble mean time series}
\description{To compute the mean of a collection of time series from 
structured or unstructured ensembles.
}
\usage{
EnsMean(Data, Groups=NULL, na.option=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{Data frame containing an observed series in column 1 and ensemble members in the remaining columns. \strong{NB} column 1 is ignored: the structure of this argument is consistent with that of the \code{Data} argument for other routines in the \pkg{TimSPEC} package.
}
  \item{Groups}{A vector or matrix specifying any structure in the ensemble, with an element (resp. row, in the matrix case) for each ensemble member and containing the codes of any grouping variables defining the structure. This corresponds to the \code{Groups} argument in, for example, \code{\link{EnsEBMtrendSmooth}} or \code{\link{EnsEBM2waytrendSmooth}}.
}
  \item{na.option}{Controls the treatment of missing ensemble members. Valid non-\code{NULL} values are \code{"na.rm"} (remove \code{NA} cases prior to the calculations), \code{"na.impute"} (see below) or \code{"na.pass"} (propagate the \code{NA}s through the calculations). The default \code{NULL} value corresponds to \code{"na.rm"} if \code{Groups} is \code{NULL}, and \code{"na.impute"} otherwise. 
}
}

\details{If \code{Groups} is \code{NULL} (the default) then this function computes the overall mean of the ensemble members (subject to the treatment of missing values according to the value of \code{na.option}). Otherwise, the result is formed from the intercepts in ANOVA models fitted separately to the ensemble members at each time point, with factors defined by the column(s) of \code{EnsStruct} and using sum-to-zero constraints: in this case, the result can be interpreted as a mean of the simulators or simulator combinations.

If \code{na.option} is \code{"na.rm"} then the missing value patterns are processed separately at each time point. This can lead to discontinuities in the calculated ensemble mean series if the ensemble members cover different periods, especially in structured ensembles where one simulator has been run for a slightly different time period from the remainder. In this case, the \code{"na.impute"} option may be preferable. With this option, means are computed after imputing any missing values using an additive ANOVA model with effects corresponding to ensemble member and to the time variable (treated as a factor). For details, see the help for \code{\link{EnsImpute}}.
}

\value{A vector containing the computed series.}
\author{
Richard E. Chandler <r.chandler@ucl.ac.uk>
}
\note{This function is intended to be called primarily from other routines in the \pkg{TimSPEC} package.}
