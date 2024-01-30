\name{EnsMean}
\alias{EnsMean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating ensemble mean time series}
\description{To compute the mean of a collection of time series from 
structured or unstructured ensembles.
}
\usage{
EnsMean(Data, EnsStruct=NULL, na.rm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{Data frame containing time in column 1, an observed series in column 2 and ensemble members in the remaining columns. \strong{NB} columns 1 and 2 is ignored: the structure of this argument is consistent with that of the \code{Data} argument for other routines in the \pkg{TimSPEC} package.
}
  \item{EnsStruct}{A vector or matrix specifying any structure in the         ensemble, with an element (resp. row, in the matrix case) for each ensemble member and containing the codes of any grouping variables defining the structure. In the matrix case, each column corresponds to a different grouping variable. The intention is that the \code{Groups} argument to, for example, \code{\link{EnsEBMtrendSmooth}} or \code{\link{EnsEBM2waytrendSmooth}} could be used here.
}
  \item{na.rm}{Remove missing ensemble members when calculating the means? If \code{TRUE} (the default), this is done separately for each time point - hence allowing for ensembles including time series of different lengths. Otherwise, the calculated means are \code{NA} at any time point with at least one missing ensemble member. 
}
}

\details{If \code{EnsStruct} is \code{NULL} (the default) then this function computes the overall mean of the ensemble members (subject to the treatment of missing values according to the value of \code{na.rm}). Otherwise, the result is formed from the intercepts in ANOVA models fitted separately to the ensemble members at each time point, with factors defined by the column(s) of \code{EnsStruct} and using sum-to-zero constraints: in this case, the result can be interpreted as a mean of the simulators or simulator combinations.}

\value{A vector containing the computed series.}
\author{
Richard E. Chandler <r.chandler@ucl.ac.uk>
}
\note{This function is intended to be called primarily from other routines in the \pkg{TimSPEC} package.}
