\name{CumulativeWeightPlot}
\alias{CumulativeWeightPlot}
\title{Plotting of importance sampling weights}
\description{To plot the cumulative contributions of a vector of importance sampling weights, sorted in descending order}
\usage{
CumulativeWeightPlot(w, main="Contributions of importance sampling weights", ...)
}
\arguments{
  \item{w}{The vector of weights}
  \item{main}{Title for the plot}
  \item{\dots}{Other arguments to \code{\link{plot.stepfun}}
}
}
\details{The function is invoked primarily to produce a plot which can be used as a diagnostic when using importance sampling. If the weights are all roughly equal then the plot will be close to a diagonal line (in which case importance sampling conveys little benefit); if they are dominated by one or two extreme values, it will pass close to the top left-hand corner of the plotting region (suggesting that the weights may have infinite variance so that importance sampling is likely to fail).
}
\value{The result from \command{plot.stepfun} is returned, invisibly.}
\references{Monahan, J. F. (2001). \emph{Numerical Methods of Statistics}, Chapter 12. Cambridge University Press, Cambridge.
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{plot.stepfun}}}
