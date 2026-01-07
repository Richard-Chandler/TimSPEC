\name{dlm.ImportanceWts}
\alias{dlm.ImportanceWts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Importance weights for DLMs}
\description{Calculates importance weights for sampling from the posterior distribution of the parameter vector in a dynamic linear model}
\usage{
dlm.ImportanceWts(samples, build, Y, prior.pars = NULL, debug = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{A matrix, containing samples from the Gaussian (or multivariate \eqn{t}) approximation to the posterior and with attributes \code{Mean}, \code{Cov.SpD} and \code{df} as produced by \code{\link{dlm.ThetaSample}}
}
  \item{build}{A function that will construct the \code{dlm} object representing the model being fitted. Typically one of the \code{.modeldef} functions in the \pkg{TimSPEC} package. 
}
  \item{Y}{The data: a matrix containing a column for each time series in the ensemble, the first column corresponding to the real climate system. 
}
  \item{prior.pars}{As in \code{\link{dlm.SafeLL}}
}
  \item{debug}{As in \code{\link[dlm]{dlmLL}}
}
  \item{\dots}{Other arguments to \code{build}. Note, however, that if \code{build} has a \code{RunsPerTS} argument then the routine will take this from the \code{RunsPerTS} attribute of \code{Y} if present e.g. arising from a call to \code{\link{CompactEns}}: in this case therefore, the user should not supply this argument separately.}
}
\details{
The data matrix \code{Y} is used to calculate the exact log-posterior for each sampled value, which is then compared with the log-density used to produce \code{samples} in order to obtain importance weights.

The routine will issue a warning if the log posterior density of the any of the sampled values have a log posterior greater than that obtained from the \code{Mean} component of \code{samples}, because this indicates that the numerical optimisation failed when computing the MAP estimates from which the samples were derived. 
}
\value{A data frame containing three columns, and a row for each row of \code{samples}. The columns are \code{log.g} containing the values of \eqn{\log(g(\theta)/g(\hat{\theta}))}{log(g(theta) / g(theta.hat))} where \eqn{\hat{\theta}}{theta.hat} is the posterior mode; \code{log.h} containing the values of \eqn{\log(h(\theta)/h(\hat{\theta}))}{log(h(theta)/h(theta.hat))}; and the corresponding normalised importance weights. Here, \eqn{g(.)} is the density used to produce \code{samples} and \eqn{h(.)} is the posterior density. 
}
\references{
Durbin, J. and Koopman, S. J. (2012). \emph{Time series analysis by state space methods (second
edition)}, Chapter 13. Oxford University Press, Oxford.
}
\seealso{\code{\link{CumulativeWeightPlot}}, \code{\link{PostPredSample}}}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>
}
