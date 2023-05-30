\name{EBMtrendSmooth}
\alias{EBMtrendSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postprocessing a single time series using an EBM-inspired dynamic linear model}
\description{Fits a model of the form defined by \code{\link{EBMtrend.modeldef}}, then applies the Kalman Smoother to obtain an estimate of the underlying trend with prediction intervals.}
\usage{
EBMtrendSmooth(Y, Xt, m0 = NULL, kappa = 1e+06, UsePhi = TRUE, 
               prior.pars = NULL, messages = TRUE, Use.dlm = FALSE, 
               debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Vector containing a univariate time series.}
  \item{Xt}{Vector containing the corresponding forcings (denoted as \eqn{f(t)} in the documentation for \code{\link{EBMtrend.modeldef}}).}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{kappa}{Initialisation variance for diffuse elements of the state vector.}
  \item{UsePhi}{Controls whether or not to use the full model including a thermal inertia parameter.}
  \item{prior.pars}{Optional \code{3*2} matrix. If present, it is taken to contain the means and standard deviations of independent Gaussian priors for the transformed model parameters (log variances, and logit thermal inertia).}
  \item{messages}{Controls whether to print progress to screen.}
  \item{Use.dlm}{Logical scalar, used to control the method used to fit the model: see \code{\link{dlm.SafeMLE}}.} 
  \item{debug}{As in \code{\link{dlmMLE}}. Only \code{TRUE} is currently supported for this model, until the memory leak in the \pkg{dlm} package is fixed.
}
}
\details{The function is a wrapper for the steps required to fit and apply a dynamic linear model. It starts by trying to find reasonable values for (penalised) maximum likelihood estimation of the model parameters; then carries out the estimation; and finally applies the Kalman Smoother to the data using the fitted model. 
}
\value{A named list with three components:
\item{Theta}{Object containing parameter estimation results; the result of a call to \code{\link{dlmMLE}} or \code{\link{dlm.SafeMLE}}}
\item{Model}{An object of class \code{dlm} containing the fitted model itself}
\item{Smooth}{An object containing the result of Kalman Smoothing the input series using the fitted model. This is produced using \code{\link{dlmSmooth}}}
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\note{\code{Use.dlm=TRUE} is currently unsupported for this model; and the routine currently contains bespoke optimisation code instead of calling \code{\link{dlm.SafeMLE}} due to a memory leak in the \code{C} code for \code{\link{dlmLL}}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{EnsEBMtrend.modeldef}}}
