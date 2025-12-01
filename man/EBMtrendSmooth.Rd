\name{EBMtrendSmooth}
\alias{EBMtrendSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postprocessing a single time series using an EBM-inspired dynamic linear model}
\description{Fits a model of the form defined by \code{\link{EBMtrend.modeldef}}, then applies the Kalman Smoother to obtain an estimate of the underlying trend with prediction intervals.}
\usage{
EBMtrendSmooth(Y, Xt, m0 = NULL, C0=NULL, kappa = 1e+06, UsePhi = TRUE, 
               prior.pars = NULL, theta=NULL, messages = TRUE, 
               Use.dlm = FALSE, debug = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Vector containing a univariate time series.}
  \item{Xt}{Vector containing the corresponding forcings (denoted as \eqn{f(t)} in the documentation for \code{\link{EBMtrend.modeldef}}).}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
    \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{Initialisation variance for diffuse elements of the state vector.}
  \item{UsePhi}{Controls whether or not to use the full model including a thermal inertia parameter.}
  \item{prior.pars}{Optional \code{3*2} matrix. If present, it is taken to contain the means and standard deviations of independent Gaussian priors for the transformed model parameters (log variances, and logit thermal inertia).}
  \item{theta}{Optional initial value for parameter vector. If \code{NULL}, the routine will auto-initialise the optimisation in the first step.}
  \item{messages}{Controls whether to print progress to screen. For details, see the help for \code{\link{dlm.SafeMLE}}.}
  \item{Use.dlm}{Logical scalar, used to control the method used to fit the model: see \code{\link{dlm.SafeMLE}}.} 
  \item{debug}{As in \code{\link[dlm]{dlmMLE}}.}
  \item{...}{Other arguments to \code{\link{dlm.SafeMLE}}.}
}
\details{The function is a wrapper for the steps required to fit and apply a dynamic linear model. It starts by trying to find reasonable values for (penalised) maximum likelihood estimation of the model parameters; then carries out the estimation; and finally applies the Kalman Smoother to the data using the fitted model. 
}
\value{A named list with three components:
\item{Theta}{Object containing parameter estimation results; the result of a call to \code{\link[dlm]{dlmMLE}} or \code{\link{dlm.SafeMLE}}}
\item{Model}{An object of class \code{dlm} containing the fitted model itself}
\item{Smooth}{An object containing the result of Kalman Smoothing the input series using the fitted model. This is produced using \code{\link[dlm]{dlmSmooth}}}
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{EBMtrend.modeldef}}, \code{\link{EnsEBMtrendSmooth}}, \code{\link{EnsEBM2waytrendSmooth}}, \code{\link{dlm.SafeMLE}}, \code{\link{summary.dlmMLE}}}

\examples{
## Load package data 
data(SSP585data)

## Prior means and standard deviations for log variances
priors <- rbind(c(-3.7, 3.45),    # log(sigsq.0)
                c(-8.3, 3.45),    # log(tausq.0)
                c(0, 5))          # logit(phi)
## And prior mean and variance for trend, drift & ERF coefficient in 1949
m0 <- c(14, 0, 1)
k <- 25

## Fit and smooth the observations, then plot. NB second element of
## PredColours isn't used here (ensemble data not used in smoothing)
EBMFit <- EBMtrendSmooth(GlobalTemps585$Obs, Xt=ERF585$NetERF,
                         kappa=k, UsePhi=TRUE, prior.pars=priors)
SmoothPlot(GlobalTemps585[,1:2], EBMFit, DatColours=c("black", NA),
           PredColours=c("darkblue", NA), ylim=c(13.5,18))

## Here's a plot of the smooth with the ensemble members as well
SmoothPlot(GlobalTemps585, EBMFit, DatColours=c("black", "coral3"),
           PredColours=c("darkblue", NA), EnsTransp=0.3)


}
