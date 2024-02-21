\name{SLLTSmooth}
\alias{SLLTSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth Local Linear Trend modelling for a single time series}
\description{Fits and applies a "smooth local linear trend" (SLLT) model to a univariate time series.}
\usage{
SLLTSmooth(Y, m0 = NULL, C0=NULL, kappa = 1e+06, prior.pars = NULL, 
           messages = TRUE, Use.dlm = FALSE, debug = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Vector containing a univariate time series.}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{Variance used to initialise diffuse elements of the state vector.}
  \item{prior.pars}{Optional \code{2*2} matrix containing means and standarddeviations of Gaussian prior distributions for the log variance parameters: in this case, maximum \emph{a posteriori} (MAP) estimation is done.}
  \item{messages}{Controls whether to print progress to screen. For details, see the help for \code{\link{dlm.SafeMLE}}.}
  \item{Use.dlm}{Logical scalar, used to control the method used to fit the model: see \code{\link{dlm.SafeMLE}}.}
  \item{debug}{As in \code{\link{dlmMLE}}.}
  \item{\dots}{Additional arguments to \code{\link{dlm.SafeMLE}}}
}
\details{The function is a wrapper for the steps required to fit and apply a dynamic linear model. It starts by trying to find reasonable values for (penalised) maximum likelihood estimation of the model parameters; then carries out the estimation; and finally applies the Kalman Smoother to the data using the fitted model.}
\value{A named list with three components:
\item{Theta}{Object containing parameter estimation results; the result of a call to \code{\link{dlmMLE}} or \code{\link{dlm.SafeMLE}}}
\item{Model}{An object of class \code{dlm} containing the fitted model itself}
\item{Smooth}{An object containing the result of Kalman Smoothing the input series using the fitted model. This is produced using \code{\link{dlmSmooth}}}
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

\seealso{\code{\link{SLLT.modeldef}}, \code{\link{EnsSLLTSmooth}}, \code{\link{dlm.SafeMLE}}, \code{\link{summary.dlmMLE}}, \code{\link{SmoothPlot}}}
\examples{
## Load package data 
data(SSP585data)

## Prior means and standard deviations for log variances
priors <- rbind(c(-3.7, 3.45),    # log(sigsq.0)
                c(-8.3, 3.45))    # log(tausq.0)

## And prior mean and variance for trend & slope in 1949
m0 <- c(14, 0)
k <- 25

## Fit and smooth the observations, then plot. NB second element of
## PredColours isn't used here (ensemble data not used in smoothing)
SLLTFit <- 
  SLLTSmooth(GlobalTemps585$Obs, m0=m0, kappa=k, prior.pars=priors)
SmoothPlot(GlobalTemps585[,1:2], SLLTFit, DatColours=c("black", NA),
           PredColours=c("darkblue", NA), ylim=c(13.5,19.5))

## Here's a plot of the smooth with the ensemble members as well
SmoothPlot(GlobalTemps585, SLLTFit, DatColours=c("black", "coral3"),
           PredColours=c("darkblue", NA), EnsTransp=0.3)

}
