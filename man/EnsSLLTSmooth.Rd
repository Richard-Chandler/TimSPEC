\name{EnsSLLTSmooth}
\alias{EnsSLLTSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postprocessing an exchangeable ensemble of time series using smooth local linear trend models}
\description{Fits a model of the form defined by \code{\link{EnsSLLT.modeldef}}, then applies the Kalman Smoother to obtain an estimate of the underlying trends with prediction intervals.}
\usage{
EnsSLLTSmooth(Y, m0 = NULL, C0=NULL, kappa = 1e+06, discrepancy = "varying", 
              Groups = NULL, UseAlpha = TRUE, prior.pars = NULL, 
              theta=NULL, constrain = TRUE, Tiny = 1/kappa, ObsSmooth, 
              Ens0Theta, messages = TRUE, Use.dlm = FALSE, debug = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Matrix containing an observed time series in its first column and exchangeable series from an ensemble in the remaining columns}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{Variance used to initialise diffuse elements of the state vector.}
  \item{discrepancy}{Either "varying" (the default) or "constant". If "constant" then the discrepancies between each ensemble member and reality are constant through time; otherwise they evolve as random walks.}
  \item{Groups}{Vector indicating group membership for each ensemble member. This should be used if an ensemble contains multiple runs from one or more simulators (climate models): in this case the elements of \code{Groups} should be integers between 1 and \code{G} where \code{G} is the total number of simulators contributing to the ensemble: \code{Groups[i]} is the index of the simulator that was used to produce the \eqn{i}th ensemble member. If \code{NULL} (the default), each member is assumed to come from a different simulator. 
}
  \item{UseAlpha}{Logical scalar indicating whether or not to include a scaling factor (\eqn{\alpha}{alpha} in the "Details" section of \code{\link{EnsEBMtrend.modeldef}}) when relating the ensemble consensus trend to that in the real climate system.  If \code{FALSE}, the value of \eqn{\alpha}{alpha} is set to 1.}
  \item{prior.pars}{An optional \code{4*2} (if \code{discrepancy} is \code{"constant"}) or \code{6*2} (if \code{discrepancy} is \code{"varying"}) matrix containing means and standard deviations of independent Gaussian priors for the log variance parameters in the model. If this is provided, maximum \emph{a posteriori} (or penalised maximum likelihood) estimation is used; otherwise standard maximum likelihood.}
  \item{theta}{Optional initial value for parameter vector. If \code{NULL}, the routine will auto-initialise the optimisation in the first step.}
  \item{constrain}{Logical scalar controlling whether to impose sum-to-zero constraints across the ensemble members where necessary to make the model identifiable.}
  \item{Tiny}{A small positive value, used to replace negative initial estimates for variances (see \code{\link{SLLT.IniPar}})}
  \item{ObsSmooth}{Optional: result of a previous call to \code{\link{SLLTSmooth}} for the first column of \code{Y}. If present, this is used to help initialise the numerical search for an MLE if a "constant" model needs to be fitted, or if the resulting smooth is needed. Otherwise the relevant model will be fitted as part of this routine.}
  \item{Ens0Theta}{Optional parameter vector for a previous fit of an "ensemble local linear trend" model with a constant discrepancy term. If present, this will be used to help initialise the numerical search for an MLE if a "varying" model needs to be fitted here. Otherwise the relevant constant-discrepancy model will be fitted as part of this routine.}
  \item{messages}{Controls whether to print progress to screen. For details, see the help for \code{\link{dlm.SafeMLE}}.}
  \item{Use.dlm}{Logical scalar, used to control the method used to fit the model: see \code{\link{dlm.SafeMLE}}.}
  \item{debug}{As in \code{\link[dlm]{dlmMLE}}.}
  \item{...}{Other arguments to \code{\link{dlm.SafeMLE}}.}
}
\details{The function is a wrapper for the steps required to fit and apply a dynamic linear model. It starts by trying to find reasonable values for (penalised) maximum likelihood estimation of the model parameters; then carries out the estimation; and finally applies the Kalman Smoother to the data using the fitted model.}
\value{A named list with three components:
\item{Theta}{Object containing parameter estimation results; the result of a call to \code{\link{dlmMLE}} or \code{\link{dlm.SafeMLE}}}
\item{Model}{An object of class \code{dlm} containing the fitted model itself}
\item{Smooth}{An object containing the result of Kalman Smoothing the input series using the fitted model. This is produced using \code{\link{dlmSmooth}}}
}

\note{If \code{m0} and / or \code{C0} are provided and \code{constrain} is \code{TRUE}, no check is done to ensure that the provided values are consistent with the specified constraints. This should not be too problematic in most cases, but users should be aware of it.}

\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{EnsSLLT.modeldef}}, \code{\link{SLLTSmooth}}, \code{\link{EnsEBMtrendSmooth}}, \code{\link{dlm.SafeMLE}}, \code{\link{summary.dlmMLE}}}

\examples{
## Load package data 
data(SSP585data)
GCMIDs <- 
  unlist(lapply(strsplit(names(GlobalTemps585)[-(1:2)], "\\\\."), 
                FUN=function(x) x[1]))

## Prior means and standard deviations for log variances
priors <- rbind(c(1, 1), # alpha    
                c(-3.7, 3.45),    # log(sigsq.0)
                c(-8.3, 3.45),    # log(tausq.0)
                c(-8.3, 3.45),    # log(tausq.w)
                c(-3.7, 3.45),    # log(sigsq.1)
                c(-8.3, 3.45))    # log(tausq.1)

## And prior mean and variance for state vector in 1949
m0 <- rep(0, 2*((length(unique(GCMIDs))+1)))
m0[1] <- 14  # Prior expectation of global mean temperature in 1949
k <- 25

## Fit and smooth the observations (this could take a little time)
EnsSLLTFit <- EnsSLLTSmooth(as.matrix(GlobalTemps585[,-1]), 
                            Groups=as.numeric(as.factor(GCMIDs)), 
                            m0=m0, kappa=k, prior.pars=priors,
                            messages=FALSE)
                            
## Examine fitted mimic, and plot
summary(EnsSLLTFit$Theta)
SmoothPlot(GlobalTemps585, EnsSLLTFit, DatColours=c("black", "coral3"),
           Groups=as.numeric(as.factor(GCMIDs)), 
           PredColours=c("darkblue", "darkgoldenrod"), PlotConsensus=TRUE, 
           EnsTransp=0.2)
}