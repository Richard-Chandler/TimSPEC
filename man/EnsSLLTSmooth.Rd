\name{EnsSLLTSmooth}
\alias{EnsSLLTSmooth}
\alias{EnsEBMtrendSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postprocessing an exchangeable ensemble of time series using dynamic linear models}
\description{Fit models of the form defined by \code{\link{EnsSLLT.modeldef}} or \code{\link{EnsEBMtrend.modeldef}}, then applies the Kalman Smoother to obtain an estimate of the underlying trends with prediction intervals.}
\usage{
EnsSLLTSmooth(Y, m0 = NULL, C0=NULL, kappa = 1e+06, discrepancy = "varying",
              Groups = NULL, compact=!is.null(Groups), UseAlpha = TRUE,
              prior.pars = NULL, theta=NULL, constrain = TRUE, 
              Tiny = 1/kappa, ObsSmooth, Ens0Theta, messages = TRUE, 
              Use.dlm = FALSE, debug = FALSE, ...)

EnsEBMtrendSmooth(Y, Xt, Groups = NULL, m0 = NULL, C0=NULL, kappa = 1e+06, 
                  compact=!is.null(Groups), prior.pars = NULL, theta=NULL,
                  UseAlpha = TRUE, UsePhi = TRUE, constrain = TRUE, 
                  messages = TRUE, Use.dlm = FALSE, debug = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Matrix containing an observed time series in its first column and exchangeable series from an ensemble in the remaining columns}
  \item{Xt}{Vector containing the corresponding forcings for each row of Y}
  \item{Groups}{Vector indicating group membership for each ensemble member. This should be used if an ensemble contains multiple runs from one or more simulators (climate models): in this case the elements of \code{Groups} should be integers between 1 and \code{G} where \code{G} is the total number of simulators contributing to the ensemble: \code{Groups[i]} is the index of the simulator that was used to produce the \eqn{i}th ensemble member. If \code{NULL} (the default), each member is assumed to come from a different simulator. 
}
  \item{compact}{Logical scalar indicating whether to compact the ensemble using \code{\link{CompactEns}} before fitting. For ensembles with many members per simulator, this may speed up the fitting without sacrificing information because the simulator-specific means are sufficient statistics. Defaults to \code{TRUE} if \code{Groups} is non-\code{NULL} and \code{FALSE} otherwise (because \code{\link{CompactEns}} just returns the original ensemble in that case).}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{Variance used to initialise diffuse elements of the state vector.}
  \item{discrepancy}{Either "varying" (the default) or "constant". If "constant" then the discrepancies between each ensemble member and reality in \code{EnsSLLTSmooth} are constant through time; otherwise they evolve as random walks.}
  \item{UseAlpha}{Logical scalar indicating whether or not to include a scaling factor (\eqn{\alpha}{alpha} in the "Details" section of \code{\link{EnsEBMtrend.modeldef}}) when relating the ensemble consensus trend to that in the real climate system.  If \code{FALSE}, the value of \eqn{\alpha}{alpha} is set to 1.}
  \item{prior.pars}{Optional 2-column matrix containing the means and standard deviations of independent Gaussian priors for transformed model parameters (log variances and, for \code{EnsEBMtrendSmooth}, logit thermal inertias). If this is provided, maximum \emph{a posteriori} (or penalised maximum likelihood) estimation is used; otherwise just standard maximum likelihood. For details of the model parameters, see \code{\link{EnsSLLT.modeldef}} and \code{\link{EnsEBMtrend.modeldef}}.}
  \item{theta}{Optional initial value for parameter vector. If \code{NULL}, the routine will auto-initialise the optimisation in the first step.}
  \item{UsePhi}{Controls whether or not to include thermal inertia parameters \eqn{\phi_{0}}{phi[0]} and \eqn{\phi_{1}}{phi[1]} in the \code{\link{EnsEBMtrend.modeldef}} model structure. If \code{FALSE}, both parameters are set to zero.}
  \item{constrain}{Logical scalar controlling whether to impose sum-to-zero constraints across the ensemble members where necessary to make the model identifiable.}
  \item{Tiny}{A small positive value, used to replace negative initial estimates for variances (see \code{\link{SLLT.IniPar}})}
  \item{ObsSmooth}{Optional: result of a previous call to \code{\link{SLLTSmooth}} for the first column of \code{Y}. If present, this is used to help initialise the numerical search for an MLE if a "constant" model needs to be fitted, or if the resulting smooth is needed. Otherwise the relevant model will be fitted as part of this routine.}
  \item{Ens0Theta}{Optional parameter vector for a previous fit of an "ensemble local linear trend" model with a constant discrepancy term. If present, this will be used to help initialise the numerical search for an MLE if a "varying" model needs to be fitted here. Otherwise the relevant constant-discrepancy model will be fitted as part of this routine.}
  \item{messages}{Controls whether to print progress to screen. For details, see the help for \code{\link{dlm.SafeMLE}}.}
  \item{Use.dlm}{Logical scalar, used to control the method used to fit the model: see \code{\link{dlm.SafeMLE}}.}
  \item{debug}{As in \code{\link[dlm]{dlmMLE}}.}
  \item{...}{Other arguments to \code{\link{dlm.SafeMLE}}.}
}
\details{The function is a wrapper for the steps required to fit and apply a dynamic linear model. It starts by trying to find reasonable values for (penalised) maximum likelihood estimation of the model parameters; then carries out the estimation; and finally applies the Kalman Smoother to the data using the fitted model.

\bold{Note:} the parameter estimates, smooths and uncertainty estimates should be identical (up to the accuracy of the numerical algorithms) regardless of the value of \code{compact}. The reported log-likelihoods will potentially differ however, because of differences in the constant terms used in the definitions of multivariate normal densities in different dimensions.}
\value{A named list with four components:
\item{Data}{List containing components \code{Y} (the data used to produce the smooths, perhaps after compacting the ensemble according to the value of \code{compact}), \code{Groups} (the supplied argument, possibly modified as appropriate if \code{compact} is \code{TRUE}. The \code{Y} component has a \code{RunsPerTS} attribute if \code{compact} is \code{TRUE} - see \code{\link{CompactEns}}.}
\item{Theta}{Object containing parameter estimation results; the result of a call to \code{\link[dlm]{dlmMLE}} or \code{\link{dlm.SafeMLE}}}
\item{Model}{An object of class \code{dlm} containing the fitted model itself}
\item{Smooth}{An object containing the result of Kalman Smoothing the input series using the fitted model. This is produced using \code{\link[dlm]{dlmSmooth}}}
}

\note{If \code{m0} and / or \code{C0} are provided and \code{constrain} is \code{TRUE}, no check is done to ensure that the provided values are consistent with the specified constraints. This should not be too problematic in most cases, but users should be aware of it.}

\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{EnsSLLT.modeldef}}, \code{\link{SLLTSmooth}}, \code{\link{EnsEBMtrend.modeldef}}, \code{\link{EBMtrendSmooth}}, \code{\link{EnsEBM2waytrendSmooth}}, \code{\link{dlm.SafeMLE}}, \code{\link{summary.dlmMLE}}}

\examples{
## Load package data 
data(SSP585data)
GCMIDs <- 
  unlist(lapply(strsplit(names(GlobalTemps585)[-(1:2)], "\\\\."), 
                FUN=function(x) x[1]))

## Prior means and standard deviations for log variances (SLLT model)
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
           EnsTransp=0.2, plot.title="SLLT: smooth from compacted ensemble")

## Now do the same but without compacting the ensemble. This will
## be slower but should give the same result except for the log-
## likelihood. To speed things up, use the estimate of theta from 
## the previous fit - and notice the difference in the reported
## log-likelihoods due to the difference in dimension of the 
## observation vector. 
EnsSLLTFit2 <- EnsSLLTSmooth(as.matrix(GlobalTemps585[,-1]), 
                            Groups=as.numeric(as.factor(GCMIDs)), 
                            compact=FALSE, m0=m0, kappa=k, 
                            prior.pars=priors, theta=EnsSLLTFit$Theta$par) 

summary(EnsSLLTFit2$Theta)
SmoothPlot(GlobalTemps585, EnsSLLTFit2, DatColours=c("black", "coral3"),
           Groups=as.numeric(as.factor(GCMIDs)), 
           PredColours=c("darkblue", "darkgoldenrod"), PlotConsensus=TRUE, 
           EnsTransp=0.2, plot.title="SLLT: smooth from original ensemble")

## Prior means and standard deviations (ensemble EBM model)
priors <- rbind(c(1, 1), # alpha    
                c(-3.7, 3.45),    # log(sigsq.0)
                c(-8.3, 3.45),    # log(tausq.0)
                c(-8.3, 3.45),    # log(tausq.w)
                c(-3.7, 3.45),    # log(sigsq.1)
                c(-8.3, 3.45),    # log(tausq.1)
                c(0, 5),          # logit(phi.0)
                c(0, 5))          # logit(phi.1)

## And prior mean and variance for state vector in 1949
m0 <- rep(0, 3*((length(unique(GCMIDs))+1)))
m0[1] <- 14  # Prior expectation of global mean temperature in 1949
k <- 25

## Fit and smooth the observations (this could take a few minutes,
## so set messages to TRUE to track progress - this will also 
## print a summary of the fitted mimic at the end)
EnsEBMFit <- 
  EnsEBMtrendSmooth(as.matrix(GlobalTemps585[,-1]), Xt=ERF585$NetERF,
                              Groups=as.numeric(as.factor(GCMIDs)), 
                              m0=m0, kappa=k, prior.pars=priors,
                              messages=TRUE)
                            
## Plot the fitted mimic
SmoothPlot(GlobalTemps585, EnsEBMFit, DatColours=c("black", "coral3"),
           Groups=as.numeric(as.factor(GCMIDs)), PlotConsensus=TRUE, 
           PredColours=c("darkblue", "darkgoldenrod"), EnsTransp=0.2, 
           plot.title="EBM-based smooth from compacted ensemble")

}