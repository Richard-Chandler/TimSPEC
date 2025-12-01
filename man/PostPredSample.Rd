\name{PostPredSample}
\alias{PostPredSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample from posterior predictive distribution of observable time series}
\description{To draw samples from the posterior predictive distribution of observable time series (real-world or ensemble members) using a fitted dynamic linear model (DLM).}
\usage{
PostPredSample(Data, ModelBundle, Build, N, Random = TRUE, 
               Quantile = TRUE, Level = ifelse(Quantile, 1.96, 0.99), 
               df = NULL, Antithetic = c("Mean", "Scale"), 
               Importance = FALSE, CheckMax = FALSE, 
               ReplaceOnFail = Random, PlotFails = FALSE, WhichEls = 1, 
               NonNeg = FALSE, ReplaceAll = FALSE, debug = FALSE, 
			   messages = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{The data: a matrix in which the first column contains an observed time series and the remainder contain ensemble members.}
  \item{ModelBundle}{A list containing at least the named components \code{Model} (an object of class \code{dlm}) and \code{Theta} (a list containing a maximimum likelihood or MAP fit of a state space model, typically the result of a call to \code{\link{dlm.SafeMLE}} or equivalent). \strong{NB} it's assumed that the model is set up in such a way that the "real-world" time series is the first element of each observation vector. The \code{*Smooth} routines in the \pkg{TimSPEC} package produce objects with the required structure. }
  \item{Build}{A function that will construct the \code{dlm} object representing the model being fitted. Typically one of the \code{.modeldef} functions in the \pkg{TimSPEC} package. }
  \item{N}{Number of posterior samples to draw. Ignored if \code{Random} is \code{FALSE} (see below). }
  \item{Random}{If \code{TRUE} (the default) then samples are drawn randomly from the (approximate) posterior. If \code{FALSE}, a systematic sample is taken. See \code{\link{dlm.ThetaSample}} for more details. }
  \item{Quantile}{Used only if \code{Random} is \code{FALSE}. See \code{\link{dlm.ThetaSample}}. }
  \item{Level}{See \code{\link{dlm.ThetaSample}}. }
  \item{df}{If non-\code{NULL} and if \code{Random} is \code{TRUE}, sampling is from a scaled and shifted multivariate \eqn{t}-distribution with \code{df} degrees of freedom. Otherwise a multivariate Gaussian is used. }
  \item{Antithetic}{Controls the use of antithetic variables to reduce sampling variation if \code{Random} is \code{TRUE}. The argument is ignored if \code{Random} is \code{FALSE}. See \code{\link{dlm.ThetaSample}} for full details. }
  \item{Importance}{Logical scalar, indicating whether or not to use importance sampling. }
  \item{CheckMax}{Logical scalar, indicating whether to check the results from the importance sampling calculations for any parameter sets with a higher (penalised) log-likelihood than that in \code{ModelBundle$Theta}. Used only if \code{Importance} is \code{TRUE}.}
  \item{ReplaceOnFail}{Logical scalar, used to control the behaviour of the routine if the posterior sampling of states fails for some parameter sets. If \code{TRUE}, these sets are replaced with alternative samples; otherwise they are retained and the corresponding samples of states and observable are set to \code{NA}. }
  \item{PlotFails}{If TRUE, produce a plot on the current graphics device (using \code{\link{ExamineFailures}}) to examine any parameter sets for which the posterior sampling of states failed. This will generate a reassuring warning message if there are no failures.}
  \item{WhichEls}{Numeric vector choosing the series (i.e. columns of \code{Data}) for which sampling is required. The default value of 1 selects just the real-world series. }
  \item{NonNeg}{If \code{TRUE}, the sampled "observations" will be drawn from gamma distributions with mean and variance matching those from the posterior predictive distribution --- this ensures that the results are all strictly positive. If \code{FALSE} (the default), they will be drawn from normal distributions. }
  \item{ReplaceAll}{If \code{TRUE}, the returned object will contain random samples for all time points. If \code{FALSE} (the default), only the values missing from \code{Data} are simulated (and the non-missing values are included in the output). }
  \item{debug}{As in \code{\link[dlm]{dlmLL}}. }
  \item{messages}{Logical scalar, controlling whether to print progress to screen.}
  \item{\dots}{Additional arguments to \code{build}. }
}

\value{The function is intended to return a list containing the following elements:

\item{Thetas}{The result of the call to \code{\link{dlm.ThetaSample}};}
\item{States}{The result of the call to \code{\link{SampleStates}};} 
\item{Obs}{The result of the call to \code{\link{SampleObs}};} 
\item{Model}{The \code{dlm} object used to generate the samples;}
\item{Weights}{Information on importance weights for each sample, if \code{Importance} is \code{TRUE}; otherwise \code{NULL}. If non-\code{NULL}, this is the result of calling \code{\link{dlm.ImportanceWts}}.}

If, however, \code{Importance} and \code{CheckMax} are both \code{TRUE} and the routine finds a sampled parameter set with a larger log-posterior / (penalised) log-likelihood than that contained in \code{ModelBundle$Theta}, then the only components of the result that are guaranteed to exist are \code{Thetas} and \code{Weights}. In this case, the \code{Thetas} component will have a \code{BestID} attribute giving the index of the parameter set with the highest log-posterior / (penalised) log-likelihood. This parameter set can then be extracted for use as the \code{theta} argument in a new call to \code{\link{dlm.SafeMLE}}, in an attempt to locate a global maximum. The examples below illustrate the use of the \code{BestID} attribute to check for problems.  
}

\details{The function is essentially a wrapper to \code{\link{dlm.ThetaSample}}, \code{\link{SampleStates}} and \code{\link{SampleObs}}, with the exception that if \code{ReplaceOnFail} is \code{TRUE} then the \code{FailedIDs} attribute of the \code{States} component of the result corresponds to the samples that have been replaced after failing initially.}

\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

\seealso{\code{\link{dlm.ThetaSample}}, \code{\link{SampleStates}}, \code{\link{SampleObs}}, \code{\link{SLLTSmooth}} (for details of example below)}

\examples{
## Start by running the example for EBMtrendSMooth to create a smooth object
example(EBMtrendSmooth)

## Now draw approximate posterior predictive samples from the created object
## (which is named EBMFit) using Laplace approximation (i.e. not using 
## importance sampling). Set seed for reproducible results.
set.seed(2000)
PPSApprox <- 
    PostPredSample(as.matrix(GlobalTemps585$Obs), ModelBundle=EBMFit, 
                 Build=EBMtrend.modeldef, N=250, m0=m0, kappa=k,
                 Xt=ERF585$NetERF)

## Plot the result. NB no ensemble members are used here, hence NAs in colours
SmoothPlot(GlobalTemps585[,1:2], Samples=PPSApprox, 
           DatColours=c("black", NA), PredColours=c("darkblue", NA), 
           SmoothPIs = TRUE, ylim=c(13, 19), 
           main="Observation-only EBM trend: Laplace approx to posterior")

## Same again with importance sampling
PPS <- 
    PostPredSample(as.matrix(GlobalTemps585$Obs), ModelBundle=EBMFit, 
                 Importance=TRUE, CheckMax=TRUE, Build=EBMtrend.modeldef, 
				 N=250, m0=m0, kappa=k, Xt=ERF585$NetERF)
				 
##  Using BestID attribute to check that there are no better parameter sets
##  than the one in EBMFit
BestID <- attr(PPS$Thetas, "BestID")
cat("Check for sampled parameter sets with higher log-posterior than EBMFit ... ")
if(is.null(BestID)) cat("None found") else {
  cat("\nSet ",BestID,"has a higher log-posterior - refitting ...\n")
  EBMFit <- EBMtrendSmooth(GlobalTemps585$Obs, Xt=ERF585$NetERF,
                         kappa=k, UsePhi=TRUE, prior.pars=priors,
                         theta=PPS$Thetas[BestID,])
  
  PPS <- 
    PostPredSample(as.matrix(GlobalTemps585$Obs), ModelBundle=EBMFit, 
                 Importance=TRUE, CheckMax=TRUE, Build=EBMtrend.modeldef, 
				 N=250, m0=m0, kappa=k, Xt=ERF585$NetERF)
}

SmoothPlot(GlobalTemps585[,1:2], Samples=PPS, 
           DatColours=c("black", NA), PredColours=c("darkblue", NA), 
           SmoothPIs = TRUE, ylim=c(13, 19), 
           main="Observation-only EBM trend: posterior with importance sampling")

## For completeness: check importance sampling weights
CumulativeWeightPlot(PPS$Weights[,3])

## Looks OK: distribution isn't uniform but weights are distributed across
## all samples. Could also visualise how the sampled parameter values are
## reweighted:
PlotImportanceWts(PPS$Thetas, PPS$Weights, oma=c(10,3,5,3))

}