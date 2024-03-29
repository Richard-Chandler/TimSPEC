\name{PostPredSample}
\alias{PostPredSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample from posterior predictive distribution of observable time series}
\description{To draw samples from the posterior predictive distribution of observable time series (real-world or ensemble members) using a fitted dynamic linear model (DLM).}
\usage{
PostPredSample(Data, ModelBundle, Build, N, Random = TRUE, 
               Quantile = TRUE, Level = ifelse(Quantile, 1.96, 0.99), 
               df = NULL, Antithetic = c("Mean", "Scale"), 
               Importance = FALSE, ReplaceOnFail = Random, 
               PlotFails = TRUE, WhichEls = 1, NonNeg = FALSE, 
               ReplaceAll = FALSE, debug = FALSE, messages = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{The data: a matrix in which the first column contains an observed time series and the remainder contain ensemble members. }
  \item{ModelBundle}{A list containing at least the named components \code{Model} (an object of class \code{dlm}) and \code{Smooth} (the result of Kalman Smoothing a dataset using the DLM stored in \code{Model}, via the \code{\link{dlmSmooth}} command). \strong{NB} it's assumed that the model is set up in such a way that the "real-world" time series is the first element of each observation vector. The \code{*Smooth} routines in the \pkg{PECTSP} package produce objects with the required structure. }
  \item{Build}{A function that will construct the \code{dlm} object representing the model being fitted. Typically one of the \code{.modeldef} functions in the \pkg{PECTSP} package. }
  \item{N}{Number of posterior samples to draw. Ignored if \code{Random} is \code{FALSE} (see below). }
  \item{Random}{If \code{TRUE} (the default) then samples are drawn randomly from the (approximate) posterior. If \code{FALSE}, a systematic sample is taken. See \code{\link{dlm.ThetaSample}} for more details. }
  \item{Quantile}{Used only if \code{Random} is \code{FALSE}. See \code{\link{dlm.ThetaSample}}. }
  \item{Level}{See \code{\link{dlm.ThetaSample}}. }
  \item{df}{If non-\code{NULL} and if \code{Random} is \code{TRUE}, sampling is from a scaled and shifted multivariate \eqn{t}-distribution with \code{df} degrees of freedom. Otherwise a multivariate Gaussian is used. }
  \item{Antithetic}{Controls the use of antithetic variables to reduce sampling variation if \code{Random} is \code{TRUE}. The argument is ignored if \code{Random} is \code{FALSE}. See \code{\link{dlm.ThetaSample}} for full details. }
  \item{Importance}{Logical scalar, indicating whether or not to use importance sampling. }
  \item{ReplaceOnFail}{Logical scalar, used to control the behaviour of the routine if the posterior sampling of states fails for some parameter sets. If \code{TRUE}, these sets are replaced with alternative samples; otherwise they are retained and the corresponding samples of states and observable are set to \code{NA}. }
  \item{PlotFails}{If TRUE and if the posterior sampling of states failed for any parameter sets, produce a plot on the current graphics device (using \code{\link{ExamineFailures}}) to examine the failed sets. }
  \item{WhichEls}{Numeric vector choosing the series (i.e. columns of \code{Data}) for which sampling is required. The default value of 1 selects just the real-world series. }
  \item{NonNeg}{If \code{TRUE}, the sampled "observations" will be drawn from gamma distributions with mean and variance matching those from the posterior predictive distribution --- this ensures that the results are all strictly positive. If \code{FALSE} (the default), they will be drawn from normal distributions. }
  \item{ReplaceAll}{If \code{TRUE}, the returned object will contain random samples for all time points. If \code{FALSE} (the default), only the values missing from \code{Data} are simulated (and the non-missing values are included in the output). }
  \item{debug}{As in \code{\link{dlmLL}}. }
  \item{messages}{Logical scalar, controlling whether to print progress to screen.}
  \item{\dots}{Additional arguments to \code{build}. }
}
\details{The function is essentially a wrapper to \code{\link{dlm.ThetaSample}}, \code{\link{SampleStates}} and \code{\link{SampleObs}}. }
\value{The function returns a list containing the following elements:

\item{Thetas}{The result of the call to \code{\link{dlm.ThetaSample}};}
\item{States}{The result of the call to \code{\link{SampleStates}};} 
\item{Obs}{The result of the call to \code{\link{SampleObs}};} 
\item{Model}{The \code{dlm} object used to generate the samples;}
\item{Weights}{Information on importance weights for each sample, if \code{Importance} is \code{TRUE}; otherwise \code{NULL}. If non-\code{NULL}, this is the result of calling \code{\link{dlm.ImportanceWts}};}  \item{Replacements}{A logical scalar indicating whether some of the initial samples were replaced --- see \code{ReplaceOnFail} argument above.} 
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

\seealso{\code{\link{dlm.ThetaSample}}, \code{\link{SampleStates}}, \code{\link{SampleObs}}, \code{\link{SLLTSmooth}} (for details of example below)}

\examples{
## Start by running the example for EBMtrendSMooth to create a smooth object
example(EBMtrendSmooth)

## Now draw approximate posterior predictive samples from the created object
## (which is named EBMFit) using Laplace approximation (i.e. not using 
## importance sampling)
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
                 Importance=TRUE, Build=EBMtrend.modeldef, N=250, 
                 m0=m0, kappa=k, Xt=ERF585$NetERF)
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