\name{EnsEBMtrendSmooth}
\alias{EnsEBMtrendSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postprocessing an exchangeable ensemble of time series using EBM-inspired dynamic linear models
}
\description{Fits a model of the form defined by \code{\link{EnsEBMtrend.modeldef}}, then applies the Kalman Smoother to obtain an estimate of the underlying trends with prediction intervals.}
\usage{
EnsEBMtrendSmooth(Y, Xt, Groups = NULL, m0 = NULL, kappa = 1e+06, 
                  prior.pars = NULL, UseAlpha = TRUE, UsePhi = TRUE, 
                  constrain = TRUE, messages = TRUE, Use.dlm = FALSE, 
                  debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Matrix containing an observed time series in its first column and exchangeable series from an ensemble in the remaining columns}
  \item{Xt}{Vector containing the corresponding forcings for each row of Y}
  \item{Groups}{Vector indicating group membership for each ensemble member. This should be used if an ensemble contains multiple runs from one or more simulators (climate models): in this case the elements of \code{Groups} should be integers between 1 and \code{G} where \code{G} is the total number of simulators contributing to the ensemble: \code{Groups[i]} is the index of the simulator that was used to produce the \eqn{i}th ensemble member. If \code{NULL} (the default), each member is assumed to come from a different simulator. 
}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{kappa}{Variance used to initialise diffuse elements of the state vector.}
  \item{prior.pars}{Optional 2-column matrix containing the means and standard deviations of independent Gaussian priors for the transformed model parameters (log variances, and logit thermal inertias). If this is provided, maximum \emph{a posteriori} (or penalised maximum likelihood) estimation is used; otherwise just standard maximum likelihood.}
  \item{UseAlpha}{Logical scalar indicating whether or not to include a scaling factor (\eqn{\alpha}{alpha} in the "Details" section of \code{\link{EnsEBMtrend.modeldef}}) when relating the ensemble consensus trend to that in the real climate system.  If \code{FALSE}, the value of \eqn{\alpha}{alpha} is set to 1.}
  \item{UsePhi}{Controls whether or not to include thermal inertia parameters \eqn{\phi_{0}}{phi[0]} and \eqn{\phi_{1}}{phi[1]} in the model structure. If \code{FALSE}, both parameters are set to zero.}
  \item{constrain}{Logical scalar controlling whether to impose sum-to-zero constraints across the ensemble members where necessary to make the model identifiable.}
  \item{messages}{Logical scalar controlling whether or not to print progress to screen.}
  \item{Use.dlm}{Logical scalar, used to control the method used to fit the model: see \code{\link{dlm.SafeMLE}}.}
  \item{debug}{As in \code{\link{dlmMLE}}. Only \code{TRUE} is currently supported for this model, until the memory leak in the \pkg{dlm} package is fixed.}
}
\details{The function is a wrapper for the steps required to fit and apply a dynamic linear model. It starts by trying to find reasonable values for (penalised) maximum likelihood estimation of the model parameters; then carries out the estimation; and finally applies the Kalman Smoother to the data using the fitted model.}
\value{A named list with three components:
\item{Theta}{Object containing parameter estimation results; the result of a call to \code{\link{dlmMLE}} or \code{\link{dlm.SafeMLE}}}
\item{Model}{An object of class \code{dlm} containing the fitted model itself}
\item{Smooth}{An object containing the result of Kalman Smoothing the input series using the fitted model. This is produced using \code{\link{dlmSmooth}}}
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{EnsEBMtrend.modeldef}}}

