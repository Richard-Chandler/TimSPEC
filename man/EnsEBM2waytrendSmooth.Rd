\name{EnsEBM2waytrendSmooth}
\alias{EnsEBM2waytrendSmooth}
\alias{EnsEBM2waytrend.NegLL}
\title{Postprocessing an two-way structured ensemble of time series using EBM-inspired dynamic linear models}
\description{Fits a model of the form defined by \code{\link{EnsEBM2waytrend.modeldef}}, then applies the Kalman Smoother to obtain an estimate of the underlying trends with prediction intervals. The motivating application is for a regional ensemble in which each member has been produced by running one of a collection of regional climate models (RCMs), with boundary conditions generated using one of a collection of global climate models (GCMs): the documentation below is framed in the context of such an application. }
\usage{
EnsEBM2waytrendSmooth(Y, Xt, m0 = NULL, kappa = 1e+06, Groups, 
                      prior.pars = NULL, interactions = "none", 
                      UseAlpha = TRUE, UsePhi = TRUE, constrain = TRUE,
                      messages = TRUE, Use.dlm = FALSE, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Matrix containing an observed time series in its first column and series from a regional ensemble (i.e. an ensemble in which each series is produced by one of \eqn{R} RCMS with boundary conditions provided by one of \eqn{G} GCMs) in the remaining columns.}
  \item{Xt}{Vector containing the corresponding forcings for each row of Y}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{kappa}{Variance used to initialise diffuse elements of the state vector.}
  \item{Groups}{Two-column matrix such that \code{Groups[i,1]} is the number of the RCM used to produce the \eqn{i}th ensemble member (i.e. column \code{i+1} of \code{Y}) and \code{Groups[i,2]} is the number of the GCM.}
  \item{prior.pars}{Optional 2-column matrix containing the means and standard deviations of independent Gaussian priors for the transformed model parameters (log variances, and logit thermal inertias). If this is provided, maximum \emph{a posteriori} (or penalised maximum likelihood) estimation is used; otherwise just standard maximum likelihood.}
  \item{interactions}{Character scalar, indicating which RCM:GCM "interaction" terms to include in the state vector. Options are \code{"none"} to exclude them completely, \code{"all"} to include all of them and \code{"available"} to include just those for which thecorresponding RCM:GCM combinations appear in \code{Groups}. Default is \code{"none"}, because the interaction terms substantially increase the dimension and are therefore likely to slow things down considerably.}
  \item{UseAlpha}{Logical scalar indicating whether or not to include a scaling factor (\eqn{\alpha}{alpha} in the "Details" section of \code{\link{EnsEBMtrend.modeldef}}) when relating the ensemble consensus trend to that in the real climate system.  If \code{FALSE}, the value of \eqn{\alpha}{alpha} is set to 1.}
  \item{UsePhi}{Controls whether or not to include thermal inertia parameters \eqn{\phi_{0}}{phi[0]} and \eqn{\phi_{1}}{phi[1]} in the model structure. If \code{FALSE}, both parameters are set to zero.}
  \item{constrain}{Logical scalar, controlling whether or not to ensure identifiability by imposing sum-to-zero constraints on relevant RCM-specific and GCM-specific effects.}
  \item{messages}{Logical scalar controlling whether or not to print progress to screen.}
  \item{Use.dlm}{Logical scalar, used to control the method used to fit the model: see \code{\link{dlm.SafeMLE}}.}
  \item{debug}{As in \code{\link{dlmMLE}}.}
}
\details{At present, the memory leak in the \code{C} code for \code{\link{dlmLL}} creates problems when calculating the log-likelihood for this DLM using \code{\link{dlm.SafeLL}}; in turn, this prevents the use of \code{\link{dlm.SafeMLE}} for this DLM. To get around this, function \code{EnsEBM2waytrend.NegLL} provides the required log-likelihood calculations. This does mean, however, that \code{Use.dlm=TRUE} is currently unsupported. 
}
\value{A named list with three components:
\item{Theta}{Object containing parameter estimation results; the structure is the same as that arising from a call to \code{\link{dlmMLE}} or \code{\link{dlm.SafeMLE}}}
\item{Model}{An object of class \code{dlm} containing the fitted model itself}
\item{Smooth}{An object containing the result of Kalman Smoothing the input series using the fitted model. This is produced using \code{\link{dlmSmooth}}}
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{EnsEBM2waytrend.modeldef}}, \code{\link{dlm.SafeMLE}}, \code{\link{summary.dlmMLE}}}

