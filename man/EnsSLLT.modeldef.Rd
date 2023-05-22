\name{EnsSLLT.modeldef}
\alias{EnsSLLT.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postprocessing an exchangeable ensemble of time series using smooth local linear trends}
\description{ets up the structure of a simple "Ensemble smooth local linear trend" state space time series model for an observed time series and exchangeable ensemble.}
\usage{
EnsSLLT.modeldef(theta, m0 = NULL, kappa = 1e+06, NRuns, Groups = NULL, discrepancy = "varying", UseAlpha = TRUE, constrain = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{Parameter vector. See "Details".}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{kappa}{Initialisation variance for diffuse elements of the state vector.}
  \item{NRuns}{Number of ensemble members}
  \item{Groups}{Vector indicating group membership for each ensemble member. This should be used if an ensemble contains multiple runs from one or more simulators (climate models): in this case the elements of \code{Groups} should be integers between 1 and \code{G} where \code{G} is the total number of simulators contributing to the ensemble: \code{Groups[i]} is the index of the simulator that was used to produce the \eqn{i}th ensemble member. If \code{NULL} (the default), each member is assumed to come from a different simulator. 
}
  \item{discrepancy}{Either "varying" (the default) or "constant". If "constant" then the discrepancies between each ensemble member and reality are constant through time; otherwise they evolve as random walks.}
  \item{UseAlpha}{Logical scalar indicating whether or not to include a scaling factor (\eqn{\alpha}{alpha} in the "Details" section of \code{\link{EnsEBMtrend.modeldef}}) when relating the ensemble consensus trend to that in the real climate system.  If \code{FALSE}, the value of \eqn{\alpha}{alpha} is set to 1.}
  \item{constrain}{Logical scalar indicating whether or not to constrain the discrepancy trends and slopes for the individual ensemble members to sum to zero. Without this constraint, the model is formally not identifiable --- although linear combinations of elements of the state vector are.}
}
\details{The function defines dynamic linear models for exchangeable ensembles of time series with a single value per year. The structure of each individual series is as defined by \code{\link{SLLT.modeldef}}; but there are shared discrepancies between the trends for each of the ensemble members and that for the real climate system. The discrepancy structure is the same as that described for \code{\link{EnsEBMtrend.modeldef}}, with appropriate modifications to account for the different trend dynamics of the SLLT model.

If \code{discrepancy="varying"} and \code{UseAlpha} is \code{TRUE} then the parameter vector \code{theta} is defined as follows:

\describe{
\item{\code{theta[1]}}{This is \eqn{\alpha}{alpha}, the scaling on the relationship between the ensemble consensus and observed trend processes.}
\item{\code{theta[2]}}{This is \eqn{\log(\sigma^{2}_{0})}{log(sigsq[0])}, where \eqn{\sigma^{2}_{0}}{sigsq[0])} is the "measurement error" variance in the real-world series.}
\item{\code{theta[3]}}{This is \eqn{\log(\tau^{2}_{0})}{log(tausq[0])}, where \eqn{\tau^{2}_{0}}{tausq[0]} is slope innovation variance in the real-world series (see documentation for \code{\link{SLLT.modeldef}}).}
\item{\code{theta[4]}}{This is \eqn{\log(\tau^{2}_{\omega})}{log(tausq[omega])}, where \eqn{\tau^{2}_{\omega}}{tausq[omega]} is the innovation variance for the random-walk process representing the shared discrepancy between the trend slopes in the ensemble members and the real world.}
\item{\code{theta[5]}}{This is \eqn{\log(\sigma^{2}_{1})}{log(sigsq[1])}, where \eqn{\sigma^{2}_{1}}{sigsq[1]} is the common "measurement error" variance for each of the ensemble members.}
\item{\code{theta[6]}}{This is \eqn{\log(\tau^{2}_{1})}{log(tausq[1])}, where \eqn{\tau^{2}_{1}}{tausq[1]} is the common slope innovation variance for each of the ensemble members' individual trend discrepancies about their consensus trend.}
}

If \code{discrepancy="constant"} and \code{UseAlpha} is \code{FALSE} then \code{theta} is a 4-element vector containing \eqn{\alpha}{alpha}, \eqn{\log(\sigma^{2}_{0})}{log(sigsq[0])}, \eqn{\log(\tau^{2}_{0})}{log(tausq[0])} and \eqn{\log(\sigma^{2}_{1})}{log(sigsq[1])}: in this case, \eqn{\tau^{2}_{\omega}}{tausq[omega]} and \eqn{\tau^{2}_{1}}{tausq[1]} are both taken as zero.

If \code{UseAlpha} is \code{FALSE} then \eqn{\alpha}{alpha} is fixed to 1 and omitted from \code{theta} --- so that in this case \code{theta[1]} is \eqn{\log(\sigma^{2}_{0})}{log(sigsq[0])}, \code{theta[2]} is \eqn{\log(\tau^{2}_{0})}{log(tausq[0])} etc. This argument makes no difference to the results if \code{discrepancy="constant"}. 
}
\value{A list defining a dynamical linear model, in the format required by routines in the \pkg{dlm} library.}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{SLLT.modeldef}}, \code{\link{EnsEBMtrend.modeldef}}}
