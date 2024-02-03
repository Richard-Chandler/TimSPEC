\name{EnsEBMtrend.modeldef}
\alias{EnsEBMtrend.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{State space model for an exchangeable ensemble of time series, inspired by Energy Balance Models}
\description{Sets up the structure of a state space model for an observed time series and exchangeable ensemble, with trend formulation inspired by the form of a simple energy balance model (EBM).}
\usage{
EnsEBMtrend.modeldef(theta, Xt, m0 = NULL, C0=NULL, kappa, NRuns, 
                     Groups = NULL, UseAlpha = TRUE, UsePhi = TRUE, 
                     constrain = TRUE)
}
\arguments{
  \item{theta}{Parameter vector. See "Details".}
  \item{Xt}{Vector containing time series of effective forcings.}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{"Large" value with which to initialise variances of diffuse components of state vector at time 0. Ignored if \code{C0} is non-\code{NULL}.}
  \item{NRuns}{Number of ensemble members}
  \item{Groups}{Vector of length \code{NRuns}, indicating group membership. This should be used if an ensemble contains multiple runs from one or more simulators (i.e. climate models): in this case the elements of \code{Groups} should be integers between \code{1} and \code{G} where \code{G} is the total number of simulators contributing to the ensemble: \code{Groups[i]} is the index of the simulator that was used to produce the \eqn{i}th ensemble member. If \code{NULL}, each member is assumed to come from a different simulator.}
  \item{UseAlpha}{See "Details" section}
  \item{UsePhi}{Controls whether or not to include thermal inertia parameters in the model structure.}
  \item{constrain}{Logical scalar indicating whether or not to constrain the discrepancy trends and slopes for the individual ensemble members to sum to zero. Without this constraint, the model is formally not identifiable --- although linear combinations of elements of the state vector are.}
}
\details{The function defines EBM-inspired dynamic linear models for exchangeable ensembles of time series with a single value per year. The structure of each individual series is as defined by \code{\link{EBMtrend.modeldef}}; but there are shared discrepancies between the trends for each of the ensemble members and that for the real climate system. Specifically, the trend for the \eqn{i}th ensemble member is represented as
\deqn{\mu_{i}(t) = \alpha\mu_{0}(t) + \tilde{\mu}_{\omega}(t) + \tilde{\mu}_{i}(t)\ ,}{mu_i(t) = alpha*mu_0(t) + mu.tilde_{omega}(t) + mu.tilde_{i}(t) ,}
where:

\itemize{
\item \eqn{\tilde{\mu}_{\omega}(\cdot)}{mu.tilde_{omega}(.)} is the shared discrepancy between the trend for the ensemble members and the real climate; 
\item \eqn{\tilde{\mu}_{i}(\cdot)}{mu.tilde_{}(.)} is the additional discrepancy between the trend for the \eqn{i}th ensemble member and the ensemble consensus \eqn{\alpha\mu_{0}(\cdot) + \tilde{\mu}_{\omega}(\cdot)}{alpha*mu_0(.) + mu.tilde_{omega}(.)}. }

If \code{UseAlpha} and \code{UsePhi} are both \code{TRUE} then the parameter vector \code{theta} is defined as follows:
\describe{
\item{\code{theta[1]}}{This is \eqn{\alpha}{alpha}, the scaling on the relationship between the ensemble consensus and observed trend processes.}
\item{\code{theta[2]}}{This is \eqn{\log(\sigma^{2}_{0})}{log(sigsq[0])}, where \eqn{\sigma^{2}_{0}}{sigsq[0])} is the "measurement error" variance in the real-world series.}
\item{\code{theta[3]}}{This is \eqn{\log(\tau^{2}_{0})}{log(tausq[0])}, where \eqn{\tau^{2}_{0}}{tausq[0]} is the innovation variance for the random walk drift approximation in the real-world series (see documentation for \code{\link{EBMtrend.modeldef}}).}
\item{\code{theta[4]}}{This is \eqn{\log(\tau^{2}_{\omega})}{log(tausq[omega])}, where \eqn{\tau^{2}_{\omega}}{tausq[omega]} is the innovation variance for the random walk drift approximation in the shared trend discrepancy process \eqn{\tilde{\mu}_{\omega}(\cdot)}{mu.tilde_{omega}(.)}.}
\item{\code{theta[5]}}{This is \eqn{\log(\sigma^{2}_{1})}{log(sigsq[1])}, where \eqn{\sigma^{2}_{1}}{sigsq[1]} is the common "measurement error" variance for each of the ensemble members.}
\item{\code{theta[6]}}{This is \eqn{\log(\tau^{2}_{1})}{log(tausq[1])}, where \eqn{\tau^{2}_{1}}{tausq[1]} is the common drift innovation variance for each of the ensemble members' individual trend discrepancies about their consensus trend.}
\item{\code{theta[7]}}{This is logit(\eqn{phi_{0}}{phi[0]}), where \eqn{phi_{0}}{\phi[0]} is the temporal dependence parameter (inertia coefficient) in the mean component of the real-world trend.}
\item{\code{theta[8]}}{This is logit(\eqn{phi_{1}}{phi[1]}), where \eqn{phi_{1}}{\phi[1]} is the common temporal dependence parameter in the mean component of the trends in all ensemble members.}
}

If \code{UseAlpha} is \code{FALSE} then the value of \eqn{\alpha}{alpha} is fixed at 1 and is omitted from \code{theta} so that \code{theta[1]} is \eqn{\log(\sigma^{2}_{0})}{log(sigsq[0]}), \code{theta[2]} is \eqn{\log(\tau^{2}_{0})}{log(tausq[0])} etc. 

If \code{UsePhi} is \code{FALSE}, then both of the final elements of \code{theta} are omitted and the values of both \eqn{\phi_{0}}{phi[0]} and \eqn{\phi_{1}}{phi[1]} are set to 0. 
}

\value{A list defining a dynamical linear model, in the format required by routines in the \pkg{dlm} library.}
\references{
REFERENCES TO BE ADDED
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{EBMtrend.modeldef}}}
