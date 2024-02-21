\name{EnsSLLT.modeldef}
\alias{EnsSLLT.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postprocessing an exchangeable ensemble of time series using smooth local linear trends}
\description{Sets up the structure of a simple "Ensemble smooth local linear trend" state space time series model for an observed time series and exchangeable ensemble.}
\usage{
EnsSLLT.modeldef(theta, m0 = NULL, C0=NULL, kappa = 1e+06, NRuns, 
                 Groups = NULL, discrepancy = "varying", 
                 UseAlpha = TRUE, constrain = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{Parameter vector. See "Details".}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{Initialisation variance for diffuse elements of the state vector. Ignored if \code{C0} is non-\code{NULL}.}
  \item{NRuns}{Number of ensemble members}
  \item{Groups}{Vector indicating group membership for each ensemble member. This should be used if an ensemble contains multiple runs from one or more simulators (climate models): in this case the elements of \code{Groups} should be integers between 1 and \code{G} where \code{G} is the total number of simulators contributing to the ensemble: \code{Groups[i]} is the index of the simulator that was used to produce the \eqn{i}th ensemble member. If \code{NULL} (the default), each member is assumed to come from a different simulator. 
}
  \item{discrepancy}{Either "varying" (the default) or "constant". If "constant" then the discrepancies between each ensemble member and reality are constant through time; otherwise they evolve as random walks.}
  \item{UseAlpha}{Logical scalar indicating whether or not to include a scaling factor (\eqn{\alpha}{alpha} in the "Details" section of \code{\link{EnsEBMtrend.modeldef}}) when relating the ensemble consensus trend to that in the real climate system.  If \code{FALSE}, the value of \eqn{\alpha}{alpha} is set to 1.}
  \item{constrain}{Logical scalar indicating whether or not to constrain the discrepancy trends and slopes for the individual ensemble members to sum to zero. Without this constraint, the model is formally not identifiable --- although linear combinations of elements of the state vector are.}
}
\details{The function defines dynamic linear models for exchangeable ensembles of time series with a single value per year. The structure of each individual series is as defined by \code{\link{SLLT.modeldef}}; but there are shared discrepancies between the trends for each of the ensemble members and that for the real climate system. Specifically, the trend for the \eqn{i}th ensemble member is represented as
\deqn{\mu_{i}(t) = \alpha\mu_{0}(t) + \tilde{\mu}_{\omega}(t) + \tilde{\mu}_{i}(t)\ ,}{mu_i(t) = alpha*mu_0(t) + mu.tilde_{omega}(t) + mu.tilde_{i}(t) ,}
where:

\itemize{
\item \eqn{\mu_{0}(\cdot)}{mu_{0}(.)} is the trend for the real climate system, represented as \eqn{\mu_{0}(t) = \mu_{0}(t-1) + \beta_{0}(t-1)}{mu_0[t] = mu_0[t-1] + beta_0[t-1]} as in \code{\link{SLLT.modeldef}}; 
\item \eqn{\tilde{\mu}_{\omega}(\cdot)}{mu.tilde_{w}(.)} is the shared discrepancy between the trend for the ensemble members and the real climate, which itself follows a SLLT process as \eqn{\tilde{\mu}_{\omega}(t) = \tilde{\mu}_{\omega}(t-1) + \tilde{\beta}_{\omega}(t-1)}{mu.tilde_w[t] = mu.tilde_w[t-1] + beta.tilde_w[t-1]}; 
\item \eqn{\tilde{\mu}_{i}(\cdot)}{mu.tilde_{}(.)} is the additional discrepancy between the trend for the \eqn{i}th ensemble member and the ensemble consensus \eqn{\alpha\mu_{0}(\cdot) + \tilde{\mu}_{\omega}(\cdot)}{alpha*mu_0(.) + mu.tilde_{omega}(.)}. This discrepancy process evolves as \eqn{\tilde{\mu}_{i}(t) = \tilde{\mu}_{i}(t-1) + \tilde{\beta}_{i}(t-1)}{mu.tilde_i[t] = mu.tilde_i[t-1] + beta.tilde_i[t-1]}.}

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

When \code{constrain} is \code{FALSE}, the state vector for the model at time \eqn{t} is defined as 

\deqn{\mathbf{S}(t) = \left(\mu_{0}(t)\ \beta_{0}(t)\ \tilde{\mu}_{\omega}(t)\ \tilde{\beta}_{\omega}(t)\ \tilde{\mu}_{1}(t)\ \tilde{\beta}_{1}(t)\ \cdots\ \tilde{\mu}_{G}(t)\ \tilde{\beta}_{G}(t)\right)'.}{S[t] = (mu_0[t] beta_0[t] mu.tilde_w[t] beta.tilde_w[t] mu.tilde_1[t] beta.tilde_1[t] ... mu.tilde_G[t] beta.tilde_G[t])',}

where \eqn{G} is the number of simulators contributing to the ensemble (see documentation for the \code{Groups} argument above). When \code{constrain} is \code{TRUE} (the default), the elements corresponding to simulator \eqn{G} are removed: in this case, they are defined in terms of the other elements by the sum-to-zero constraint. 


}
\value{A list defining a dynamical linear model, in the format required by routines in the \pkg{dlm} library.}

\note{If \code{m0} and / or \code{C0} are provided and \code{constrain} is \code{TRUE}, no check is done to ensure that the provided values are consistent with the specified constraints. This should not be too problematic in most cases, but users should be aware of it.}

\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}

\seealso{\code{\link{EnsSLLTSmooth}} for a convenience wrapper function; \code{\link{SLLT.modeldef}}; \code{\link{EnsEBMtrend.modeldef}}, for a physically motivated alternative for application to climate time series.}
