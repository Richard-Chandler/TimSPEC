\name{EnsEBM2waytrend.modeldef}
\alias{EnsEBM2waytrend.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{State space model for a two-way structured ensemble of time series, inspired by Energy Balance Models}
\description{Sets up the structure of a state space model for an observed time series and two-way structured ensemble, with trend formulation inspired by the form of a simple energy balance model (EBM). The motivating application is for a regional ensemble in which each member has been produced by running one of a collection of regional climate models (RCMs), with boundary conditions generated using one of a collection of global climate models (GCMs): the documentation below is framed in the context of such an application. 
}
\usage{
EnsEBM2waytrend.modeldef(theta, Xt, m0 = NULL, C0=NULL, kappa=1e6, Groups, 
                         interactions = "none", UseAlpha = TRUE, 
                         UsePhi = TRUE, constrain = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{Parameter vector. See "Details".}
  \item{Xt}{Vector containing time series of effective forcings.}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{"Large" value with which to initialise variances of diffuse components of state vector at time 0. Ignored if \code{C0} is non-\code{NULL}.}
  \item{Groups}{Two-column matrix specifying the group membership structureof the ensemble. E.g. if the ensemble contains runs from some combination of RCMs \eqn{1, 2, 3, ..., R} driven by GCMs \eqn{1, 2, 3, ..., G} then \code{Groups[i,1]} and \code{Groups[i,2]} should respectively contain the numbers of the RCM and GCM used for the \eqn{i}th ensemble member. The values of \eqn{R} and \eqn{G} will be taken from the largest numbers found in columns 1 and 2 respectively.}
  \item{interactions}{Character scalar, indicating which RCM:GCM "interaction" terms to include in the state vector. Options are \code{"none"} to exclude them completely, \code{"all"} to include all of them and \code{"available"} to include just those for which thecorresponding RCM:GCM combinations appear in \code{Groups}. Default is \code{"none"}, because the interaction terms substantially increase the dimension and are therefore likely to slow things down considerably.}
  \item{UseAlpha}{See "Details" section of help for \code{\link{EnsEBMtrend.modeldef}}}
  \item{UsePhi}{Logical scalar, controlling whether or not to include thermal inertia parameters in the model structure.}
  \item{constrain}{Logical scalar, controlling whether or not to ensure identifiability by imposing sum-to-zero constraints on relevant RCM-specific and GCM-specific effects.}
}
\details{The dynamic linear models defined by this function are very similar to those defined by \code{\link{EnsEBMtrend.modeldef}}, but with a more complex dependence structure among the ensemble members due to the two-way structure. The trend for an ensemble member generated using the \eqn{r}th RCM driven by the \eqn{g}th GCM is represented as
\deqn{\mu_{rg}(t) = \alpha\mu_{0}(t) + \tilde{\mu}_{\omega}(t) + \tilde{\mu}_{r\cdot}(t) + \tilde{\mu}_{\cdot g}(t) + \tilde{\mu}_{rg}(t)\ ,}{mu_{rg}(t) = alpha*mu_0(t) + mu.tilde_{omega}(t) + mu.tilde_{r.}(t)  + mu.tilde_{.g}(t) + mu.tilde_{rg}(t),}
incorporating discrepancy processes attributable to the individual RCMs and GCMs as well as to their interaction (the final term in the expression above, omitted from the analysis if \code{interactions="none"}). 

The parameters for these models are the same as those for \code{\link{EnsEBMtrend.modeldef}}. In particular, common "measurement error" and innovation variances are assumed for all ensemble members.

Again as in \code{\link{EnsEBMtrend.modeldef}}, an EBM-inspired representation is used to describe the evolution of each contribution to \eqn{\mu_{rg}(\cdot)}{mu_{rg}[.]} in the model definition above. For example, the real-world trend \eqn{\mu_{0}(\cdot)}{mu_0[.]} evolves as \eqn{\mu_{0}(t) = \phi_{0}\mu_{0}(t-1) + \gamma_{0} f(t) + \beta_{0}(t-1)}{mu_0[t] = phi_0*mu_0[t-1] + gamma_0*f[t] + beta_0[t-1]}, with corresponding expressions for all other contributions. 

When \code{constrain} is \code{FALSE} and \code{interactions} is \code{"all"}, the state vector for the model at time \eqn{t} is defined as 

\deqn{\mathbf{S}(t) = \left(\mathbf{S}'_{0}(t)\ \tilde{\mathbf{S}}'_{\omega}(t)\ \tilde{\mathbf{S}}'_{1\cdot}(t)\ \cdots\ \tilde{\mathbf{S}}'_{R\cdot}(t)\ \tilde{\mathbf{S}}'_{\cdot 1}(t)\ \cdots\ \tilde{\mathbf{S}}'_{\cdot G}(t)\ \tilde{\mathbf{S}}'_{11}(t)\ \cdots\ \tilde{\mathbf{S}}'_{RG}(t)\right)\  ,}{S[t] = (S'[0](t) S.tilde'[w](t) S.tilde'[1.](t) ... S.tilde'[R.](t) S.tilde'[.1](t) ... S.tilde'[.G](t) S.tilde'[11](t) ... S.tilde'[RG](t))' ,}

where:

\itemize{
\item{\eqn{\mathbf{S}_{0}(t) = \left(\mu_{0}(t)\ \beta_{0}(t)\ \gamma_{0}\right)'}{S[0](t) = (mu[0](t) beta[0](t) gamma[0])'}};
\item{\eqn{\tilde{\mathbf{S}}_{*}(t)= \left(\tilde{\mu}_{*}(t)\ \tilde{\beta}_{*}(t)\ \tilde{\gamma}_{*}\right)'}{S.tilde[*](t) = (mu.tilde[*](t) beta.tilde[*](t) gamma.tilde[*])'} (the asterisk \eqn{*} denoting a generic subscript);}
\item \eqn{R} and \eqn{G} are the numbers of RCMs and GCMs contributing to the ensemble (see documentation for the \code{Groups} argument above).
}

When \code{constrain} is \code{TRUE} (the default), the elements corresponding to RCM \eqn{R} and GCM \eqn{G} are removed: in this case they are defined in terms of the other elements by the sum-to-zero constraints. Similarly, the last group of terms is omitted when \code{interactions} is \code{"none"} (the default).

}
\value{A list defining a dynamical linear model, in the format required by routines in the \pkg{dlm} library.}
\references{
REFERENCES TO BE ADDED
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{EnsEBM2waytrendSmooth}} for a convenience wrapper function; \code{\link{EBMtrend.modeldef}}; \code{\link{EnsEBMtrend.modeldef}}}

