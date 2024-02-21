\name{SLLT.modeldef}
\alias{SLLT.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth local linear trend model for a single time series}
\description{Sets up the structure of a "smooth local linear trend" state space
model for a single time series.}
\usage{
SLLT.modeldef(theta, m0 = NULL, C0=NULL, kappa = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{Vector of length 2: \code{theta[1]} is the log of the "measurement error" variance, and \code{theta[2]} is the log of the innovation variance on the slope process.}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{C0}{Optional covariance matrix for initialising the state vector. If \code{NULL} (the default) this is determined automatically from the model structure and from \code{kappa} (see below).}
  \item{kappa}{Initialisation variance for diffuse elements of the state vector. Ignored if \code{C0} is non-\code{NULL}.}
}
\details{The model is defined by the equations
\deqn{Y(t) = \mu(t) + \varepsilon(t)}{Y(t) = mu[t] + e[t]}
\deqn{\mu(t) = \mu(t-1) + \beta(t-1)}{mu(t) = mu(t-1) + beta(t-1)}
and
\deqn{\beta(t) = \beta(t-1) + z[t]\ ,}{beta(t) = beta(t-1) + z[t] ,}
where \eqn{E[\varepsilon(t)]=E[z(t)]=0}{E[e(t)]=E[z(t)]0}, \eqn{Var[\varepsilon(t)]=\sigma^{2}}{Var[e(t)]=sigma^2} and \eqn{Var[z(t)]=\tau^{2}}{Var[z(t)]=tau^2}, independently of each other and for each \eqn{t}.

At time \eqn{t}, the state vector for the model is defined as 

\deqn{\mathbf{S}(t) = \left(\mu(t)\ \beta(t)\right)'.}{S[t] = (mu[t] beta[t])'.}
}


\value{A list defining a dynamical linear model, in the format required by routines in the \pkg{dlm} library.}
\references{Chandler, R. E. and Scott, E. M. (2011). \emph{Statistical methods for trend detection and analysis in the environmental sciences}, Section 5.5. John Wiley & Sons, Chichester.}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{SLLTSmooth}} for a convenience wrapper function; \code{\link{EnsSLLT.modeldef}} for an extension to climate model ensembles; \code{\link{EBMtrend.modeldef}}, for a physically motivated alternative for application to climate time series.}
