\name{EBMtrend.modeldef}
\alias{EBMtrend.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{State space model for a single time series, inspired by Energy Balance Models}
\description{Sets up the structure of a state space model for a single time series with one value per year, inspired by the structure of a simple energy balance model (EBM)}
\usage{
EBMtrend.modeldef(theta, Xt, m0 = NULL, kappa = 1e+06, UsePhi = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{Parameter vector. If \code{UsePhi} is \code{TRUE} (see below) then this contains three elements: the log of the "measurement error" variance (\eqn{\sigma^{2}}{sigma^2} in the notation below), the log of the "trend drift" innovation variance \eqn{\tau^{2}}{tau^2}, and the logit of the dependence parameter \eqn{\phi}{phi}. If \code{UsePhi} is \code{FALSE} then \code{theta} contains just the first two of these elements.}
  \item{Xt}{Vector containing a time series of effective forcings}
  \item{m0}{Optional vector of initial values for the state vector. If \code{NULL} (the default) this is determined automatically from the model structure.}
  \item{kappa}{Large" value giving the variance used to initialise diffuse elements of the state vector.}
  \item{UsePhi}{Controls whether or not to include a thermal inertia parameter in the model structure.}
}
\details{The model is defined by the equation
\deqn{Y(t) = \mu(t) + \varepsilon(t)}{Y(t) = mu[t] + e[t]}
with \eqn{E[\varepsilon(t)]=0}{E[e(t)]=0} and \eqn{Var[\varepsilon(t)]=\sigma^{2}}{Var[e(t)]=sigma^2}, and where the trend \eqn{\mu(t)}{mu(t)} follows the EBM-inspired dynamics
\deqn{\mu(t) = \phi\mu(t-1) + \gamma f(t) + \beta(t-1)\ .}{mu(t) = phi*mu(t-1) + gamma*f(t) + beta(t-1)\ .}
Here, \eqn{f(t)} represents the "forcing" at time \eqn{t} and \eqn{(\beta(\cdot))} is a random walk "drift" with innovation variance \eqn{\tau^{2}}{tau^2} representing the error arising from the use of the EBM approximation to represent the trend. The state vector at time \eqn{t} is defined as
\deqn{\bm{S}(t) = \left(\mu(t)\ \beta(t)\ \gamma\right)'}{S(t) = (mu(t)  beta(t)  gamma)'}

If \code{UsePhi} is \code{FALSE}, the value of \eqn{\phi}{phi} in the model definition is set to zero.
}
\value{A list defining a dynamical linear model, in the format required by routines in the \pkg{dlm} library.
}
\references{
REFERENCES TO BE ADDED
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
