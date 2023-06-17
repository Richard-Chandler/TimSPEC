\name{dlm.SafeLL}
\alias{dlm.SafeLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{"Safe" computation of a (penalised) log-likelihood for a dynamic linear model. 
}
\description{This is a wrapper for the \code{\link{dlmLL}} routine in the \pkg{dlm} 
library. It (a) handles situations where \code{\link{dlmLL}} returns an error or an infinite  value (b) is designed for use in calls to \code{\link{nlm}} rather than \code{\link{optim}} when estimating parameters (c) allows for maximisation of a \emph{penalised} likelihood which can be regarded as an approximation to a Bayesian posterior.
}
\usage{
dlm.SafeLL(theta, Y, build, prior.pars = NULL, BigVal = 1e+12, 
           debug = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{The parameter vector}
  \item{Y}{The data matrix}
  \item{build}{A function that will construct the \code{dlm} object representing the model for which the (penalised) log-likelihood is required. Typically one of the \code{.modeldef} functions in the \pkg{TimSPEC} package.}
  \item{prior.pars}{If non-\code{NULL}, a two-column matrix with number of rows equal to the length of \code{theta}: the columns contain respectively the means and standard deviations of independent Gaussian priors on the elements of \code{theta}. 
  
  If \code{NULL} (the default), no prior information is used. 
}
  \item{BigVal}{Value to return in case \code{\link{dlmLL}} fails. This should be much larger than any reasonably expected value of the negative (penalised) log-likelihood, but not so large as to cause numerical instability when evaluating gradients etc.
}
  \item{debug}{As in \code{\link{dlmLL}}
}
  \item{\dots}{Other arguments to \code{build}}
}
\value{If \code{prior.pars} is \code{NULL} then the function returns the result of \code{\link{dlmLL}} when that routine returns a finite number, otherwise it returns \code{BigVal}. If \code{prior.pars} is non-null, the function returns the negative log posterior (or penalised negative log-likelihood if you prefer) when \code{\link{dlmLL}} returns a finite number. In this case the result also has an "\code{Unpenalised}" attribute containing the unpenalised negative log-likelihood: this can be helpful for comparing (informally) fits of models with different parameters & hence different priors. 
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>
}

\seealso{\code{\link{dlmLL}}, \code{\link{dlm.SafeMLE}}
}
