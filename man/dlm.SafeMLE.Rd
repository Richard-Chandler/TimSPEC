\name{dlm.SafeMLE}
\alias{dlm.SafeMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{"Safe" maximum likelihood estimation for a dynamic linear model}
\description{This replaces the \code{\link{dlmMLE}} routine in the \pkg{dlm} library, which can be unstable. It uses \code{\link{nlm}} rather than \code{\link{optim}} for the optimisation.
}
\usage{
dlm.SafeMLE(theta.init, Y, build, debug = FALSE, Use.dlm = FALSE, 
            par.names = NULL, prior.pars = NULL, messages = TRUE, 
            hessian = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta.init}{Initial value for the parameter vector}
  \item{Y}{The data}
  \item{build}{A function that will construct the \code{dlm} object representing the model being fitted. Typically one of the \code{.modeldef} functions in the \pkg{TimSPEC} package.}
  \item{debug}{As in \code{\link{dlmLL}}}
  \item{Use.dlm}{If \code{TRUE}, the function will use \code{\link{dlmMLE}} in place of \code{\link{nlm}} for the optimisation. This allows easy comparison of the results from the two optimisation procedures.}
  \item{par.names}{Optional character vector, the same length as \code{theta.init}, that will be used to label the parameters in the output}
  \item{prior.pars}{If non-\code{NULL}, a two-column matrix with number of rows equal to the length of \code{theta.init}: the columns contain respectively the means and standard deviations of independent Gaussian priors on the individual parameters.} 
  \item{messages}{Controls output of messages to screen during fitting}
  \item{hessian}{Controls whether to compute the hessian on exit from \code{\link{nlm}}. The internal calculations within \code{\link{nlm}} itself can be inaccurate, so this is done via a call to \code{link{hessian}} from the \pkg{numDeriv} package.}
  \item{\dots}{ther arguments passed to \code{build} and \code{\link{nlm}}}
}
\details{If \code{Use.dlm} is \code{FALSE} and \code{prior.pars} is non-\code{NULL}, the routine minimises the negative log posterior density (equivalently the negative penalised log-likelihood) for the parameters; in this case, the value of \code{prior.pars} is also included in the list result (see below) for use by routines such as \code{\link{PrintDLMFit}}. Otherwise it minimises the unpenalised log-likelihood. 

If \code{use.dlm} is \code{TRUE} then non-\code{NULL} values of \code{prior.pars} are ignored, with a warning. 
}
\value{If \code{Use.dlm} is \code{TRUE} then the routine returns the result of a call to \code{\link{dlmMLE}}, with various control parameters set to stabilise the optimisation as far as possible. If \code{Use.dlm} is \code{FALSE}, the routine returns the list result from a call to \code{\link{nlm}} (see below for what is optimised in this case), but with the names of key components modified to match the results of a call to \code{\link{optim}} (which is used by \code{\link{dlmMLE}}. The exception is the \code{code} component of the result, which is retained and a \code{convergence} component is added which is compatible with the convergence codes from \code{\link{optim}}. This ensures that the routine can be run interchangeably with \code{\link{dlmMLE}}. A component \code{HessFail} is also added, to record whether the Hessian had to be tweaked to make it positive definite. 
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>
}
\note{Users will rarely need to call this routine directly: it is called by the wrapper \code{...Smooth} functions.}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dlmMLE}}, \code{\link{dlm.SafeLL}}}
