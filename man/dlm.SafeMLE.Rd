\name{dlm.SafeMLE}
\alias{dlm.SafeMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{"Safe" maximum likelihood estimation for a dynamic linear model}
\description{This replaces the \code{\link{dlmMLE}} routine in the \pkg{dlm} library, which can be unstable. It uses \code{\link[stats]{nlm}} rather than \code{\link{optim}} for the optimisation.
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
  \item{Use.dlm}{If \code{TRUE}, the function will use \code{\link{dlmMLE}} in place of \code{\link[stats]{nlm}} for the optimisation. This allows easy comparison of the results from the two optimisation procedures.}
  \item{par.names}{Optional character vector, the same length as \code{theta.init}, that will be used to label the parameters in the output}
  \item{prior.pars}{If non-\code{NULL}, a two-column matrix with number of rows equal to the length of \code{theta.init}: the columns contain respectively the means and standard deviations of independent Gaussian priors on the individual parameters.} 
  \item{messages}{Either a logical or numeric scalar, controlling the output of messages to screen during fitting when \code{Use.dlm} is \code{FALSE}. If \code{FALSE} or less than 1, no messages will be produced. If \code{TRUE} or equal to 1, basic messages will be produced. If numeric and greater than 1, the value of \code{messages-1} will be passed as the \code{print.level} argument to \code{\link[stats]{nlm}}.
  
  \strong{Note:} if \code{messages} is numeric and greater than 1, it's quite likely that the routine will issue repeated messages `\code{Iteration limit exceeded. Algorithm failed.}'. This is nothing to worry about: it arises because it calls \code{\link[stats]{nlm}} repeatedly with a relatively small number of iterations, at each stage updating the \code{\link[stats]{nlm}} control parameters to try and stabilise the subsequent calculations. The `Iteration limit exceeded' message comes from \code{\link[stats]{nlm}} on exit from each of these intermediate calls. 
  }
  \item{hessian}{Controls whether to compute the hessian on exit from \code{\link[stats]{nlm}}. The internal calculations within \code{\link[stats]{nlm}} itself can be inaccurate, so this is done via a call to \code{\link[numDeriv]{hessian}} from the \pkg{numDeriv} package.}
  \item{\dots}{Other arguments passed to \code{build} and \code{\link[stats]{nlm}}}
}
\details{If \code{Use.dlm} is \code{FALSE} and \code{prior.pars} is non-\code{NULL}, the routine minimises the negative log posterior density (equivalently the negative penalised log-likelihood) for the parameters; in this case, the value of \code{prior.pars} is also included in the list result (see below) for use by routines such as \code{\link{summary.dlmMLE}}. Otherwise it minimises the unpenalised log-likelihood. 

If \code{use.dlm} is \code{TRUE} then non-\code{NULL} values of \code{prior.pars} are ignored, with a warning. 

It should be noted that as with any numerical optimisation process, this routine is not guaranteed to find a global maximum of the (penalised) log-likelihood. It incorporates some levels of protection but, particularly for more complex mimics and with large ensembles, it can sometimes fail to locate the global maximum without giving any obvious warning. Some further protection can be gained by taking posterior samples from the Laplace approximation (e.g. using \code{\link{dlm.ThetaSample}}) and then calculating importance weights using \code{\link{dlm.ImportanceWts}} which will issue a warning if it finds any samples with higher (penalised) log-likelihoods. Equivalently, the \code{CheckMax} argument to \code{\link{PostPredSample}} can be used to automate this check when the \code{Importance} argument is \code{TRUE}.
}
\value{If \code{Use.dlm} is \code{TRUE} then the routine returns the result of a call to \code{\link{dlmMLE}}, with various control parameters set to stabilise the optimisation as far as possible. If \code{Use.dlm} is \code{FALSE}, the routine returns the list result from a call to \code{\link[stats]{nlm}} (see below for what is optimised in this case), but with the names of key components modified to match the results of a call to \code{\link{optim}} (which is used by \code{\link{dlmMLE}}. The exception is the \code{code} component of the result, which is retained and a \code{convergence} component is added which is compatible with the convergence codes from \code{\link{optim}}. This ensures that the routine can be run interchangeably with \code{\link{dlmMLE}}. A component \code{HessFail} is also added, to record whether the Hessian had to be tweaked to make it positive definite. 

The result is assigned a class of \code{dlmMLE}, in order that the appropriate \code{\link[=summary.dlmMLE]{summary}} method can be used (see below).

}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dlmMLE}}, \code{\link{dlm.SafeLL}}, and \code{\link{summary.dlmMLE}} for the \code{summary} method.}
