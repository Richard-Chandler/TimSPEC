\name{summary.dlmMLE}
\alias{summary.dlmMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarise a fitted model object}
\description{Produces a summary table of a dynamic linear model fit obtained using \code{\link{dlmMLE}} or \code{\link{dlm.SafeMLE}}. }
\usage{
\method{summary}{dlmMLE}(fit, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{An object of class \code{dlmMLE}, typically the result of a call to \code{\link{dlmMLE}},  \code{\link{dlm.SafeMLE}} or equivalent. \strong{Note} that this object class is defined within the \pkg{TimSPEC} package, not within \pkg{dlm}.}
  \item{print}{If \code{TRUE} (the default), results are written to screen. }
}
\value{A list containing the following components:

  \item{model.table}{Table of estimates and, if \code{fit} contains a \code{hessian} component (the result of calling \code{\link{dlmMLE}} or \code{\link{dlm.SafeMLE}} with argument \code{hessian=TRUE}) approximate standard errors}
  \item{logLik}{The unpenalised log-likelihood for the fit}
  \item{logLik.Pen}{The maximised log-likelihood for the fit. If penalised likelihood fitting was used e.g. by supplying a non-null \code{prior.pars} argument to \code{\link{dlm.SafeMLE}}, then this is the maximum \emph{penalised} log-likelihood; otherwise it's the maximum unpenalised version (which, in this case, is the same as the \code{logLik} component of the result).}
  \item{Corr}{If \code{fit} contains a \code{hessian} component, approximate correlation matrix of the parameter estimates. }
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
