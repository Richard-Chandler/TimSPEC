\name{summary.dlmMLE}
\alias{print.dlmMLE}
\alias{summary.dlmMLE}
\alias{print.summary.dlmMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarise a fitted model object}
\description{Methods for dynamic linear model fits obtained using \code{\link[dlm]{dlmMLE}} or \code{\link{dlm.SafeMLE}}. }
\usage{
\method{print}{dlmMLE}(x, ...)
\method{summary}{dlmMLE}(object, ...)
\method{print}{summary.dlmMLE}(x, digits=4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, object}{An object of class \code{dlmMLE} or \code{summary.dlmMLE}, typically the result of a call to \code{\link[dlm]{dlmMLE}},  \code{\link{dlm.SafeMLE}} or the resulting \code{summary}. \strong{Note} that these object classes are defined within the \pkg{TimSPEC} package, not within \pkg{dlm}.}
  \item{digits}{Controls the number of significant digits to display in the output when printing summaries.}
  \item{...}{Additional arguments to the \code{\link[base]{print}} and \code{\link[base]{summary}} generics}
}
\value{The \code{print} method returns \code{object}, invisibly. The \code{summary} method returns a list containing the following components:

  \item{model.table}{Table of estimates and, if \code{fit} contains a \code{hessian} component (the result of calling \code{\link[dlm]{dlmMLE}} or \code{\link{dlm.SafeMLE}} with argument \code{hessian=TRUE}) approximate standard errors}
  \item{logLik}{The unpenalised log-likelihood for the fit}
  \item{logLik.Pen}{The maximised log-likelihood for the fit. If penalised likelihood fitting was used e.g. by supplying a non-null \code{prior.pars} argument to \code{\link{dlm.SafeMLE}}, then this is the maximum \emph{penalised} log-likelihood; otherwise it's the maximum unpenalised version (which, in this case, is the same as the \code{logLik} component of the result).}
  \item{Corr}{If \code{fit} contains a \code{hessian} component, approximate correlation matrix of the parameter estimates. }
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
