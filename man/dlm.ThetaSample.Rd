\name{dlm.ThetaSample}
\alias{dlm.ThetaSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate posterior sampling for dynamic linear models}
\description{Generates samples from the approximate posterior distribution for the parameter vector of a dynamic linear model}
\usage{
dlm.ThetaSample(Fit, N, Random = TRUE, Quantile = TRUE, 
                Level = if (!Quantile) 1.96 else 0.99, 
                  df = NULL, Antithetic = c("Mean", "Scale"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fit}{A list containing the results of a numerical maximisation of a log-posterior. Typically the result of a call to \code{\link{dlm.SafeMLE}} or equivalent.}
  \item{N}{Number of samples to draw. Ignored if \code{Random} is \code{FALSE} (see below).}
  \item{Random}{If \code{TRUE} (the default) then samples are drawn randomly from the approximate posterior. Otherwise points are taken from the principal axes of the ellipsoid defined by the posterior covariance matrix: see "Details" below.
}
  \item{Quantile}{Logical scalar, used only if \code{Random} is \code{FALSE}. In this case it determines the behaviour of \code{Level} as described in the  "Details" section.}
  \item{Level}{Numeric scalar, used to control output when \code{Random} is \code{FALSE}. See "Details" section.}
  \item{df}{If non-\code{NULL} and if \code{Random} is \code{TRUE}, sampling will be from a scaled and shifted multivariate \eqn{t}-distribution with \code{df} degrees of freedom. In this case, the "mean" and "covariance matrix" of the posterior distribution are taken as the location and dispersion matrices of the \eqn{t}-distribution.}
  \item{Antithetic}{Controls use of antithetic sampling to reduce sampling variability if \code{Random} is TRUE. If \code{NULL} then independent random sampling is used); otherwise it should be a character vector containing a subset of \code{c("Mean", "Scale", "Corr")} to balance for location, scale and dependence respectively --- see "Details" section.The argument is ignored if \code{Random} is \code{FALSE}.  
  %%     ~~Describe \code{Antithetic} here~~
}
}
\details{Most users won't need to call this function directly: use \code{\link{PostPredSample}} instead, which is a wrapper to all of the function calls needed to produce a posterior predictive sample. However, some of the arguments to this function are passed through \code{\link{PostPredSample}}; some key points are summarised here.

\subsection{Behaviour for non-random sampling}{

When \code{Random} is \code{FALSE}, the routine returns the values \eqn{m}, \eqn{m \pm A*u[1]*v[,1], m \pm A*u[2]*v[,2], ... , m  \pm A*u[d]*v[,d]}{m+/-A*u[1]*v[,1], m+/-A*u[2]*v[,2], ... , m+/-A*u[d]*v[,d]}, where \eqn{m} is theposterior mode, \eqn{d} is the number of parameters in the model, \eqn{u} is a vector containing the square roots of the singular values of the posterior covariance matrix and \eqn{v} is a matrix containing the corresponding singular vectors in its columns. \eqn{A} is set using the values of \code{Quantile} and \code{Level}. Specifically: 

\itemize{
 \item If \code{Quantile} is \code{FALSE} then \code{Level} is itself the value of \eqn{A} and the default value of 1.96 delineates a 95\% componentwise interval along each of the principal axes of the underlying covariance matrix. 
 \item If \code{Quantile} is \code{TRUE} then \eqn{A} is chosen so that the  sampled points all lie on the ellipsoid defining the approximate joint (Gaussian) central credible region for the parameter vector: in this case, the default value of 0.99 specifies the 99\% credible region.
}
}
\subsection{Use of antithetic sampling}{

When \code{Random} is \code{TRUE}, the argument \code{Antithetic} is used to control the use of antithetic variables which may help to reduce sampling variability in the posterior. If \code{NULL}, independent random sampling is used. Otherwise:

\itemize{
\item If the argument contains \code{"Mean"} and not \code{"Corr"} then, for every sampled value \eqn{\theta}{theta}, the routine also returns \eqn{2m-\theta}{2m-theta} which is the reflection about the posterior mode: this ensures that the mean of the sampled values is equal to the posteriormean. 
\item If the argument contains \code{"Corr"} then the routine produces \eqn{2^d-2} structured samples for each independent value, by projecting the initial value onto the principal axes of variation about the mode and then multiplying each centred component by \eqn{\pm 1}{+/-1} (the single antithetic sample generated using "Mean" corresponds to the case when all the multipliers are 1): this distributes samples around the ellipsoid with joint density equal to that ofthe initial value, thus ensuring that the correlation shape of the underlying distributionis preserved. 
\item Finally, if the argument contains \code{Scale} then, for each sample \eqn{\theta}{theta} (after potentially generating the \code{Mean} and \code{Corr} antithetics), the value \eqn{m + \sqrt{q2/q1}(\theta-m)}{m + sqrt(q2/q1)(theta-m)} is also returned where \eqn{q1} is the quadratic form giving the chi-squared quantile for the joint distribution of \eqn{\theta-m}{theta-m}, and \eqn{q2} is the opposite tail quantile for the same distribution: this ensures that samples have the correct overall dispersion about the posterior mode.
}
}

If the antithetic options are used, the total number of samples produced is a multiple of a power of 2: if this is greater than \code{N} then a random subset of \code{N} samples is selected. 
}
\value{A matrix in which each row contains a posterior sample. The matrix has additional attributes \code{Mean}, \code{Cov.SpD} and \code{df} giving the mean vector, covariance (dispersion) matrix spectral decomposition and degrees of freedom for the underlying distribution. 
}
\references{Durbin, J. and Koopman, S. J. (2012). \emph{Time series analysis by state space methods (second edition)}, Chapter 11. Oxford University Press, Oxford.}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>
}
