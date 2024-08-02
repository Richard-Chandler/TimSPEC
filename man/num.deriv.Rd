\name{num.deriv}
\alias{num.deriv}
\alias{num.hess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Numerical differentiation of a function}
\description{To compute the first and second derivatives of a function with respect to its first argument.
}
\usage{
num.deriv(f, x0, xlo=-Inf, xhi=Inf, steps=NULL, shrink=1.4,
					maxit=10, relerr=1e-6, QuitEarly=NULL, ...)
num.hess(f, x0, xlo=-Inf, xhi=Inf, steps=NULL, shrink=1.4, 
				 maxit=10, relerr=1e-6, QuitEarly=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{The function to be differentiated. This should be an R function with a scalar result: differentiation is with respect to its (vector) first argument.}
  \item{x0}{The point at which to evaluate the derivative. This is a vector, corresponding to the first argument of the function \code{f}.}
  \item{xlo}{Vector of lower bounds on the first argument of \code{f}. This can be either a scalar, in which case it is replicated across all elements, or a vector the same length as \code{x0}. Defaults to \code{-Inf}, indicating that the elements of the function argument are unbounded below.}
  \item{xhi}{Vector of upper bounds on the first argument of \code{f}, similarly. Defaults to \code{Inf}, indicating that the elements of the function argument are unbounded above.}
  \item{steps}{Vector of largest step sizes to consider initially when evaluating the derivatives. For \code{num.deriv} this defaults to the larger of \code{abs(x0)/10^4} and \code{0.001}; for \code{num.hess} it is approximately the square root of this quantity.}
  \item{shrink}{Amount by which to shrink the step size at each iteration (see details below), to obtain an improved estimate.}
  \item{maxit}{Maximum number of iterations.}
  \item{relerr}{Scalar indicating the relative error required for each element of the result.}
  \item{QuitEarly}{Either \code{NULL} (the default) or a list containing two elements, in which case it defines conditions for the routine to exit early if a "surprising" value of \code{f} is found. In the latter case, the named elements of the list should be \code{reference} (a numeric value against which to compare the values of \code{f}) and \code{greater} (a logical scalar indicating whether to quit if a value \code{f>reference} is found - if \code{FALSE}, the function will quit if a value \code{f<reference} is found).}
  \item{...}{Other arguments to \code{f}. None of these should have names matching other arguments of \code{num.deriv} or \code{num.hess}.}
}

\details{These routines use Ridders' implementation of Richardson's method - see Press et al. (1992), Section	5.7 for example. The implementation here is basically that	given as routine \code{dfridr} by Press et al., translated from \code{Fortran} to \code{R} and adapted to handle boundaries (in which case a right- or left-hand derivative is evaluated as appropriate) and vector arguments. One-sided derivatives are calculated by reflecting the function about the boundary, which is equivalent to using function evaluations on one side of the boundary only - and correspondingly less accurate. For one-sided derivatives it may therefore be necessary to decrease the default value of steps to obtain the required accuracy.
}

\value{A list, containing elements \code{gradient} or \code{hessian} (for \code{num.deriv} and \code{num.hess} respectively), \code{error} and \code{details}. \code{gradient} and \code{hessian} are estimates of the gradient vector and hessian matrix respectively, while \code{error} is an estimate of the magnitude of the error in each element of \code{gradient} (resp. \code{hessian}). \code{details} is a data frame recording the function evaluations used in the computations: it contains \eqn{p+2} columns where \eqn{p} is the length of \code{x0}. The first column is the iteration number in the algorithm (note that each iteration involves several function evaluations); the second is the value of the function \code{f}; and the remaining columns contain the value of the function argument. If \code{x0} has names, the corresponding columns of \code{details} will have the same names; otherwise they will be names as \code{x[1]}, \code{x[2]} etc.}
\author{
Richard E. Chandler <r.chandler@ucl.ac.uk>
}
\references{
Press, W., Teukolsky, S., Vetterling, W., and Flannery, B. (1992). \emph{Numerical Recipes in FORTRAN} (second edition). Cambridge University Press.
}
\examples{
QuadFun <- function(xx, coef) {
 #
 # Define a quadratic function of two variables
 #
 x <- xx[1]
 y <- xx[2]
 z <- (coef[1]*x^2) + (coef[2]*y^2) + (coef[3]*x*y) + 
      (coef[4]*x) + (coef[5]*y) + (coef[6])
 z
}
#
# Set some coefficients and evaluate the derivatives at 
# the point (1,1)
#
coef <- c(6, 3, -2, 7, -12, 4)
num.deriv(QuadFun, c(1,1), coef=coef)
num.hess(QuadFun, c(1,1), coef=coef)
}
