\title{Update history for package "TimSPEC"}
\name{NEWS}
\section{Changes in version 0.0-4 (2024-01-24)}{
  \itemize{
  \item{Added \code{EnsMean} routine with \code{EnsStruct} and \code{na.rm} options; amended \code{PlotEnsTS} and other routines to use this, thus allowing more flexible calculation of ensemble means.}
  \item{Improved calculation of initial values for maximum-likelihood / maximum a posteriori estimation in \code{Ens...} mimics.}
  \item{Added \code{C0} argument to \code{...modeldef} and smoothing routines, so that users can initialise the full covariance matrix for the Kalman Filter / Smoother directly (as well as the mean (previously any diffuse elements of the state vector were allocated the same initial variance given by the \code{kappa} argument, which lacks flexibility in a Bayesian interpretation whereby \code{m0} and \code{C0} are considered as prior specifications.)}
  \item{Updated documentation throughout, for consistency with CRAN conventions.}
  \item{Renamed some functions and their arguments, most notably those relating to \code{summary.dlmMLE} - again, for consistency with CRAN conventions.}
  \item{Added this \code{News.Rd} file.}
  }
}

\section{Changes in version 0.0-3}{

\itemize{
\item{And 0.0-2: fixed some initial teething problems}
}
}

\section{Changes in version 0.0-1 (2023-04-01)}{

\itemize{
\item{Initial version of package}
}

}