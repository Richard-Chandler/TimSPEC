\title{Update history for package "TimSPEC"}
\name{NEWS}

\section{Changes in version 0.0-9 (2024-06-13)}{
  \itemize{
    \item{Added \code{theta} argument to the \code{*Smooth} wrapper functions, allowing users to supply initial values for the optimisation routine in the first step.}
    \item{Minor bug fixes in \code{SmoothPlot}, relating to plotting of confidence intervals when \code{NonNeg} is \code{TRUE} and some intervals go negative.}
    \item{Corrected code for calculation of ensemble consensus in \code{Ens*modeldef} functions - length of \code{Consensus} component of result wasn't always correct.}
	\item{Improved documentation of some functions, correcting typos and clarifying mathematical notation.}
  }
}

\section{Changes in version 0.0-8 (2024-05-28)}{
  \itemize{
    \item{Replaced out-of-date package acronyms in documentation with \code{TimSPEC} throughout.}
    \item{Added \code{Legend} argument to \code{SmoothPlot}, allowing users to suppress the default addition of legends to plots.}
    \item{Minor change to the output of \code{PostPredSample}, allowing the user to identify the entries of the \code{States} component that were replaced if \code{ReplaceOnFail} is \code{TRUE}.}
    \item{Amended the \code{ExamineFailures} function so that it still
    produces a plot (with a warning message) even if there are no
    failures.}
  }
}

\section{Changes in version 0.0-7 (2024-02-21)}{
  \itemize{
    \item{Corrected definition of \code{Consensus} elements of ensemble mimic definitions (functions \code{EnsSLLT.modeldef}, \code{EnsEBMtrend.modeldef} and \code{EnsEBM2waytrend.modeldef})}
    \item{Minor change to initialisation of parameters in estimation of mimics for ensembles.}
    \item{Amended definition of \code{messages} argument to \code{dlm.SafeMLE} and other functions, to allow access to \code{print.level} control in \code{nlm} when fitting models. This can be helpful for diagnosing problems in estimation.}
    \item{Added details of mimic definitions and state vectors to all of the \code{.modeldef} documentation files.}
    \item{Updated \code{seealso} sections in documentation to improve navigation.}
  }
}

\section{Changes in version 0.0-6 (2024-02-09)}{
  \itemize{
    \item{Minor bug fixes and amendments to documentation.}
  }
}

\section{Changes in version 0.0-5 (2024-02-07)}{
  \itemize{
    \item{Amended \code{EnsMean} function to cope with a broader range of missing value patterns, and renamed some of its arguments to reflect this and for consistency with object names elsewhere in the package.}
  }
}

\section{Changes in version 0.0-4 (2024-01-29)}{
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