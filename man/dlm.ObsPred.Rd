\name{dlm.ObsPred}
\alias{dlm.ObsPred}
\title{Prediction for observable quantities in a dynamic linear model
}
\description{To calculate the predicted values, and associated error variances, for \emph{observable} quantities (i.e. components of the observation vector) according to a state space time series model that has been used to run a Kalman Smoother.}
\usage{
dlm.ObsPred(ModelBundle)
}
\arguments{
  \item{ModelBundle}{A list containing at least the named components \code{Model} (an object of class \code{dlm}) and \code{Smooth} (the result of Kalman Smoothing a dataset using the dlm stored in \code{Model}, via the \code{\link{dlmSmooth}} command). It is assumed that the model is set up in such a way that the "real world" value is the first element of each observation vector.}
}
\value{The routine returns a list containing the predictions and their standard errors, each returned as a matrix in which the rows represent time points and columns represent series.
}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
