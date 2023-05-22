\encoding{utf8}
\name{GraphAnnot}
\alias{GraphAnnot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add copyright information to a plot}
\description{Adds copyright information and acknowledgements in the lower corners of a plot}
\usage{
GraphAnnot(Copyright = "© UCL 2022", Acknowl = "Funded by the UK Climate Resilience Programme")
}
\arguments{
  \item{Copyright}{Text to be added at bottom left of plot, typically a copyright statement.}
  \item{Acknowl}{Text to be added at bottom right, typically acknowledgements.}
}
\value{The function merely adds the required text to an existing plot, and returns a \code{NULL} result. }
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
