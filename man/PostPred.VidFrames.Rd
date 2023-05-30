\name{PostPred.VidFrames}
\alias{PostPred.VidFrames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Animated samples from a posterior predictive distribution for a time series}
\description{To produce the video frames for an animation of samples from a postprocessed ensemble and, optionally, use the frames to produce an animated GIF image}
\usage{
PostPred.VidFrames(Data, Samples, AntiLog = FALSE, DatColours, DatTypes = c(1, 1), 
                   PredColours, Nsamples = 20, Nsteps = 5, replace = TRUE, Folder, 
                   FileNameRoot = Folder, WarnExisting = TRUE, PNG.args = list(), 
                   MakeGif = TRUE, DelFrames = FALSE, Annotation = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{Matrix or data frame containing time (i.e. year) in column 1, real-world series in column 2 and ensemble members in the remaining columns.}
  \item{Samples}{A list containing at least named components \code{Obs} and \code{Weights} --- as produced by, for example, \code{\link{PostPredSample}}. }
  \item{AntiLog}{If \code{TRUE}, samples will be antilogged before plotting (appropriate when models were fitted on a log scale) }
  \item{DatColours}{Vector of colours used to plot the ensemble members: see \code{\link{PlotEnsTS}}. }
  \item{DatTypes}{Vector of line types, similarly.}
  \item{PredColours}{Vector of two colours: the first is the base colour for plotting objects relating to the real-world trend \eqn{\mu_{0}(\cdot)}{mu[0](.)}, and the second for plotting the ensemble consensus. Only the first element is used here: the second element is retained so that (for example) the same object can be used for this routine as for calls to \code{\link{SmoothPlot}}. }
  \item{Nsamples}{Number of "primary" samples to use to construct the video frames. }
  \item{Nsteps}{Number of interpolation points between each primary sample. }
  \item{replace}{Logical scalar, controlling whether or not to draw primary samples with replacement. }
  \item{Folder}{Folder within which to store the output graphics files. }
  \item{FileNameRoot}{Root of filename for output graphics files. E.g. if \code{FileNameRoot} is \code{"fish"} then the files will be named \code{fish001.png}, \code{fish002.png}, ... . Defaults to the name of \code{Folder}. }
  \item{WarnExisting}{Logical scalar indicating whether to warn if \code{Folder} exists already. }
  \item{PNG.args}{List of arguments to the \code{\link{png}} command, used to control (e.g.) the size of the images produced. }
  \item{MakeGif}{If \code{TRUE}, a GIF animation will also be created. }
  \item{DelFrames}{If \code{TRUE} and if \code{MakeGif} is \code{TRUE}, the files containing the individual video frames will be deleted after the GIF is created. }
  \item{Annotation}{If \code{TRUE}, plots / frames will be annotated with copyright and funding information using \code{\link{GraphAnnot}}. }
  \item{\dots}{Other arguments to \code{\link{SmoothPlot}} }
}
\details{This function implements the suggestion of Bowman (2019), to visualise uncertainty in complex objects using animations. Each individual frame in an animationshows a sample from the posterior distribution (strictly speaking, from a joint distribution with the same mean and covariance as the posterior); but the samples are consgtructed in such a way that the transitions between frames are smooth. The smoothness of the transitions helps the viewer to process the information as it is presented. 

The animations are produced by sampling a small number of 'primary' frames from the posterior predictive distribution for quantities of interest (here, typically a real-world time series), and then interpolating between these primary frames in such a way that each interpolation also has the correct mean and covariance structure. 
}
\value{A matrix containing the time series used for each frame of the animation: the column names are \code{Frame1.0}, \code{Frame1.1}, ..., \code{Frame2.0} etc., where frames ending in 0 are the 'primary' frames and the others are the interpolations.}
\references{Bowman, A. W. (2019). Graphics for uncertainty. Appl. Statist., 182(2):403-418. doi:
\doi{10.1111/rssa.12379}}
\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\seealso{\code{\link{SmoothPlot}} }
