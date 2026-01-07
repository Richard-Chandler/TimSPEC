######################################################################
######################################################################
######################################################################
#
# Mimic definitions and utility routines for time series outputs from
# climate model ensembles. The routines in this script are as follows:
#
# CompactEns            Replaces a matrix of time series from an 
#                       ensemble with a smaller matrix in which 
#                       members from the same source have been 
#                       averaged. 
# CumulativeWeightPlot  Plots the cumulative contributions of a vector 
#                       of weights, sorted in descending order.
# dlm.ImportanceWts     Calculate importance weights for a sample
#                       from the approximate Gaussian posterior 
#                       distribution of parameters for a dlm. 
# dlm.ObsPred           Calculates predicted values and error variances for
#                       the observed quantities in a dynamic linear model
#                       that has been used to run a Kalman Smoother. 
# dlm.SafeLL            "Safe" log-likelihood computation for a dynamic 
#                       dlm library. 
# dlm.SafeMLE           "Safe" maximum-likelihood fitting, which can use
#                       dlmMLE or alternatives. Also offers the option
#                       to use maximum a posterior fitting with 
#                       independent Gaussian priors on the parameters
#                       (equivalently: penalised likelihood fitting).
# dlm.ThetaSample       Samples from Gaussian approximation to posterior
#                       for dlm parameters. Includes options for antithetic
#                       and deterministic sampling. 
# dmvnorm.SpD, dmvt.SpD Calculate density functions of multivariate 
#                       normal and t-distributions in which the 
#                       covariance / dispersion matrices are specified
#                       via their spectral / eigendecompositions. 
# EBMtrend.modeldef,    The modeldef function sets up a dynamic linear model 
# EBMtrendSmooth        structure for a single time series, following 
#                       dynamics "inspired by" a simple energy balance  
#                       model. The Smooth function fits the model and 
#                       implements the Kalman Smoother.
# EnsEBMtrend.modeldef, The same as the EBMtrend[...] functions, but 
# EnsEBMtrendSmooth     now for simultaneous analysis of an "observed"
#                       time series together with a set of exchangeeable
#                       ensemble members.
# EnsImpute             Impute missing values in an ensemble
# EnsMean               Calculates mean series for an ensemble, possibly
#                       accounting for ensemble structure. 
# EnsSLLT.modeldef,     Analysis of observed time series and exchangeable
# EnsSLLT0.modeldef,    ensemble members, using a smooth local linear
# EnsSSLLTSmooth,       trend model. The difference between the two modeldef
#                       constructor functions is that SLLT0 constrains the
#                       discrepancies between ensemble and observations
#                       to be constant in time. The Smooth function handles
#                       both model versions
# EnsEBM2waytrend.modeldef, Define and implement an EBM-inspired dynamic
# EnsEBM2waytrendSmooth  linear model for simultaneous analysis of 
#                        observations and simulations from a regional
#                       climate model ensemble, in which each ensemble
#                       member is produced using one of a set of 
#                       regional climate models with boundary conditions
#                       determined by one of a set of global models. 
# ExamineFailures       Produces a pairs plot of samples from the 
#                       posterior distribution of the parameter vector
#                       for a dlm, highlighting cases for which the
#                       sampling of the corresponding state vectors
#                       failed. 
# GraphAnnot            Adds copyright info etc. to a graphics file
# make.SSinits          Finds initial estimates for the state vector in a
#                       dynamic linear model. This is from the software 
#                       provided with Chandler & Scott (2011).
# NApatterns            Returns a vector coding for the patterns of 
#                       missingness in an ensemble of time series. 
# plot.ImportanceWts    Pairs plot of posterior parameter samples, 
#                       with points shaded / coloured according to 
#                       importance sampling weights. 
# PlotEnsTS             Plots annual observed time series, together
#                       with those from each member of an ensemble. 
# PostPred.VidFrames    Produces video frames for an animation of 
#                       samples from a postprocessed ensemble. 
# PostPredSample        Draws samples from the posterior predictive
#                       distribution of observable time series using
#                       a fitted dynamic linear model. 
# print.dlmMLE          print method for dlmMLE objects.
# print.summary.dlmMLE  And for summary.dlmMLE objects. 
# SampleObs             Draws samples of "observed" time series in 
#                       a dlm, conditional on corresponding samples
#                       of parameters and states.
# SampleStates          Draws samples of dlm state vectors, conditional
#                       on corresponding samples of parameter vectors. 
# SLLT.IniPar           Finds initial parameter estimates to initialise
#                       maximum likelihood fitting of a smooth local linear
#                       rend model. 
# SLLT.modeldef,        Define and implement a smooth local linear trend model
# SLLTSmooth            for a single time series.
# SmoothPlot            To plot the smoothing results from a dlm, 
#                       together with the ensembl data used to produce
#                       them. 
# summary.dlmMLE        Summarises the maximum likelihood estimates from
#                       the results of a call to dlmMLE or equivalent
# which.diffuse         Identifies which components of the state vector
#                       require diffuse initialisation, in a dynamic linear
#                       model. From the software for Chandler & Scott (2011). 
#
######################################################################
######################################################################
######################################################################
covxfrm <- function(covmat,A) {
  #
  #   To compute the covariance matrix of A %*% Y, where the 
  #   covariance matrix of Y is given. Arguments:
  #
  #   covmat  Covariance matrix of Y
  #   A       The matrix A
  #
  A %*% covmat %*% t(A) 
}
######################################################################
num.deriv <- function(f, x0, xlo=-Inf, xhi=Inf, steps=NULL,
                      shrink=1.4, maxit=10, relerr=1e-6,
                      QuitEarly=NULL, ...) {
  ######################################################################
  #   Numerical differentiation. Arguments:
  #
  #	f	      The function to be differentiated
  #	x0	    The point at which to evaluate the gradient
  #	lb	    Vector of lower bounds on elements of x (either
  #		      a scalar, or a vector the same length as x0)
  #	ub	    Vector of upper bounds
  #	steps   Vector of largest step sizes to consider when
  #		      evaluating the gradient. Defaults to the larger
  #		      of abs(x0)/10^4 and 0.001.
  #	shrink	Amount by which to shrink the stepsize at each
  #		      iteration, to obtain an improved estimate
  # maxit   Maximum number of iterations
  # relerr  Relative error required
  # QuitEarly	Either NULL (the default) or a list containing 
  #				  two elements, in which case it defines conditions 
  #				  for the routine to exit early if a "surprising" value of 
  #         f is found. In the latter case, the named elements of 
  #				  the list should be "reference" (a numeric value against which to 
  #				  compare the function evaluations in this routine) and 
  #         "greater" (a logical scalar indicating whether to quit if 
  #         a value f>reference is found - if FALSE, the function will quit
  #         if a value f<reference is found).
  #	...	    Any other arguments to f
  #
  #   This uses Ridders' implementation of Richardson's method - 
  #   see Press et al. (1992), Section5.7 for example. The 
  #   implementation here is basically that given by Press et al., 
  #   translated from FORTRAN to R and adapted to handle boundaries 
  #   (in which case a right- or left-hand derivative is evaluated 
  #   as appropriate) and vector arguments. 1-sided derivatives are 
  #   calculated by reflecting the function about the boundary, which 
  #	  is equivalent to using function evaluations on one side of the
  #	  boundary only - and correspondingly less accurate. For 1-sided 
  #	  derivatives it may therefore be necessary to decrease the 
  #	  default value of steps to obtain the required accuracy.
  #
  #	The value of the function is a list containing three elements,
  #	as follows:
  #
  #	gradient	vector of length length(x0), containing
  #				    the calculated derivatives
  #	error		  Estimate of the errors in the computed
  #				    gradients
  # details		Data frame containing details of the calculations:
  #      			column 1 is the element being computed, column 2 is 
  #				    the iteration number, column 3 the value of the 
  #				    function f and the remaining columns are the elements 
  #				    of the function argument.
  #	
  ######################################################################
  np <- length(x0); d <- numeric(np)
  nb <- length(xlo)
  if (nb > 1 & nb != np) stop("xlo and x0 must be the same length")
  nb <- length(xhi)
  if (nb > 1 & nb != np) stop("xhi and x0 must be the same length")
  if (is.null(steps)) {
    h <- pmax(abs(x0)*1e-4,1e-3) 
  } else {
    h <- steps
  }
  right <- (x0 == xlo); left <- (x0 == xhi)
  f0 <- f(x0,...)
  details <- data.frame(matrix(c(0, 0, f0, x0), nrow=1))
  if (!is.null(names(x0))) {
    names(details) <- c("Element", "Iteration", "Value", names(x0))
  } else {
    names(details) <- 
      c("Element", "Iteration", "Value", paste("x[", 1:np, "]", sep=""))
  }
  if (!is.null(QuitEarly)) {
    if (!isTRUE(all(c("reference", "greater") %in% names(QuitEarly)))) {
      stop(paste("If non-NULL, QuitEarly must be a list containing",
                 "elements\n'reference' and 'greater'", sep=" "))
    }
  }
  shr2 <- shrink^2; err <- rep(Inf,np)
  delta <- .Machine$double.eps
  #
  # Set initial step sizes for each parameter. 
  #
  hh <- pmin(h,(x0-xlo-delta),(xhi-delta-x0))
  hh[left] <- pmin(h[left],(xhi-xlo-delta)[left])
  hh[right] <- pmin(h[right],(xhi-xlo-delta)[right])
  tight <- right & left; d[tight] <- NA
  wanted.els <- (1:np)[!tight]
  Done <- FALSE
  for (p in wanted.els) {
    if (Done) break # In case we bailed out during the last iteration
    a <- matrix(nrow=maxit,ncol=maxit)
    i <- 1; cnvrge <- FALSE
    while (i <= maxit & !cnvrge & !Done) {
      xx <- x0
      if (!right[p]) {
        xlo <- x0[p] - hh[p]; xx[p] <- xlo; f1 <- f(xx,...)
        details <- rbind(details, c(p, i, f1, xx))
        if (!is.null(QuitEarly)) Done <- QuitNow(f1, QuitEarly)
        if (Done) break
        #
        #	For left derivatives, extrapolate symmetrically beyond the
        #	upper limit
        #
        if (left[p]) f2 <- (2*f0) - f1
      }
      if (!left[p]) {
        xhi <- x0[p] + hh[p]; xx[p] <- xhi; f2 <- f(xx,...)
        details <- rbind(details, c(p, i, f2, xx))
        if (!is.null(QuitEarly)) Done <- QuitNow(f2, QuitEarly)
        if (Done) break
        #
        #	And for right derivatives, extrapolate beyond the lower
        #	limit
        #
        if (right[p]) f1 <- (2*f0) - f2
      }
      a[1,i] <- (f2 - f1) / (2*hh[p])
      #
      #	This is the j-loop in Numerical Recipes. Although looping
      #	is normally to be avoided in R, experimentation suggests
      #	that in this case it's actally quicker than the corresponding
      #	vector-based calculation - partly do do with the fact that 
      #	in that case you have do do things like identify the 
      #	location of the smallest value in the column etc., I guess
      #
      if (i > 1) {
        fac <- shr2
        for (j in 2:i) {
          a[j,i] <- (a[j-1,i]*fac-a[j-1,i-1]) / (fac-1)
          fac <- shr2*fac
          errt <- max(abs(a[j,i]-a[j-1,i]),abs(a[j,i]-a[j-1,i-1]))
          if (is.na(errt)) errt <- Inf
          if (errt < err[p]) { 
            err[p] <- errt
            d[p] <- a[j,i]
          }
        }
        cnvrge <- (abs(a[i,i]-a[i-1,i-1]) > 2*err[p])
        if (is.na(cnvrge)) cnvrge <- FALSE
        if (abs(d[p]) > relerr) {
          cnvrge <- cnvrge & (abs(err[p]) < abs(relerr * d[p]))
        }
      }
      hh[p] <- hh[p] / shrink; i <- i+1
    }
  }
  list(gradient=d, error=err, details=details)
}
######################################################################
num.hess <- function(f, x0, xlo=-Inf, xhi=Inf, steps=NULL, shrink=1.4,
					maxit=10, relerr=1e-6, QuitEarly=NULL, ...) {
######################################################################
#       Numerical evaluation of a Hessian. Arguments:
#
#	  f	    The function to be differentiated
#	  x0	  The point at which to evaluate the Hessian
#	  xlo	  Vector of lower bounds on elements of x (either
#		      a scalar, or a vector the same length as x0)
#	  xhi	  Vector of upper bounds
#	  steps	Vector of largest step sizes to consider for 
#         finite-difference approximations. Default is
#         approximately the square root of than that in 
#		      num.deriv, because the denominator in each 
#		      approximation is steps^2)
#	  shrink	Amount by which to shrink the stepsize at each
#		      iteration, to obtain an improved estimate
#   maxit Maximum number of iterations
#   relerr  Relative error required
#   QuitEarly	Either NULL (the default) or a list containing 
#				  two elements, in which case it defines conditions 
#				  for the routine to exit early if a "surprising" value of 
#         f is found. In the latter case, the named elements of 
#				  the list should be "reference" (a numeric value against which to 
#				  compare the function evaluations in this routine) and 
#         "greater" (a logical scalar indicating whether to quit if 
#         a value f>reference is found - if FALSE, the function will quit
#         if a value f<reference is found).
#	...	Any other arguments to f
#
# This uses polynomial extrapolation to the limit, in exactly
# the same way as num.deriv(). The result is regarded is the
# limit of a sequence of finite-difference approximations, 
# evaluated at a sequence of points e1, e2, e3, ... 
# which represent progressively smaller multiples of the 
# initial value of steps. The finite difference approximation 
#	actually used is, for off-diagonal elements,
#
#	{f(x1+h1,x2+h2)+f(x1-h1,x2-h2)-[f(x1+h1,x2-h2)+f(x1-h1,x2+h2)]}
#       _______________________________________________________________
#				4*h1*h2
#
#	or (f1 + f2 - f3 - f4) / (4*h1*h2), where h1 and h2 are 
#	multiples of the e's. For diagonal elements the approximation
#	is
#
#	{f(x+h) - f(x) - f(x) + f(x-h)} / h
#
#	- NB the use of h rather than 4h is *correct* here - repeated
#   application of the first-order finite difference calculation
#   yields a step size of 2h, which has been halved in this case.
#
#	The value of the function is a list containing two elements,
#	as follows:
#
#	hessian	matrix containing the calculated Hessian. For
#         all variables where the lower and upper bounds 
#			    xlo and xhi are equal, a value of zero is returned. 
#	error   Estimate of the errors in the computed hessian
# details	Data frame containing details of the calculations:
#      		columns 1 and 2 are the indices of the element being
#				  computed, column 3 is the iteration number, column 4 the
#         value of the function f and the remaining columns
#				  are the elements of the function argument.
#	
######################################################################
 np <- length(x0); d <- matrix(rep(0,np*np),nrow=np,ncol=np)
 nb <- length(xlo)
 if (nb > 1 & nb != np) stop("xlo and x0 must be the same length")
 nb <- length(xhi)
 if (nb > 1 & nb != np) stop("xhi and x0 must be the same length")
 if (is.null(steps)) {
  h <- pmax(abs(x0)*1e-2,1e-2) 
 } else {
  h <- steps
 }
 right <- (x0 == xlo); left <- (x0 == xhi)
 tight <- right & left
 if (any((right | left)[!tight])) {
   stop("I don't know how to evaluate Hessian on a boundary")
 }
 if (!is.null(QuitEarly)) {
   if (!isTRUE(all(c("reference", "greater") %in% names(QuitEarly)))) {
     stop(paste("If non-NULL, QuitEarly must be a list containing",
                "elements\n'reference' and 'greater'", sep=" "))
   }
 }
 shr2 <- shrink^2; err <- matrix(rep(Inf,np*np),nrow=np)
 delta <- .Machine$double.eps
 oppsign <- c(1,-1)
 f0 <- f(x0,...)
 details <- data.frame(matrix(c(0, 0, 0, f0, x0), nrow=1))
 if (!is.null(names(x0))) {
   names(details) <- c("Idx1", "Idx2", "Iteration", "Value", names(x0))
 } else {
   names(details) <- 
       c("Idx1", "Idx2", "Iteration", "Value", paste("x[", 1:np, "]", sep=""))
 }
 Done <- FALSE
 for (p1 in (1:np)[!tight]) {
   if (Done) break # In case we bailed out during the last iteration
   for (p2 in (1:p1)[!tight[1:p1]]) {
     if (Done) break
     #
     # Reset initial step sizes for each parameter. No need to set 
     # tight values to NA here, because d contains NA on initialisation
     #
     hh <- pmin(h,(x0-xlo-delta),(xhi-delta-x0))
     a <- matrix(nrow=maxit,ncol=maxit); els <- c(p1,p2)
     i <- 1; cnvrge <- FALSE
     while (i <= maxit & !cnvrge & !Done) {
       xx <- x0
       #
       # Different finite difference approximations required on and off
       # the diagonal
       #
       if (p1 == p2) {
         xnew <- x0[p1] - hh[p1]; xx[p1] <- xnew; f1 <- f(xx,...)
         details <- rbind(details, c(p1, p2, i, f1, xx))
         if (!is.null(QuitEarly)) Done <- QuitNow(f1, QuitEarly)
         if (Done) break
         xnew <- x0[p1] + hh[p1]; xx[p1] <- xnew; f2 <- f(xx,...)
         details <- rbind(details, c(p1, p2, i, f2, xx))
         if (!is.null(QuitEarly)) Done <- QuitNow(f2, QuitEarly)
         if (Done) break
         a[1,i] <- ((f1 - f0) - (f0 - f2)) / (hh[p1]*hh[p1])
       } else {
         xnew <- x0[els] + hh[els]; xx[els] <- xnew; f1 <- f(xx,...)
         details <- rbind(details, c(p1, p2, i, f1, xx))
         if (!is.null(QuitEarly)) Done <- QuitNow(f1, QuitEarly)
         if (Done) break
         xnew <- x0[els] - hh[els]; xx[els] <- xnew; f2 <- f(xx,...)
         details <- rbind(details, c(p1, p2, i, f2, xx))
         if (!is.null(QuitEarly)) Done <- QuitNow(f2, QuitEarly)
         if (Done) break
         xnew <- x0[els] + (hh[els]*oppsign); xx[els] <- xnew; f3 <- f(xx,...)
         details <- rbind(details, c(p1, p2, i, f3, xx))
         if (!is.null(QuitEarly)) Done <- QuitNow(f3, QuitEarly)
         if (Done) break
         xnew <- x0[els] - (hh[els]*oppsign); xx[els] <- xnew; f4 <- f(xx,...)
         details <- rbind(details, c(p1, p2, i, f4, xx))
         if (!is.null(QuitEarly)) Done <- QuitNow(f4, QuitEarly)
         if (Done) break
         a[1,i] <- (f1 + f2 - f3 - f4) / (4*hh[p1]*hh[p2])
       }
       #
       #	This is the j-loop in Numerical Recipes. Although looping
       #	is normally to be avoided in R, experimentation suggests
       #	that in this case it's actally quicker than the corresponding
       #	vector-based calculation - partly do do with the fact that 
       #	in that case you have do do things like identify the 
       #	location of the smallest value in the column etc., I guess
       #
       if (i > 1) {
         fac <- shr2 
         for (j in 2:i) {
           a[j,i] <- (a[j-1,i]*fac-a[j-1,i-1]) / (fac-1)
           fac <- shr2*fac
           errt <- max(abs(a[j,i]-a[j-1,i]),abs(a[j,i]-a[j-1,i-1]))
           if (is.na(errt)) errt <- Inf
           if (errt < err[p1,p2]) { 
             err[p1,p2] <- errt
             d[p1,p2] <- a[j,i]
           }
         }
         cnvrge <- (abs(a[i,i]-a[i-1,i-1]) > 2*err[p1,p2]) 
         if (is.na(cnvrge)) cnvrge <- FALSE
         if (abs(d[p1,p2]) > relerr) {
           cnvrge <- cnvrge & (abs(err[p1,p2]) < abs(relerr * d[p1,p2]))
         }
       }
       hh[els] <- hh[els] / shrink; i <- i+1
     }
     d[p2,p1] <- d[p1,p2]
     err[p2,p1] <- err[p1,p2]
   }
 }
 list(hessian=d,error=err, details=details)
}
######################################################################
QuitNow <- function(x, QuitEarly) { 
  # 
  # Utility function for use in num.deriv and num.hess
  #
  SignMult <- 2*as.numeric(QuitEarly$greater)-1 # 1 or -1
  isTRUE(SignMult * (x-QuitEarly$reference) > 0)
}
######################################################################
CompactEns <- function(Data, Groups, impute=FALSE) {
  #
  #  Produces a compact representation of an ensemble of time series,
  #  in which multiple members from the same source are replaced with
  #  a single series which is the mean of those members. Arguments:
  #
  #  Data   Matrix or data frame containing observed time series in its
  #         first column and ensemble members in the remainder.
  #         NB column 1 is ignored: it's part of the argument
  #         for consistency with other routines in the TimSPEC
  #         package. May have an attribute RunsPerTS (see the value
  #         of the function below), if the ensemble has already been
  #         partially compacted.
  #  Groups A vector (numeric, character or factor) specifying the 
  #         grouping structure in the ensemble, of length 
  #         ncol(Data)-1 corresponding to the ensemble-derived
  #         series in Data. Series with the same value of Groups
  #         are considered to come from the same source (e.g. the
  #         same GCM but under different initial conditions).
  #  impute Controls whether to attempt to impute missing values
  #         using EnsImpute.
  #
  #  The function returns a matrix or data frame analogous to Data, but 
  #  containing the observed series in column 1 and the group mean 
  #  series in subsequent columns. The result has an attribute 
  #  RunsPerTS, which is a vector containing the numbers of series
  #  contributing to each of the group means. If Data itself 
  #  had such an attribute, then the RunsPerTS attribute of the 
  #  result is obtained by summing the corresponding elements and 
  #  the "group mean" series are weighted accordingly; otherwise 
  #  it's just a count of the numbers of series being averaged. 
  #
  ######################################################################
  y <- Data[,-1] # Initial data matrix: each column is a member
  NGroups <- length(unique(Groups))
  z <- y[,1:NGroups] # Same class as y but with a column per group
  #
  #   Can have time points with no ensemble members: exclude
  #   these from calculations
  #
  NoData <- apply(y, MARGIN=1, FUN=function(x) all(is.na(x)))
  #
  #   Do imputation if requested. 
  #
  if (impute) y <- EnsImpute(y)
  #
  #  Set number of runs contributing to each column to 1 if necessary
  #
  RunsPerTS <- attr(Data, "RunsPerTS")
  if (is.null(RunsPerTS)) RunsPerTS <- rep(1, ncol(y))
  if (length(RunsPerTS) != ncol(y)) {
    stop(paste("RunsPerTS attribute of Data should have", ncol(y), "elements"))
  }
  #
  #   Now we have the data from which the group mean series can be
  #   calculated. Here are functions to calculate numerators and
  #   denominators for group means in a single row. 
  #
  WtSum <- function(x, RunsPerTS, Groups) { # Weighted sums by group
    n <- RunsPerTS; n[is.na(x)] <- 0
    xx <- x; xx[is.na(x)] <- 0 # sum(xx)/sum(n) will give correct means
    z <- tapply(xx*n, INDEX=Groups, FUN=sum)
    z
  } 
  Numers <- apply(y[!NoData,], MARGIN=1, FUN=WtSum, 
                  RunsPerTS=RunsPerTS, Groups=Groups)
  NSum <- function(x, RunsPerTS, Groups) { # No. of runs per group
    n <- RunsPerTS; n[is.na(x)] <- 0
    z <- tapply(n, INDEX=Groups, FUN=sum)
    z
  } 
  Denoms <- apply(y[!NoData,], MARGIN=1, FUN=NSum, 
                  RunsPerTS=RunsPerTS, Groups=Groups)
  #
  #   Check that there are no rows with incomplete but non-empty groups
  #
  DenomCheck <- apply(Denoms, MARGIN=1, FUN=function(x) {
    max(x[x>0], na.rm=TRUE)-min(x[x>0], na.rm=TRUE) < 
      10*.Machine$double.neg.eps
  })
  if (!all(DenomCheck)) {
    stop("Some time points have non-empty but incomplete groups")
  }
  z[NoData,] <- NA
  TotRuns <- apply(Denoms, MARGIN=1, FUN=max, na.rm=TRUE)
  z[!NoData,] <- t(Numers / Denoms)
  z <- cbind(Data[,1], z) # Add the observed series back
  attr(z, "RunsPerTS") <- TotRuns
  colnames(z) <- c(colnames(Data)[1], row.names(Numers))
  z
}
######################################################################
EnsImpute <- function(Y, method="lm") {
  #
  #  To impute missing values in an ensemble of time series. 
  #  Missing values are replaced by the fitted values from an 
  #  additive ANOVA model with effects corresponding to ensemble 
  #  member and to the time variable (treated as a factor). 
  #  The argument Y is a matrix or data frame, with each column 
  #  corresponding to an ensemble member and with each row
  #  corresponding to a time point. The function returns an
  #  object like Y, with the NAs replaced as described. 
  #
  #  With method="lm" (the default), the ANOVA model is fitted
  #  using least squares. With method="approx", an approximate
  #  fit is obtained using formulae for a complete 2-way layout:
  #  this is much cheaper than the least-squares approach and
  #  should give an accurate approximation if the proportion of 
  #  missing values is very small. 
  #
  ######################################################################
  z <- Y # Storage for final result
  #
  #   Can have time points with no ensemble members: exclude
  #   these from calculations
  #
  NoData <- apply(Y, MARGIN=1, FUN=function(x) all(is.na(x)))
  y <- Y[!NoData,]
  #
  #   Do imputation if requested. The approximate method is 
  #   based on the algebraic computation for a 2-way additive 
  #   ANOVA with complete data which is 
  #
  #        Y.. + (Yi.-Y..) + (Y.j-Y..) = Yi. + (Y.j-Y..),
  #
  #   but estimating each term using the most complete set of
  #   data possible.
  #
  if (method=="approx") {
    MbrEffect <- colMeans(y, na.rm=TRUE) # Ensemble member effects
    TimeEffect <- as.numeric(scale(rowMeans(y, na.rm=TRUE), scale=FALSE))
    Fitted <- outer(TimeEffect, MbrEffect, "+") 
  } else if (method=="lm") {
    Row <- as.factor(row(y)); Col <- as.factor(col(y))
    Model <- lm(as.numeric(as.matrix(y)) ~ Row + Col)
    Fitted <- matrix(predict(Model, newdata=data.frame(Row=Row, Col=Col)),
                     nrow=nrow(y))
  } else {
    stop("'method' must be either 'approx' or 'lm'")
  }
  y[is.na(y)] <- Fitted[is.na(y)]
  z[!NoData,] <- y
  z
}
######################################################################
NApatterns <- function(Y) {
  #
  #   To identify the patterns of missingness in rows of a matrix
  #   or data frame (Y, the input). The function returns a character
  #   vector with an element for each row of y, encoding the 
  #   positions of the NA values in that row. For example, if 
  #   the values in columns 1 and 3 are NA then the corresponding
  #   element of the result will be "1,3". 
  #
  apply(Y, MARGIN=1, FUN=function(x) paste(which(is.na(x)), collapse=","))
}
######################################################################
NA.AllOr0 <- function(Y) {
  #
  #   To check whether the rows of a data frame or matrix are either
  #   all missing or all non-missing. Returns TRUE if so, and FALSE if 
  #   at least one row is partially missing.  
  #
  isTRUE(
    all(apply(Y, MARGIN=1, FUN=function(x) { 
      sum(is.na(x)) %in% c(0,length(x)) } ))
  )
}
######################################################################
EnsMean <- function(Data, Groups=NULL, na.option=NULL) {
  #
  #  Calculates an ensemble mean series from its members. Arguments: 
  #
  #  Data   Data frame containing observed time series in its
  #         first column and ensemble members in the remainder.
  #         NB column 1 is ignored: it's part of the argument
  #         for consistency with other routines in the TimSPEC
  #         package.
  #  Groups A vector or matrix specifying any structure in the
  #         ensemble, with a row for each ensemble member and
  #         (in the matrix case) columns representing the 
  #         codes of any grouping variables defining the structure.
  #         This corresponds to the Groups argument in, for
  #         example, EnsEBMtrendSmooth or EnsEBM2waytrendSmooth.
  #  na.option  Controls the treatment of missing ensemble members.
  #         Options are "na.rm" (remove NA cases prior to the 
  #         calculations), "na.impute" (see below) or "na.pass"
  #         (propagate the NAs through the calculations). 
  #         mean? The default behaviour corresponds to "na.rm"
  #         if Groups is NULL, and "na.impute" otherwise.
  #
  #  The function returns a vector containing the computed series.
  #  If Groups is NULL (the default) then this series is the 
  #  overall mean of the ensemble members (subject to the 
  #  treatment of missing values according to the value of 
  #  na.option). Otherwise, the series corresponds to the 
  #  sequence of intercepts in ANOVA models fitted to the 
  #  ensemble members at each time point, with factors defined
  #  by the column(s) of Groups and using sum-to-zero 
  #  constraints to ensure that the intercept has the desired
  #  interpretation as a mean of the simulators or simulator 
  #  combinations.
  #
  #  If na.option is "na.rm" then the missing value patterns
  #  are processed separately at each time point. This can 
  #  lead to discontinuities in the calculated ensemble mean
  #  series if the ensemble members cover different periods, 
  #  especially in structured ensembles where one simulator
  #  has been run for a slightly different time period from
  #  the remainder. In this case, the "na.impute" option may 
  #  be preferable. With this option, means are computed
  #  after imputing any missing values using an additive ANOVA
  #  model with effects corresponding to ensemble member
  #  and to the time variable (treated as a factor). 
  #
  ######################################################################
  y <- as.matrix(Data[,-1]) # Initial data matrix: each column is a member
  z <- rep(NA, nrow(y)) # Storage for final result
  #
  #   Can have time points with no ensemble members: exclude
  #   these from calculations
  #
  NoData <- apply(y, MARGIN=1, FUN=function(x) all(is.na(x)))
  y <- y[!NoData,]
  #
  #   Do imputation if requested. 
  #
  if ( (is.null(Groups) & "na.impute" %in% na.option) |
       (!is.null(Groups) & ("na.impute" %in% na.option | is.null(na.option))) ) {
    y <- EnsImpute(y)
  }
  #
  #   Now we have the data from which the means can be calculated
  #
  if (is.null(Groups)) { # Unstructured case: "na.pass" will lead to NAs
    na.rm <- (is.null(na.option) | "na.rm" %in% na.option)
    z[!NoData] <- rowMeans(y, na.rm=na.rm)
  } else {
    #
    #   Now deal with other missing value patterns. The precise
    #   treatment depends on na.option, but in the first instance
    #   it's helpful to record the patterns of missingness at
    #   each time point. Then loop over the patterns needed
    #   (which, at this stage, are just the "no missing members"
    #   pattern "" if na.option is either na.impute or na.pass,
    #   and all the available patterns if it's na.rm); in each
    #   case construct the design matrix for the ANOVA model
    #   fitted to the required subset of data, and then calculate
    #   the intercept directly (avoids use of lm() which produces
    #   a whole load of other stuff that isn't needed).
    #
    NAs <- NApatterns(y)
    PatternsToProcess <- ""
    if ("na.rm" %in% na.option) PatternsToProcess <- unique(NAs)
    Intercepts <- rep(NA, nrow(y))
    for (i in 1:length(PatternsToProcess)) {
      Times <- (NAs==PatternsToProcess[i]) # Time points with this pattern
      Members <- # Members contributing to estimates for this pattern
        !( (1:ncol(y)) %in% unlist(strsplit(PatternsToProcess[i], ",")) )
      GrpFac <- apply(as.matrix(Groups)[Members,,drop=FALSE], 
                      MARGIN=2, FUN=as.factor, simplify=FALSE)
      names(GrpFac) <- paste("V",1:length(GrpFac),sep="")
      GrpFac <- list2DF(GrpFac)
      OldOpt <- options()$contrasts
      options(contrasts = c("contr.sum", "contr.poly"))
      X <- model.matrix(~ ., GrpFac) 
      options(contrasts = OldOpt)
      XXInv <- solve(crossprod(X))[1,,drop=FALSE] # Just need first row
      XY <- crossprod(X, t(y[Times, Members])) # NB standard formula is
                                               # X'Y but that's when Y 
                                               # is a column vector: here
                                               # the vectors are rows of y
      Intercepts[Times] <- XXInv %*% XY
    }
    z[!NoData] <- Intercepts
  }
  z
}
######################################################################
PlotEnsTS <- function(Data, Colours, Groups=NULL, na.option=NULL, EnsTransp=1, 
                      Types=c(1,1), Units=expression(degree*"C"), 
                      plot=TRUE, add=FALSE, ...) {
  #
  #  Plots observed annual time series, individual ensemble members and 
  #  ensemble mean. Arguments: 
  #
  #  Data   Data frame containing time in column 1, observed series
  #         in column 2 and ensemble members in the remaining
  #         columns. 
  #  Groups See documentation for EnsMean
  #  na.option  Ditto
  #  Colours Vector containing the colour specifications for 
  #         the observations and ensemble members respectively.
  #         If of length 2 then the two colours will be used
  #         to represent the observations and ensemble 
  #         members, respectively. If of length ncol(Data)
  #         then the first element will be used for the 
  #         observations, the second for the ensemble mean
  #         and the remainder for the individual ensemble
  #         members. 
  #  EnsTransp Transparency (alpha) value to use when plotting
  #         the ensemble members.
  #  Types  Vector containing the line types to use for 
  #         observations and ensemble members respectively. 
  #         Specification as for Cols, defaults to c(1,1)
  #         i.e. solid lines throughout. 
  #  Units  Length-1 character vector, used to label the
  #         vertical axis of the plot. Default value 
  #         is appropriate for plots of temperature in 
  #         degrees Celsius.
  #  plot   If TRUE (the default) then plot the results;
  #         otherwise just set up the plot region but don't
  #         actually plot anything.
  #  add    If TRUE, add lines to an existing plot; if
  #         FALSE (the default) then create a new plot. 
  #  ...    Other arguments to matplot(). 
  #
  #  The function returns a NULL value: its purpose is to produce
  #  a plot, not to perform calculations for subsequent use. 
  #
  ncols <- length(Colours)
  if (ncols==2) {
    ObsCol <- Colours[1]
    EnsMeanCol <- Colours[2]
    EnsCols <- Colours[2]
  } else if (ncols == ncol(Data)) {
    ObsCol <- Colours[1]
    EnsMeanCol <- Colours[2]
    EnsCols <- Colours[-(1:2)]
  } else {
    stop("Length of 'Colours' doesn't match dimensions of 'Data'")
  }
  if (EnsTransp < 1) {  # Set transparency for ensemble members  
    EnsCols <- rgb(t(col2rgb(EnsCols)) / 255, alpha=EnsTransp)
  }
  ntypes <- length(Types)
  if (ntypes==2) {
    ObsType <- Types[1]
    EnsMeanType <- Types[2]
    EnsTypes <- Types[2]
  } else if (ntypes == ncol(Data)) {
    ObsType <- Types[1]
    EnsMeanType <- Types[2]
    EnsTypes <- Types[-(1:2)]
  } else {
    stop("Length of 'Types' doesn't match dimensions of 'Data'")
  }
  if (!add) matplot(Data[,1], Data[,-1], 
                    type="n", xlab="Year", ylab=Units, ...)
  if (plot) {
    MeanLine <- EnsMean(Data[,-1, drop=FALSE], Groups, na.option)
    matlines(Data[,1], Data[,-(1:2)], col=EnsCols, lty=EnsTypes)
    lines(Data[,1], MeanLine, col=EnsMeanCol, lty=EnsMeanType, lwd=5)
    lines(Data[,1:2], col=ObsCol, lty=ObsType, lwd=5)
  }
  invisible(NULL)
}
######################################################################
SmoothPlot <- 
  function(Data, Smooth=NULL, Samples=NULL, AntiLog=FALSE, NonNeg=FALSE,
           PlotMu0=TRUE, PlotConsensus=FALSE, Weights=NULL, DatColours, 
           DatTypes=c(1,1), PredColours, alpha=c(0.6,0.1), 
           SmoothPIs=FALSE, Units=expression(degree*"C"), 
           plot.title="", Legend=TRUE, LegPos="topleft", ...) {
  #
  #   To plot the smoothing results from a dlm, together with
  #   the data used to produce them. Arguments:
  #
  #   Data    Data frame containing time in column 1, observed
  #           series in column 2 and ensemble members in the 
  #           remaining columns. 
  #   Smooth  A list containing at least the named components
  #           "Model" (an object of class dlm) and "Smooth"
  #           (the result of Kalman Smoothing a dataset using 
  #           the dlm stored in Model, via the dlmSmooth command.
  #           Typically produced via one of the *Smooth()
  #           functions in this script. Ignored (with a 
  #           warning) if Samples is non-NULL. NB though: 
  #           if used, it is assumed that the "observed"
  #           time series forms the first element of each
  #           data vector Y[t].
  #   Samples List containing at least the named elements
  #           States (an array containing posterior 
  #           samples of the state vector, produced using the
  #           SampleStates() function) and Obs (a matrix
  #           containing samples from the posterior predictive
  #           distribution of the observed series, produced
  #           using SampleObs()). If non-NULL, this is 
  #           used to compute prediction intervals and
  #           posterior means of relevant quantities. 
  #   PlotMu0 Logical scalar controlling whether or not
  #           to plot the estimated trend in the observations. 
  #   PlotConsensus  Controls whether to add the ensemble 
  #           consensus to the plot.
  #   AntiLog If TRUE, all results will be exponentiated 
  #           (with appropriate adjustments for the posterior
  #           mean) prior to plotting. This allows for situations
  #           where the DLM has been fitted to logged data
  #           (e.g. to ensure that postprocessed precipitation
  #           samples are guaranteed to be non-negative).
  #   NonNeg  If FALSE (the default), prediction intervals will
  #           be computed in the usual way as mean +/- 
  #           1.96 standard deviations. If TRUE, the quantiles
  #           of a gamma distribution will be used instead to
  #           ensure that the 
  #   Weights Optional vector of importance weights 
  #           associated with the samples in Samples.
  #           Ignored if Samples itself contains an
  #           element named Weights (as, for example, 
  #           in a call to PostPredSample()).
  #   DatColours, Vectors of colours and line types used to 
  #   DatTypes    plot the ensemble members: for details, 
  #               see header to PlotEnsTS(). 
  #   PredColours Vector of two colours: the first is the
  #               base colour for plotting estimates of
  #               the observed trend mu[0], and the 
  #               second for plotting the ensemble consensus.
  #   alpha   Vector of two transparency values (between 0 & 1)
  #           used for plotting estimates and prediction
  #           intervals respectively.
  #   SmoothPIs   If TRUE and if Samples is provided, the
  #           prediction intervals on the plot will be 
  #           smoothed using loess. Otherwise (the default)
  #           the raw quantiles from the sampled series
  #           will be plotted. 
  #   Units   Length-1 character vector, used to label the
  #           vertical axis of the plot. Default value 
  #           is appropriate for plots of temperature in 
  #           degrees Celsius.
  #   plot.title  Self-explanatory, hopefully. 
  #   Legend  Logical controlling whether or not a legend is
  #           added to the plot
  #	  LegPos  Keyword controlling position of the legend - see
  #           help for legend
  #   ...     Other arguments to PlotEnsTS(). 
  #
  #   The function returns, invisibly, a data frame containing 
  #   the values of the estimated "real" and "ensemble consensus"
  #   trends, together with the upper and lower 95% limits for 
  #   the real trend. The interpretation of the limits depends
  #   on the inputs: for a "Smooth" object, they're pointwise
  #   confidence bands while for a "Samples" object they are
  #   pointwise posterior predictive intervals.
  #
  if (!is.null(Smooth) & !is.null(Samples)) {
    warning("'Smooth' argument is ignored when 'Samples' is non-NULL")
  }
  if (!is.null(Samples)) { # Set weights to be used if Samples is non-NULL
    WtsToUse <- Samples$Weights$w
    if (is.null(WtsToUse)) WtsToUse <- Weights
  }
  PlotData <- Data
  if (AntiLog) PlotData[,-1] <- exp(PlotData[,-1])
  #
  #   Make transparent colours for estimates and prediction intervals
  #
  pred.cols <- CI.cols <- t(col2rgb(PredColours)) / 255
  pred.cols <- rgb(pred.cols, alpha=alpha[1]) # Transparent versions for preds & CIs
  CI.cols <- rgb(CI.cols, alpha=alpha[2])
  if (is.null(Samples)) { # Take plugin MAP estimates of everything
    Preds <- dlm.ObsPred(Smooth)
    Mu0.hat <- Smooth$Smooth$s[-1, 1]
    if (NonNeg) { # For non-negative quantities, use quantiles of gamma dbn
      Shape <- (Mu0.hat / Preds$SE[,1])^2
      Rate <- Mu0.hat / (Preds$SE[,1]^2)
      CI.Lims <- matrix(NA, nrow=2, ncol=length(Mu0.hat))
      OK <- Mu0.hat > 0 # Trap negative means
      CI.Lims[1,OK] <- qgamma(0.025, shape=Shape[OK], rate=Rate[OK])
      CI.Lims[2,OK] <- qgamma(0.975, shape=Shape[OK], rate=Rate[OK])
      if (!all(OK)) {
        warning("NonNeg is TRUE but not all means are strictly positive")
        CI.Lims[,!OK] <- 0 # To avoid problems in subsequent plotting with polygon
      }
    } else {
      CI.Lims <- rbind(Mu0.hat - (1.96 * Preds$SE[,1]),
                       Mu0.hat + (1.96 * Preds$SE[,1])) 
    }
    if (AntiLog) {
      Mu0.hat <- exp(Mu0.hat + (Preds$SE[,1]^2)/2)
      CI.Lims <- exp(CI.Lims)
    }
    if (PlotConsensus) {
      if (!("Consensus" %in% names(Smooth$Model))) {
        stop("Can't plot ensemble consensus without a 'Consensus' component of Model")
      }
      A <- matrix(Smooth$Model$Consensus, nrow=1)
      MuCons <- Smooth$Smooth$s[-1,] %*% t(A)
      if (AntiLog) {
        state.vars <- 
          simplify2array(dlmSvd2var(Smooth$Smooth$U.S,
                                    Smooth$Smooth$D.S)[-1]) # Covariance matrices
        MuCons.Var <- apply(X=state.vars, MARGIN=3, FUN=covxfrm, A=A)
        #
        #   Next line should be a roughly unbiased estimator of the
        #   true consensus on the original scale
        #
        MuCons <- exp(MuCons - MuCons.Var/2) 
      }
    }
  } else { # Or use posterior samples
    Mu0.hat <- Samples$States[,,1] # Mu0[t] is first element of state vector
    if (AntiLog) Mu0.hat <- exp(Mu0.hat)
    Obs <- Samples$Obs
    if (AntiLog) Obs <- exp(Obs)
    if (PlotConsensus) {
      if (!("Consensus" %in% names(Samples$Model))) {
        stop("Can't plot ensemble consensus without a 'Consensus' component of Model")
      }
      EnsCons.Sampled <- # Calculate individual consensus estimates
        apply(Samples$States, MARGIN=1:2, 
              FUN=function(x) x %*% Samples$Model$Consensus)
      if (AntiLog) EnsCons.Sampled <- exp(EnsCons.Sampled)
    }
    if (is.null(WtsToUse)) {
      Mu0.hat <- apply(Mu0.hat, MARGIN=2, FUN=mean, na.rm=TRUE) 
      CI.Lims <- apply(Obs, MARGIN=1, FUN=quantile, 
                     probs=c(0.025,0.975),na.rm=TRUE)
      if (PlotConsensus) MuCons <- colMeans(EnsCons.Sampled, na.rm=TRUE)
    } else {
      Mu0.hat <-
        apply(Mu0.hat, MARGIN=2, FUN=weighted.mean, 
              w=WtsToUse, na.rm=TRUE) 
      #
      #   Turn off warnings here: if one or two of the weights
      #   are strongly dominant then wtd.quantile will get
      #   very upset and issue a warning for every row of
      #   Samples$Obs. Instead, issue a single warning if
      #   90% of the total weight comes from 10% or fewer
      #   of the samples. 
      #
      owarn <- options()$warn
      options(warn=-1)
      CI.Lims <- apply(Obs, MARGIN=1, FUN=wtd.quantile, 
                       probs=c(0.025,0.975),
                       weights=WtsToUse, normwt=TRUE, na.rm=TRUE)
      options(warn=owarn)
      SortedWts <- sort(WtsToUse, decreasing=TRUE)
      NWts <- length(WtsToUse)
      if (PlotConsensus) {
        MuCons <- apply(EnsCons.Sampled, MARGIN=2, FUN=weighted.mean, 
                        w=WtsToUse, na.rm=TRUE)
      }
    }
    if (SmoothPIs) { # Smooth the prediction intervals if requested
      CI.Lims <- t(apply(CI.Lims, MARGIN=1, 
                       FUN=function(x) loess(x ~ I(1:length(x)))$fitted))
    }
  }
  #
  #   Now plot everything
  #
  PlotEnsTS(PlotData, Colours=DatColours, Types=DatTypes, Units=Units, ...)
  YrLims <- par("usr")[1:2]
  abline(h=pretty(par("usr")[3:4], n=10), col=grey(0.8), lty=2)
  abline(v=pretty(par("usr")[1:2], n=diff(YrLims)/10),col=grey(0.8),lty=2)
  if (PlotConsensus) lines(Data[,1], MuCons, lwd=5, col=pred.cols[2])
  polygon(c(Data[,1], rev(Data[,1])), c(CI.Lims[1,], rev(CI.Lims[2,])),
          border=NA, col=CI.cols[1])
  if (PlotMu0) lines(Data[,1], Mu0.hat, lwd=5, col=pred.cols[1])
  title(main=plot.title)
  #
  #   Set up legend depending on what's being plotted
  #
  if (Legend) {
    leg.cols <- c(DatColours[c(1,2,2)], PredColours[c(1,2)], NA)
    leg.lwd <- c(5,5,1,5,5,NA)
    leg.fill <- c(rep(NA,5), CI.cols[1])
    leg.text <- c("Observations", "Ensemble mean", "Ensemble members", 
                  expression(hat(mu)[0](t)),
                  "Ensemble consensus",
                  expression("95% interval for "*Y[0](t)))
    Wanted <- rep(TRUE, length(leg.cols))
    if (!PlotConsensus) Wanted[5] <- FALSE
    if (!PlotMu0) Wanted[4] <- FALSE
    if (ncol(PlotData)<3) Wanted[c(2,3,5)] <- FALSE
    legend(LegPos, col=leg.cols[Wanted], lwd=leg.lwd[Wanted], 
           fill=leg.fill[Wanted], border=NA, ncol=2, bg="white",
           legend=leg.text[Wanted])
  }
  #
  #   Assemble return value
  #
  z <- data.frame(Time=Data[,1], Mu0=Mu0.hat, LL95=CI.Lims[1,],
                  UL95=CI.Lims[2,])
  if (PlotConsensus) z$Consensus <- MuCons
  invisible(z)
}
######################################################################
GraphAnnot <- function(Copyright="\uA9 UCL 2022",
                       Acknowl="Funded by the UK Climate Resilience Programme") {
  #
  #   To add copyright info etc. to a plot.
  #   Arguments: 
  #
  #   Copyright     Text to be added at bottom left of
  #                 plot, typically a copyright statement. 
  #   Acknowl       Text to be added at bottom right. 
  # 
  #   The function merely adds the required text to an existing
  #   plot, and returns a NULL result. 
  #
  TextSize <- par("din")[2]*0.12
  mtext(Copyright, outer=TRUE, side=1, line=-1, 
        at=0.01, adj=0, padj=0, cex=TextSize)
  mtext(Acknowl, outer=TRUE, side=1, line=-1, 
        at=0.99, adj=1, padj=0, cex=TextSize)
  invisible(NULL)
}
######################################################################
make.SSinits <- function(GG,W,kappa) {
  ######################################################################
  #
  # Initialises a time-invariant state space model with transition 
  # matrix GG and innovation covariance matrix W. Stationary elements 
  # of the state vector are initialised with the stationary mean 
  # which is zero; nonstationary elements are also initialised with 
  # zero. The corresponding covariance matrices are initialised with 
  # the stationary covariance matrix (for the stationary parts) and 
  # kappa times the identity matrix for the nonstationary part. 
  # The covariance matrix for the stationary part is calculated 
  # according to equation (3.3.21) of Harvey (1989). 
  #
  # Value: a list containing components m0 and C0: m0 is the estimate
  # of the state vector prior to seeing any data and C0 is the associated 
  # error covariance matrix.
  #
  ######################################################################
  if (!is.matrix(GG)) stop("GG should be a matrix")
  if (!is.matrix(W)) stop("W should be a matrix")
  p <- nrow(GG)
  if (ncol(GG) != p) stop("GG should be a square matrix")
  if (ncol(W) != p | nrow(W) != p) stop("W should be the same size as GG")
  
  m0 <- rep(0,p)
  C0 <- matrix(0,nrow=p,ncol=p)
  diff.els <- which.diffuse(GG)
  diag(C0)[diff.els] <- kappa
  p <- sum(!diff.els)
  if (p>0) {
    GG.stat <- GG[!diff.els,!diff.els,drop=FALSE]
    W.stat <- W[!diff.els,!diff.els,drop=FALSE]
    GG.GG <- GG.stat %x% GG.stat
    II <- diag(rep(1,p*p))
    Cvec <- solve(II-GG.GG,as.vector(W.stat))
    C0[!diff.els,!diff.els] <- Cvec
  }
  list(m0=m0,C0=(C0+t(C0))/2)	# Ensure C0 is symmetric
}
######################################################################
which.diffuse <- function(GG) {
  ######################################################################
  #
  # Takes the transition matrix GG of a state space model and identifies
  # which elements of the state vector are nonstationary. No checking is
  # done; this is designed to be called from other functions.
  #
  # Value: a logical vector that is TRUE for the nonstationary components
  #        of the state vector and FALSE for the stationary components.
  #
  ######################################################################
  p <- nrow(GG)
  dep.marks <- (GG!=0)
  n.depend <- rowSums(dep.marks)
  no.self.depend <- 1-diag(dep.marks) 
  #
  # The next line creates an ordering of the state elements, 
  # sorted firstly by whether or not they depend directly on
  # their own previous values, and within that by the number
  # of dependencies. Thus we start by checking those elements
  # of the state vector that only depend on their own 
  # previous values, which is easy; as we discover nonstationary
  # components, we can then look to see which other components
  # depend on them because these must also be nonstationary. 
  # As we go through, whenever we get to a component that hasn't
  # been marked as nonstationary we need to check whether the 
  # subsystem upon which it depends is stationary; if all else
  # fails, this is done by examining the largest eigenvalue of
  # that subsystem. Hopefully however, the reordering will 
  # reduce the number of required eigendecompositions to the 
  # bare minimum. 
  #
  state.code <- (p*no.self.depend)+n.depend
  direct.dep <- state.code == 1
  state.order <- order(state.code)
  z <- rep(FALSE,p); checked <- z
  for (s in (1:p)[direct.dep[state.order]]) {
    z[s] <- abs(diag(GG)[s]) >= 1
    checked[s] <- TRUE
    if (z[s]) {
      z[dep.marks[,s]] <- TRUE
      checked[dep.marks[,s]] <- TRUE
    }
  }
  for (s in (1:p)[!direct.dep[state.order]]) {
    if (!checked[s]) {
      subsys.deps <- dep.marks[s,]
      subsys.deps[s] <- TRUE
      subsys.mat <- GG[subsys.deps,subsys.deps]
      lambda <- max(abs(eigen(subsys.mat,symmetric=FALSE,only.values=TRUE)$values))
      z[s] <- lambda >= 1
    }
    checked[s] <- TRUE
    if (z[s]) {
      z[dep.marks[,s]] <- TRUE
      checked[dep.marks[,s]] <- TRUE
    }
  }
  z
}
######################################################################
dlm.SafeLL <- function(theta, Y, build, prior.pars=NULL, 
                       BigVal=1e12, debug=FALSE, verbose=FALSE, ...) {
  #
  # "Safe" computation of a (penalised) log-likelihood for a dynamic 
  # linear model. This is a wrapper for the dlmLL routine in the dlm 
  # library, which (a) handles situations where dlmLL returns an error 
  # or an infinite  value (b) is designed for use in calls to nlm() 
  # rather than optim() (c) allows for maximisation of a penalised
  # likelihood which can be regarded as an approximation to a 
  # Bayesian posterior. Arguments:
  #
  # theta   The parameter vector
  # Y       The data
  # build   A function that will construct the dlm object representing
  #         the model for which the log-likelihood is required. Typically
  #         one of the .modeldef functions in this script. 
  # prior.pars  If non-null, a two-column matrix with number of rows
  #         equal to the length of theta: the columns contain
  #         respectively the means and standard deviations of 
  #         independent Gaussian priors on the elements of theta.
  #         If NULL (the default), no prior information is used. 
  # BigVal  Value to return in case dlmLL fails. This should be
  #         much larger than any reasonably expected value of the 
  #         negative log-likelihood, but not so large as to 
  #         cause numerical instability when evaluating gradients
  #         etc. 
  # debug   As in dlmLL
  # verbose If TRUE, print the parameter value and calculated 
  #         (penalised) log-likelihood
  # ...     Other arguments passed to build
  #
  # If prior.pars is NULL then the function returns the result 
  # of dlmLL when that routine returns a finite number, otherwise 
  # it returns BigVal. If prior.pars is non-null, the function
  # returns the negative log posterior (or penalised negative 
  # log-likelihood if you prefer) when dlmLL returns a finite number. 
  # In this case the result also has an "Unpenalised" attribute
  # containing the unpenalised negative log-likelihood: this can 
  # be helpful for comparing (informally) fits of models with
  # different parameters & hence different priors. 
  #
  # The next couple of lines remove any arguments that are not needed
  # by the build() function (which may be passed via the ... argument
  # e.g. when using nlm, and which cause errors). Also ensure that
  # any defaults are passed explicitly. 
  #
  ArgsGiven <- c(as.list(environment()), list(...))
  ArgsNeeded <- formals(build)
  BuildArgs <- c(ArgsGiven[names(ArgsGiven) %in% names(ArgsNeeded)],
                 ArgsNeeded[!(names(ArgsNeeded) %in% names(ArgsGiven))])
  Model <- do.call(build, BuildArgs)
  LL <- tryCatch(dlmLL(y=Y, mod=Model, debug=debug), 
                 error=function(e) { z <- Inf; attr(z, "error") <- e; z})
  if (!is.finite(LL)) {
    warning("dlmLL failed: returning BigVal with error attribute", immediate.=verbose)
    err <- attr(LL, "error")
    LL <- BigVal
    attr(LL, "error") <- err
  } else if (!is.null(prior.pars)) {
    if (nrow(prior.pars)!=length(theta) | ncol(prior.pars) !=2) 
      stop("prior.pars should be a 2-column matrix with a row per element of theta")
    LL <- LL - sum(dnorm(theta, mean=prior.pars[,1], sd=prior.pars[,2], log=TRUE))
  }
  if (verbose) {
   cat("(Penalised) log-likelihood from dlm.SafeLL: ", -round(LL,2),"\n")
   cat("Value of theta for this penalised log-likelihood:\n")
   print(theta)
   cat("\n")
  }
  LL
}
######################################################################
dlm.SafeMLE <- function(theta.init, Y, build, debug=FALSE, 
                        Use.dlm=FALSE, par.names=NULL, prior.pars=NULL, 
                        messages=TRUE, hessian=FALSE, ...) {
  #
  # "Safe" maximum likelihood estimation for a dynamic linear model. This
  # replaces the dlmMLE routine in the dlm library, which can be unstable. 
  # It uses nlm() rather than optim() for the optimisation. Arguments:
  #
  # theta.init    Initial value for the parameter vector
  # Y             The data
  # build         A function that will construct the dlm object 
  #               representing the model being fitted. Typically
  #               one of the .modeldef functions in this script. 
  # debug         As in dlmLL
  # Use.dlm       If TRUE, the function will use dlmMLE in place
  #               of nlm for the optimisation. This allows easy
  #               comparison of the results from the two 
  #               optimisation procedures. 
  # par.names     An optional character vector, the same length as
  #               theta.init, that will be used to label the 
  #               parameters in the output
  # prior.pars    If non-null, a two-column matrix with number of rows
  #               equal to the length of theta.init: the columns contain
  #               respectively the means and standard deviations of 
  #               independent Gaussian priors on the elements of theta.
  #               If NULL (the default), no prior information is used. 
  # messages      Either a logical or numeric scalar, controlling 
  #               the output of messages to screen during fitting
  #               when Use.dlm is FALSE. If FALSE or less than 1, no 
  #               messages will be produced. If TRUE or equal to 1, 
  #               basic messages will be produced. If numeric and 
  #               greater than 1, the value of messages-1 will be 
  #               passed as the print.level argument to nlm.
  #               
  # hessian       Controls whether to compute the hessian on exit
  #               from nlm(). The internal calculations within nlm()
  #               itself can be inaccurate, so this is done via a 
  #               call to the num.hess() routine (see above).
  # ...           Other arguments passed to build and nlm
  #
  # If Use.dlm is TRUE then the routine returns the result of a call
  # to dlmMLE(), with various control parameters set to stabilise the 
  # optimisation as far as possible. If Use.dlm is FALSE (the default), 
  # the routine returns the list result from a call to nlm (see below 
  # for what is optimised in this case), but with the names of key 
  # components modified to match the results of a call to optim 
  # (which is used by dlmMLE). The exception is the "code" component
  # of the result, which is retained and a "convergence" component is
  # added which is compatible with the convergence codes from optim. 
  # This ensures that the routine can be run interchangeably with 
  # dlmMLE. A component "HessFail" is also added, to record whether
  # the Hessian had to be tweaked to make it positive definite. 
  #
  # If Use.dlm is FALSE and prior.pars is non-NULL, the routine 
  # minimises the log posterior density (equivalently the 
  # minimimum penalised log-likelihood) for the parameters; 
  # in this case, the value of prior.pars is also included
  # in the list result for use by routines such as summary.dlmMLE().
  # otherwise it minimises the unpenalised log-likelihood. 
  # If use.dlm is TRUE then non-NULL values of prior.pars are 
  # ignored, with a warning. 
  #
  # Find the log-likelihood for the initial value:
  #
  if (Use.dlm) {
    if (!is.null(prior.pars)) warning("prior.pars is ignored when Use.dlm is TRUE")
    LL.init <- dlm.SafeLL(theta=theta.init, Y=Y, build=build, debug=debug, ...)
    if (as.numeric(messages)>=1) {
      cat("Initial value for theta:\n")
      print(signif(theta.init, 5))
      cat("Log-likelihood at initial value:", round(-LL.init,2),"\n")
    }
    z <- dlmMLE(y=Y, parm=theta.init, build=build, 
                control=list(fnscale=abs(LL.init), parscale=pmax(0.1, abs(theta.init))), ...)
  } else {
    #
    #   Here is a careful use of nlm. To protect against inappropriate 
    #   starting values, it's worth doing it in batches and, at the
    #   end of each batch, updating the parameter and objective
    #   scaling. The batch size of 20 ensures that the extra cost is
    #   about 5% which is probably tolerable. It's worth limiting 
    #   stepmax, because otherwise nlm() can fly off to far reaches 
    #   of the parameter space (noting that most parameters are 
    #   estimated on the log scale) - and increase it if necessary. 
    #   NB also: set hessian=FALSE here because we'll evaluate it 
    #   afterwards. And limit the size of "BigVal" (the value 
    #   returned if dlmLL fails) so that the algorithm will
    #   avoid it but it won't kill the gradient computations. 
    #
    NotDone <- TRUE
    LL.init <- NULL
    BatchSize <- 20
    StepMax <- 1
    StepTol <- 1e-6
    InitErr <- FALSE
    num.msg <- max(0, as.numeric(messages)-1)
    while(NotDone) {
      #
      #	  Need to handle ... arguments carefully here, because they could 
      #   relate to either dlm.SafeLL, build or nlm
      #
      ArgsGiven <- c(as.list(environment()), list(...))
      names(ArgsGiven)[names(ArgsGiven)=="theta.init"] <- "theta"
      ArgsNeeded <- c(formals(dlm.SafeLL), formals(build))
      ArgsNeeded <- ArgsNeeded[!duplicated(names(ArgsNeeded))]
      ArgsNeeded <- ArgsNeeded[!(names(ArgsNeeded)=="...")]
      LLArgs <- c(ArgsGiven[names(ArgsGiven) %in% names(ArgsNeeded)],
                  ArgsNeeded[!(names(ArgsNeeded) %in% names(ArgsGiven))])
      if (is.null(LL.init)) LL.init <- do.call(dlm.SafeLL, LLArgs)
      if (!is.null(par.names)) names(theta.init) <- par.names
      if (as.numeric(messages) > 0) {
        cat("Initial value for theta:\n")
        print(signif(theta.init, 5))
        cat("(Penalised) log-likelihood at initial value:", round(-LL.init,2),"\n")
      }
      #
      #   Bail out if initial value calculation failed
      #
      if (!is.null(attr(LL.init, "error"))) {
        warning(paste("dlm.SafeLL failed at initial value of theta -", 
                      "error message was\n ", 
                      attr(LL.init, "error")$message), immediate.=TRUE)
        z <- list(par=theta.init, value=Inf, code=-1, counts=0)
        if (!is.null(par.names)) names(z$par) <- par.names
        if (!(Use.dlm | is.null(prior.pars))) z$prior.pars <- prior.pars
        class(z) <- "dlmMLE"
        return(z)
      }
      z <- try(nlm(dlm.SafeLL, p=theta.init, Y=Y, build=build,
                   prior.pars=prior.pars, BigVal=LL.init+1, 
                   typsize=pmax(abs(theta.init), 0.1), 
                   fscale=abs(LL.init), stepmax=StepMax, steptol=StepTol,
                   hessian=FALSE, iterlim=BatchSize, gradtol=1e-6, 
                   print.level=num.msg, ...), silent=TRUE)
      if (isTRUE(class(z)=="try-error")) { # When nlm throws an Inf
        warning(paste("nlm failed: error message was", as.character(z), sep="\n"),
                immediate.=TRUE)
        BatchSize <- floor(BatchSize / 2)
        if (BatchSize > 0) {
          if (as.numeric(messages) > 0) {
            cat("Trying again, reducing number of nlm iterations for each batch ...\n")
          }
        } else {
          NotDone <- FALSE # Bail out
          warning("Unable to fix nlm failure: result is probably not global optimum",
                  immediate.=TRUE)
        }
      } else if (z$code==3) {
        theta.init <- z$estimate
        LL.init <- z$minimum
        StepTol <- StepTol / 10
        if (as.numeric(messages) > 0) {
          cat("Reducing StepTol and retrying optimisation from best value so far ...\n")
        }	  
      } else if (z$code %in% 4:5) { # Iteration limit exceeded: move theta.init & try again
        theta.init <- z$estimate
        LL.init <- z$minimum
        if (as.numeric(messages) > 0) {
          cat("Retrying optimisation from best value found so far ...\n")
          StepTol <- 1e-6
          StepMax <- 1
          BatchSize <- 20
        }
        if (z$code==5) {
          theta.init <- z$estimate
          LL.init <- z$minimum
          StepMax <- 1.2*StepMax # Increase StepMax if it was too small
          if (as.numeric(messages) > 0) {
            cat("Increasing value of StepMax in optimisation ...\n")
          }
        }
      } else NotDone <- FALSE
      if (hessian & !NotDone) { # Only do this next bit if we think we've finished
        #
        #   Calculate hessian as part of the "while" loop, because the extra Function
        #   evaluations may find better parameter values. While doing this, note that
        #   dlmLL sometimes throws infinite values which screw up the Hessian
        #   calculations. Operational workaround: replace these infinite values 
        #   with something slightly larger (i.e. a bit worse but not much) than 
        #   the optimum. This will tend to suggest limited curvature, hence 
        #   increasing standard errors etc.
        #
        #   Note also the QuitEarly argument: if the hessian calculations
        #   find an "appreciably better" log-likelihood then bail out and
        #   re-optimise. There's not much point in bailing out at an
        #   intermediate stage if the hessian finds a parameter value
        #   with a very small improvement, because in that case 
        #   re-optimising isn't likely to change things very much and
        #   we'll just be wasting time to get back where we started.
        #
        if (as.numeric(messages) > 0) cat("Now calculating hessian ... \n")
        HessCalc <- 
          num.hess(dlm.SafeLL, x0=z$estimate, Y=Y, build=build, 
                   prior.pars=prior.pars, BigVal=z$minimum+1e-3, 
                   verbose=(num.msg>1), relerr=1e-3, 
                   QuitEarly=list(reference=z$minimum-0.01, greater=FALSE), 
                   ...)
        if (any(HessCalc$details$Value < z$minimum)) {
          NotDone <- TRUE # OK, we haven't finished after all
          theta.init <- as.numeric(HessCalc$details[which.min(HessCalc$details$Value),-(1:4)])
          LL.init <- min(HessCalc$details$Value)
          if (as.numeric(messages) > 0) {
            cat("Found better parameters while calculating hessian: now re-optimising ...\n")
          } else {
            if (as.numeric(messages) > 0) cat("Hessian calculation complete.\n")
          }
        }
      }
    }
    names(z)[names(z)=="minimum"] <- "value"
    names(z)[names(z)=="estimate"] <- "par"
    names(z)[names(z)=="iterations"] <- "counts"
    z$convergence <- ifelse(z$code %in% 1:2, 0, ifelse(z$code==4, 1, 2))
  }
  if (z$convergence == 0) {
    if (as.numeric(messages) > 0) cat("Numerical optimiser reports successful convergence\n")
  } else {
    warning("numerical optimiser reports possible convergence problem", immediate.=TRUE)
  }
  if (!is.null(par.names)) {
    names(z$par) <- par.names
  }
  if (hessian) {
    z$hessian <- HessCalc$hessian
    #
    #   Check the hessian for positive definiteness; if it isn't PD
    #   then tweak it so that it *is*. Do the check using both 
    #   chol and eigen, because both may be used later (and they
    #   give slightly different results!)
    #
    z$HessFail <- isTRUE(class(try(chol(z$hessian), silent=TRUE)) == "try-error")
    if (!z$HessFail) z$HessFail <- !all(eigen(z$hessian)$values > 0)
    if (z$HessFail) {
      HessPD <- z$hessian
      EigenVals <- eigen(HessPD)$values
      if (sum(EigenVals)<0) {
        warning("Hessian isn't positive definite: check optimisation results", 
                immediate.=TRUE)
      } else {
        warning("Hessian isn't positive definite: shrinking eigenvalues towards their mean",
                immediate.=TRUE)
        while(min(EigenVals) < 0) {
          HessPD <- 0.95*HessPD + 
            ( (0.05*sum(EigenVals) / nrow(HessPD)) * diag(rep(1, nrow(HessPD))) )
          EigenVals <- eigen(HessPD)$values
        }
        z$hessian <- HessPD
      }
    }
    if (!is.null(par.names)) row.names(z$hessian) <- colnames(z$hessian) <- par.names
  }
  if (!(Use.dlm | is.null(prior.pars))) z$prior.pars <- prior.pars
  class(z) <- "dlmMLE"
  z
}
######################################################################
print.dlmMLE <- function(x, ...) {
  ######################################################################
  #
  # print method for dlm fits obtained using dlmMLE.
  #
  ######################################################################
  cat("Object of class dlmMLE, containing the following components:\n") 
  print(str(x), ...)
  cat("Use summary() method for information on the model fit.\n")
  invisible(x)
}
######################################################################
summary.dlmMLE <- function(object, ...) {
  ######################################################################
  #
  # summary method for a dlm fit obtained using dlmMLE.
  # Arguments:
  #
  # object		  The result of a call to dlmMLE
  # ...         For consistency with the summary() generic
  #
  # Value: a list containing the following components:
  #
  # model.table	Table of estimates and, if fit contains a
  #			"hessian" component (the result of calling
  #			dlmMLE with argument hessian=TRUE) 
  #			approximate standard errors
  # logLik		The log-likelihood for the fit
  # Corr		If fit contains a "hessian" component,
  #			approximate correlation matrix of the
  #			parameter estimates. 
  #
  ######################################################################
  have.hess <- "hessian" %in% names(object)
  if (have.hess) {
    covmat <- try(solve(object$hessian), silent=TRUE)
    if (isTRUE(class(covmat) == "try-error")) {
      Std.Errs <- cormat <- NULL
      cat("\nHessian of log-likelihood is singular: standard errors unavailable\n")
    } else {
      Std.Errs <- sqrt(diag(covmat))
      cormat <- cov2cor(covmat)
    }
  } else {
    cat("\nNo hessian in model fit: standard errors unavailable\n")
    Std.Errs <- cormat <- NULL
  }
  model.table <- rbind(object$par,Std.Errs)
  rownames(model.table) <- 
    if (!is.null(Std.Errs)) c("Estimate","SE") else "Estimate"
  logLik <- -object$value
  if (!have.hess) {
  }
  if (is.null(object$prior.pars)) {
    z <- list(model.table=model.table,Corr=cormat,logLik=logLik)
  } else {
    log.prior <- 
      sum(dnorm(object$par, object$prior.pars[,1], object$prior.pars[,2], log=TRUE))
    z <- list(model.table=model.table,Corr=cormat,
              logLik.Pen=logLik, logLik=logLik-log.prior)
  }
  class(z) <- "summary.dlmMLE"
  z
}
######################################################################
print.summary.dlmMLE <- function(x, digits=4, ...) {
  ######################################################################
  #
  # print method for objects of class summary.dlmMLE. Arguments:
  #
  # x	      	  The result of a call to summary.dlmMLE()
  # digits      Controls number of significant digits in output.
  # ...         For consistency with the summary() generic
  #
  ######################################################################
  cat("Parameter estimates:\n\n")
  print(signif(x$model.table, digits))
  if (!is.null(x$Corr)) {
    cat("\nCorrelation matrix of parameter estimates:\n\n")
    print(signif(x$Corr,digits-1))
  }
  if (is.null(x$logLik.Pen)) {
    cat(paste("\nLog-likelihood:",round(x$logLik,digits-2),"\n\n"))
  } else {
    cat(paste("\nPenalised log-likelihood (log posterior):",round(x$logLik.Pen, digits-2)))
    cat(paste("\nUnpenalised log-likelihood:",
                round(x$logLik, 2),"\n\n"))
  }
  invisible(x)
}
######################################################################
dlm.ObsPred <- function(ModelBundle) {
  #
  #  To calculate the predicted values, and associated error variances,
  #  for the *observed* quantities according to a state space time
  #  series model that has been used to run a Kalman Smoother.
  #  Arguments: 
  #
  #  ModelBundle  A list containing at least the named components
  #               "Model" (an object of class dlm) and "Smooth"
  #               (the result of Kalman Smoothing a dataset using 
  #               the dlm stored in Model, via the dlmSmooth 
  #               command). NB it's assumed that the model is
  #               set up in such a way that the "observed" 
  #               time series is the first element of each
  #               observation vector. 
  #
  #  The routine returns a list containing the predictions 
  #  and their standard errors, each returned as a matrix
  #  in which the rows represent time points and columns 
  #  represent series.
  # 
  Model <- ModelBundle$Model; Smooth <- ModelBundle$Smooth
  state.vars <- simplify2array(dlmSvd2var(Smooth$U.S,Smooth$D.S)[-1]) # Covariance matrices
  #
  #   Note that Smooth$s contains the estimated state vectors as *rows*
  #   rather than columns: so we postmultiply by t(FF) to get the 
  #   predicted values of the observed quantities instead of 
  #   premultiplying by FF.
  #
  preds <- (Smooth$s %*% t(Model$FF))[-1,,drop=FALSE]
  #
  #   In the next line, the apply() command usually produces a matrix
  #   with a row for every component of the observation vector and a 
  #   column for every time point; the entries are the variances of
  #   the state estimates. Then diag(Model$V) is a vector with an 
  #   element for every component; this gets recycled columnwise
  #   in the sum. 
  #
  #   The exception is when the observation time series is scalar-
  #   valued: in this case sapply() transposes the result. I can't
  #   figure out an elegant way to prevent this: the inelegant way
  #   is the if() in the next line but one. 
  #
  preds.var <- 
    apply(X=state.vars, MARGIN=3, 
          FUN=function(x,A,V) diag(covxfrm(x,A))+diag(V), 
          A=Model$FF, V=Model$V)
  if (nrow(Model$V)>1) {
    preds.var <- t(preds.var)
  } else preds.var <- matrix(preds.var, ncol=1)
  list(Pred=preds, SE=sqrt(preds.var))
}
######################################################################
dlm.ThetaSample <- function(Fit, N, Random=TRUE, Quantile=TRUE, 
                            Level=if (!Quantile) 1.96 else 0.99, 
                            df=NULL, Antithetic=c("Mean", "Scale")) {
  #
  #   Generates samples from the approximate posterior distribution
  #   for the parameter vector of a dynamic linear model. Arguments: 
  #
  #   Fit         A list containing the results of a numerical maximisation
  #               of a log-posterior. Typically the result of a call to 
  #               dlm.SafeMLE or equivalent. 
  #   N           Number of samples to draw. Ignored if Random is FALSE
  #               (see below).
  #   Random      If TRUE (the default) then samples are drawn 
  #               randomly from the approximate posterior. If
  #               FALSE, the routine returns the values m, m+/-A*u[1]*v[,1], 
  #               m+/-A*u[2]*v[,2], ... , m+/-A*u[d]*v[,d], where m is the
  #               posterior mode, d is the number of parameters in the model, 
  #               u is a vector containing the square roots of the singular 
  #               values of the posterior covariance matrix and v is a 
  #               matrix containing the corresponding singular vectors in its 
  #               columns. A is set using the values of Quantile and Level
  #               (see below). 
  #   Quantile    Used only if Random is FALSE. In this case, Quantile is
  #               a logical scalar determining how the "Level" argument 
  #               (see below) is used to calculate the value of A is
  #               calculated. 
  #   Level       A numeric scalar, used to calculate the value of A when
  #               Random is FALSE. In this case, if Quantile is also FALSE
  #               then Level is itself the value of A and the default 
  #               value of 1.96 delineates a 95% componentwise interval
  #               along each of the principal axes of the underlying
  #               covariance matrix. By contrast, if Quantile is TRUE then 
  #               A is chosen so that the sampled sampled points all lie 
  #               on the ellipsoid defining the approximate joint (Gaussian) 
  #               central credible region for the parameter vector: in this
  #               case, the default value of 0.99 specifies the 99% 
  #               credible region. 
  #
  #               This argument is ignored if Random is TRUE. 
  #   df          If non-NULL and if Random is TRUE, sampling will be 
  #               from a shifted t-distribution with df degrees of 
  #               freedom. In this case, the "mean" and "covariance matrix"
  #               of the posterior distribution are taken as the location
  #               and dispersion matrices of the t=distribution. 
  #   Antithetic  If TRUE, antithetic variable sampling will be used
  #               to reduce sampling variability if Random is TRUE. 
  #               The argument is ignored if Random is FALSE. The 
  #               argument can be either NULL (in which case 
  #               independent random sampling is used) or a 
  #               character vector containing a subset of 
  #               c("Mean", "Scale", "Corr"). If the argument
  #               contains "Mean" and not "Corr" then, for 
  #               every sampled value theta, the routine also 
  #               returns 2m-theta which is the reflection about 
  #               the posterior mode: this ensures that the 
  #               sample mean vector is equal to the posterior
  #               mean. If the argument contains "Corr" then
  #               the routine produces 2^d-2 structured samples
  #               for each independent value, by projecting the
  #               the initial value onto the principal axes of
  #               variation about the mode and then multiplying 
  #               each centred component by +/-1 (the single
  #               antithetic sample generated using "Mean"
  #               corresponds to the case when all the multipliers
  #               are 1): this distributes samples around the 
  #               ellipsoid with joint density equal to that of
  #               the initial value, thus nsuring that the 
  #               correlation shape of the underlying distribution
  #               is preserved. Finally, if the argument contains
  #               "Scale" then, for each sample theta (after
  #               potentially generating the "Mean" and "Corr"
  #               antithetics), the value m + sqrt(q2/q1)(theta-m)
  #               is also returned where q1 is the quadratic form
  #               giving the chi-squared quantile for the joint
  #               distribution of theeta-m, and q2 is the opposite
  #               tail quantile for the same distribution: this
  #               ensures that samples have the correct overall
  #               dispersion about the posterior mode.
  #
  #               If either of the antithetic options are used,
  #               the total number of samples produced is thus a 
  #               multiple of a power of 2: if this is greater than 
  #               N then a random subset of N samplee is selected. 
  #
  #   The routine returns a matrix in which each row contains a
  #   posterior sample; and has attributes "Mean", "Cov.SpD"
  #   and "df" giving the mean vector, covariance (dispersion) 
  #   matrix spectral decomposition and degrees of freedom for 
  #   the underlying distribution. 
  #
  #   Start by computing the spectral decomposition of the returned 
  #   Hessian matrix, check that it's non-negative definite and compute 
  #   the spectral decomposition of the covariance matrix, working 
  #   with a generalised inverse if necessary. 
  #
  if (is.null(Fit$hessian)) stop("Fit doesn't contain a hessian component")
  Cov.SpD <- eigen(Fit$hessian, symmetric=TRUE)
  if (any(Cov.SpD$values<0)) {
    stop("Fit$hessian isn't non-negative definite")
  }
  Zeroes <- Cov.SpD$values==0
  Cov.SpD$values[!Zeroes] <- 1/Cov.SpD$values[!Zeroes]
  d <- length(Cov.SpD$values) # length of parameter vector
  design <- matrix(rep(1, d), nrow=1) # Basic "design matrix" for a single sample
  if (!Random) { # Design points for deterministic samples
    design <- as.matrix(expand.grid(rep(list(c(-1,1)), d)))
  }
  if ("Corr" %in% Antithetic) {
    design <- as.matrix(expand.grid(rep(list(c(-1,1)), d)))
    if ("Mean" %in% Antithetic) {
      warning("'Mean' antithetic option ignored when 'Corr' is specified")
    }
  } else if ("Mean" %in% Antithetic) {
    design <- rbind(design, -design)
  }
  design.scale <- rep(1, nrow(design))  # To allow for antithetic scale sampling
  if ("Scale" %in% Antithetic) {
    design <- design[rep(1:nrow(design), each=2),]
    design.scale <- c(1,-1) * rep(design.scale, each=2) # Uses recycling of c(1,-1)
  }
  if (Random) {
    NReps <- nrow(design) # Number of correlated samples for each independent one
    NInd <- ceiling(N / NReps)
    z <- matrix(rnorm(NInd*d), ncol=d) # Matrix filled with standard normals
    if (!is.null(df)) {
      denom <- sqrt(rchisq(NInd, df=df) / df)
      z <- z / denom # Uses columnwise recycling, so each row gets the same divisor
    }
    if (!is.null(Antithetic)) { # Create "near-replicates" if antithetic sampling is used
      z <- z[rep(1:nrow(z), each=nrow(design)),] # Replicate generated samples
      design <- design[rep(1:nrow(design),NInd),]
      z <- design*z # Elementwise multiplication does the reflections
      if ("Scale" %in% Antithetic) {
        q1 <- rowSums(z*z)[design.scale==-1] # chi-squared quantiles for each sample
        q2 <- qchisq(pchisq(q1, d, lower.tail=TRUE), d, lower.tail=FALSE)
        z[design.scale==-1,] <- # This line recycles elements of design.scale
          sqrt(q2/q1) * z[design.scale==-1,]
      }
      z <- z[sample(1:nrow(z), N),,drop=FALSE] # Thin to leave N samples
    }
  } else { # Non-random sample; need value of A, stored in the object NSDs
    if (Quantile) {
      NSDs <- sqrt(qchisq(Level, d) / d)
    } else {
      NSDs <- Level
    }
    z <- rbind(rep(0,d), design*NSDs)
  }
  #
  # Postmultiply by square root of covariance matrix (post- rather than
  # premultiply because results are stored as rows i.e. the transposes
  # of the usual vectors). The "z*rep(...)" part uses columnwise 
  # recycling to avoid a matrix multiplication. 
  #
  e <- (z*rep(sqrt(Cov.SpD$values), each=nrow(z))) %*% t(Cov.SpD$vectors)
  #
  # Add mean vector to each sample
  #
  mu.mat <- matrix(rep(Fit$par, nrow(e)), nrow=nrow(e), byrow=TRUE)
  z <- mu.mat + e
  attr(z, "Mean") <- Fit$par
  attr(z, "Cov.SpD") <- Cov.SpD
  attr(z, "df") <- df
  colnames(z) <- names(Fit$par)
  z
}######################################################################
dmvnorm.SpD <- function(x, mu, Cov.SpD, logged=FALSE) {
  #
  #   To calculate the density of a multivariate normal distribution
  #   parameterised in terms of its mean vector and the spectral
  #   decomposition (eigendecomposition) of its covariance matrix. 
  #   Arguments: 
  #
  #   x       An n*d matrix, where each row represents a point in R^d
  #           for which the density is required.
  #   mu      Vector of length d: the mean vector of the distribution
  #   Cov.SpD A list, with components "values" and "vectors", 
  #           defining the spectral decomposition of the covariance
  #           matrix of the distribution. Typically the result of
  #           a call to eigen(). The eigenvalues must all be 
  #           strictly positive.
  #   logged  Return the density or its logarithm?
  #
  #   Value: a vector of length n, containing the (log) densities.
  #
  #   Output from this routine has been cross-checked against
  #   that from dmvnorm in the mvtnorm library. 
  #
  if (any(Cov.SpD$values <= 0)) stop("Covariance matrix is not positive definite")
  #
  #   Evaluate the log density first, starting with the quadratic form which, 
  #   for an individual d*1 vector x, is -(x-mu)'VL^{-1/2}L^(-1/2)V'(x-mu) / 2 
  #   where V & L are the matrices of eigenvectors and values respectively. 
  #   In this function, x is a matrix with *rows* containing the evaluation
  #   points. 
  #
  n <- nrow(x); d <- ncol(x)
  mu.mat <- matrix(rep(mu, n), nrow=n, byrow=TRUE)
  ZScores <- (x - mu.mat) %*% Cov.SpD$vectors / # Division here exploits
    rep(sqrt(Cov.SpD$values), each=n)           # columnwise recycling
  f <- -( (d*log(2*pi)) + sum(log(Cov.SpD$values)) + rowSums(ZScores*ZScores) ) / 2
  if (!logged) f <- exp(f)
  f
}
######################################################################
dmvt.SpD <- function(x, mu, Cov.SpD, df, logged=FALSE) {
  #
  #   To calculate the density of a shifted multivariate t distribution
  #   parameterised in terms of its median vector, spectral decomposition
  #   (eigendecomposition) of its dispersion matrix and degrees of freedom. 
  #   Arguments: 
  #
  #   x       An n*d matrix, where each row represents a point in R^d
  #           for which the density is required.
  #   mu      Vector of length d: the mean vector of the distribution
  #   Cov.SpD A list, with components "values" and "vectors", 
  #           defining the spectral decomposition of the dispersion
  #           matrix of the distribution. Typically the result of
  #           a call to eigen(). The eigenvalues must all be 
  #           strictly positive.
  #   df      The degrees of freedom of the distribution. 
  #   logged  Return the density or its logarithm?
  #
  #   Value: a vector of length n, containing the (log) densities. 
  #
  #   Output from this routine has been cross-checked against
  #   that from dmvt in the mvtnorm library. 
  #
  if (any(Cov.SpD$values <= 0)) stop("Dispersion matrix is not positive definite")
  #
  #   Evaluate the log density first, starting with the quadratic form which, 
  #   for an individual d*1 vector x, is -(x-mu)'VL^{-1/2}L^(-1/2)V'(x-mu) / 2 
  #   where V & L are the matrices of eigenvectors and values respectively. 
  #   In this function, x is a matrix with *rows* containing the evaluation
  #   points. 
  #
  n <- nrow(x); d <- ncol(x)
  mu.mat <- matrix(rep(mu, n), nrow=n, byrow=TRUE)
  ZScores <- (x - mu.mat) %*% Cov.SpD$vectors / # Division here exploits
    rep(sqrt(Cov.SpD$values), each=n)           # columnwise recycling
  f <- -(df + d) * log(1 + ( rowSums(ZScores*ZScores) / df )) / 2
  f <- f - ( (d*log(df*pi)) + sum(log(Cov.SpD$values)) ) / 2
  f <- f + lgamma((df+d)/2) - lgamma(df/2)
  if (!logged) f <- exp(f)
  f
}
######################################################################
dlm.ImportanceWts <- function(samples, build, Y, prior.pars=NULL, 
                              debug=FALSE, ...) {
  #
  #   Calculates importance weights for sampling from the
  #   posterior distribution of the parameter vector in a
  #   dynamic linear model. Arguments: 
  #
  #   samples A matrix, containing samples from the Gaussian
  #           approximation to the posterior and with attributes
  #           "Mean", "Cov.SpD" and "df" as produced by 
  #           dlm.ThetaSample(). 
  #   build   A function that will construct the dlm object 
  #           representing the model being fitted. Typically
  #           one of the .modeldef functions in this script. 
  #   Y       The data: a matrix containing a column for each
  #           time series. This will be used to calculate the
  #           exact log-posterior for each element of theta. 
  #   prior.pars  As in dlm.SafeLL.
  #   debug   As in dlmLL. 
  #   ...     Other arguments to build. 
  #
  #   Value: a data frame containing three columns and a row for each
  #   row of theta. The columns are log.g containing the values
  #   of log(g(theta)/g(theta.hat)) where theta.hat is the 
  #   posterior mode; log(h(theta)/h(theta.hat)); and the 
  #   corresponding normalised importance weights. 
  #
  if (is.null(attr(samples, "df"))) {
    log.g <- dmvnorm.SpD(samples, attr(samples, "Mean"), 
                         attr(samples, "Cov.SpD"), logged=TRUE)
    log.ghat <- dmvnorm.SpD(matrix(attr(samples, "Mean"),nrow=1), 
                            attr(samples, "Mean"), 
                            attr(samples, "Cov.SpD"), logged=TRUE)
  } else {
    log.g <- dmvt.SpD(samples, attr(samples, "Mean"), 
                      attr(samples, "Cov.SpD"), attr(samples, "df"),
                      logged=TRUE)
    log.ghat <- dmvt.SpD(matrix(attr(samples, "Mean"),nrow=1), 
                         attr(samples, "Mean"), attr(samples, "Cov.SpD"),
                         attr(samples, "df"), logged=TRUE)
  }
  #
  #   Use shifted data to calculate log-likelihoods, for improved
  #   accuracy; and take care over arguments, in case user has
  #   passed RunsPerTS in ... when Y already has it as an 
  #   attribute.
  #
  Ybar <- mean(Y[,1], na.rm=TRUE)
  DotArgs <- list(...)
  if (!is.null(attr(Y, "RunsPerTS"))) {
    if (!is.null(DotArgs$RunsPerTS)) {
      warning("Ignoring argument RunsPerTS and using data attribute instead")
      DotArgs <- DotArgs[-which(names(DotArgs)=="RunsPerTS")]
    }
  }
  DotArgs$RunsPerTS <- attr(Y, "RunsPerTS")
  log.h <- 
    -do.call(apply, c(list(X=samples, MARGIN=1, FUN=dlm.SafeLL,
                           Y=Y-Ybar, build=build, prior.pars=prior.pars,
                           Shift=Ybar, debug=debug), DotArgs))
  log.hhat <- 
    -do.call(dlm.SafeLL, c(list(theta=attr(samples, "Mean"), Y=Y-Ybar,
                                build=build, prior.pars=prior.pars, 
                                Shift=Ybar, debug=debug), DotArgs))
  #
  #   The unnormalised weights are now exp(log.h-log.g); these can 
  #   be large due to the arbitrary constants omitted from the
  #   log-likelihood / log posterior. Normalise therefore, and
  #   stabilise the calculations by subtracting log.ghat and log.hhat
  #   respectively (all this does is to scale the unnormalised
  #   weights by a constant, which doesn't matter given that we're
  #   going to normalise them anyway). Subtracting log.hhat from
  #   the log posteriors also enables us to identify the important
  #   region of parameter space, later. 
  #
  log.g <- log.g-log.ghat
  log.h <- log.h-log.hhat
  if (any(log.h>0)) {
    warning(paste("Log posterior for some samples exceeds that for MAP estimate:",
    "\n  optimisation must have failed"))
  }
  w <- exp(log.h - log.g)
  w <- w / sum(w) # Normalise for convenience
  z <- data.frame(log.g=log.g, log.h=log.h, w=w)
  z
}
######################################################################
CumulativeWeightPlot <- 
  function(w, main="Contributions of importance sampling weights", ...) {
  #
  #   To plot the cumulative contributions of a vector of
  #   weights, sorted in descending order. Arguments:
  #
  #   w     The vector of weights
  #   main  Title for the plot
  #   ...   Other arguments to plot.stepfun(). 
  #
  NWts <- length(w)
#  z <- plot(stepfun(100*(1:NWts)/NWts, 
#               cumsum(100*c(0,sort(w, decreasing=TRUE)))),
#       xval=100*(0:NWts)/NWts, xlim=c(0,100), lwd=2, 
#       col="darkblue", do.points=FALSE, verticals=TRUE,
#       xlab="Cumulative % of samples", 
#       ylab="Cumulative % of total weight", main=main, ...)
  z <- stepfun(100*(1:NWts)/NWts, cumsum(100*c(0,sort(w, decreasing=TRUE))))
  plot(c(0,100), c(0,100), type="n", xlim=c(0,100), ylim=c(0,100),
       xlab="Cumulative % of samples", 
       ylab="Cumulative % of total weight", main=main, ...)
  abline(h=seq(0,100,20), col=grey(0.8), lty=2)
  abline(v=seq(0,100,20), col=grey(0.8), lty=2)
  lines(z, xval=100*(0:NWts)/NWts, col="darkblue", 
  do.points=FALSE, verticals=TRUE, lwd=2)
  invisible(z)
}
######################################################################
PlotImportanceWts <- function(samples, weights, colours=NULL, 
                                   oma=c(7,3,5,3), ...) {
  #
  #   To produce a pairs plot of sampled parameter values from 
  #   a dlm, with points shaded according to the value of 
  #   their importance sampling weights. Arguments: 
  #
  #   samples   A matrix containing sampled parameter values, 
  #             with attributes as set by (e.g.) dlm.ThetaSample()
  #   weights   A 3-column matrix in which the third column 
  #             contains importance weights for each row of 
  #             samples, as produced by (e.g.) dlm.ImportanceWts.
  #   colours   A vector of colours to use when shading the 
  #             plots. If NULL (the default), a blue-red scale 
  #             is used, with dark blue denoting points with 
  #             very low importance weights and dark red denoting
  #             those with very high weights. 
  #   oma       This is a fudge, used to define the "outer margins"
  #             of the plot area in order to put a legend at the
  #             bottom (the underlying pairs() command doesn't 
  #             allow for this). The default value may need to be
  #             tweaked, depending on the dimensions of the output
  #             device. 
  #
  if (nrow(weights) != nrow(samples)) {
    stop("weights and samples have differing numbers of rows")
  }
  if (is.null(colours)) {
    transp <- 0.2+(abs( 4*((1:11)-6))/25)
    colours <- hcl.colors(11, palette="Blue-Red")
  }
  breakpts <- quantile(weights$w, probs=ppoints(length(colours)))
  colour.bands <- findInterval(weights$w, vec=breakpts, all.inside=TRUE)
  point.cols <- colours[colour.bands]
  opar <- par(no.readonly=TRUE)
  par(xpd=NA)
  pairs(samples, col=point.cols, pch=1, cex=0.8, oma=oma, ...)
  legend.quantiles <- c(0,0.25,0.5,0.75,1)
  leg.cols <- findInterval(quantile(weights$w, probs=legend.quantiles), 
                           vec=breakpts, rightmost.closed=TRUE, all.inside=TRUE)
  leg.cols <- colours[leg.cols]
  legend(x=0.5, y=0, xjust=0.5, yjust=1, fill=leg.cols, 
         legend=legend.quantiles, horiz=TRUE, cex=0.8, bty="n",
         title="Quantiles of importance weight distribution")
  par(opar)
}
######################################################################
SampleStates <- function(Thetas, build, Y, NonNeg=NULL, 
                         debug=FALSE, messages=FALSE, ...){
  #
  # To sample from the posterior distribution of the state 
  # vector for a dlm. This is essentially a wrapper for
  # dlmBSample(). Arguments:
  #
  #   Thetas  A matrix containing sampled parameter values, 
  #           as produced by (e.g.) dlm.ThetaSample()
  #   build   A function that will construct the dlm object 
  #           representing the model being fitted. Typically
  #           one of the .modeldef functions in this script. 
  #   Y       The data: a matrix containing a column for each
  #           time series. This will be used to calculate the
  #           exact log-posterior for each element of theta. 
  #   NonNeg  An optional vector containing indices of state
  #           elements that should not be negative (e.g. 
  #           when analysing untransformed precipitation series).
  #          
  #   debug   As in dlmLL.
  #   messagesControls whether to print progress to screen. 
  #   ...     Additional arguments to build().
  #
  # The function returns an array of dimension c(N, T, p)
  # where N is the number of samples, T the number of time
  # points in Y, and p the dimension of the state vector. 
  # Samples for which the Kalman Filter failed are left as
  # NA, with a warning, as are samples containing negative
  # values in the columns selected by NonNeg if this is 
  # non-NULL. NOTE that the routine removes the initial
  # "time zero" values inserted by the dlm library. 
  #
  if (!is.matrix(Y)) stop("Y must be a matrix")
  #
  #   The smart way to do this would be via apply(), but
  #   I can't figure out how to pass the arguments through
  #   to build() when doing that; another potential 
  #   disadvantage is that it doesn't allow printing of
  #   progress to screen. In any case, looping here
  #   probably isn't too inefficient given that each 
  #   iteration takes a bit of time. So: start by defining
  #   an array to store the samples - for which we need to
  #   find the dimension of the state vector.
  #
  DotArgs <- list(...)
  if (!is.null(attr(Y, "RunsPerTS"))) {
    if (!is.null(DotArgs$RunsPerTS)) {
      warning("Ignoring argument RunsPerTS and using data attribute instead")
      DotArgs <- DotArgs[-which(names(DotArgs)=="RunsPerTS")]
      DotArgs$RunsPerTS <- attr(Y, "RunsPerTS")
    }
  }
  Model <- do.call(build, c(list(theta=Thetas[1,]), DotArgs))
  if (ncol(Y) != nrow(Model$FF)) {
    stop("Dimension mismatch between Y and model definition in build")
  }
  NSamp <- nrow(Thetas)
  Nt <- nrow(Y)
  p <- ncol(Model$FF)
  res <- array(dim=c(NSamp, Nt, p),
               dimnames=list(Sample=1:NSamp, Time=1:Nt, StateEl=1:p))
  FailedIDs <- numeric(0)
  NegWarned <- FALSE
  #
  #   Now ready to go. For numerical stability, do the calculations
  #   on shifted data (subtract the mean of the first column of Y),
  #   then add back the necessary adjustments to all state vectors.
  #
  Ybar <- mean(Y[,1], na.rm=TRUE)
  for (i in 1:NSamp) {
    if (as.numeric(messages)>0) {
      cat(paste("Processing sample",i,"of",NSamp,"...\r"))
    }
    Model <- 
      do.call(build, c(list(theta=Thetas[i,], Shift=Ybar), DotArgs))
    Mod.Filtered <- try(dlmFilter(y=Y-Ybar, mod=Model, debug=debug),
                        silent=TRUE)
    if (!inherits(Mod.Filtered, "try-error")) {
      NegOK <- FALSE
      Attempt <- 0
      while (!NegOK) { # If NonNeg is TRUE, repeat until no -ve values are obtained; if it's FALSE then just go through once. 
        Sample <- try(dlmBSample(Mod.Filtered), silent=TRUE)
        if (!inherits(Sample, "try-error")) {
          Sample <-  # Add the shift back
            t(t(Sample) + attr(Model, "Shift"))
          res[i,,] <- NA # May need resetting in "while" loop
          res[i,,] <- Sample[-1,] # First row is "time zero"
          if (!isTRUE(all(abs(Sample[-1,])<1e12))) FailedIDs <- c(FailedIDs, i) # Check for exploding values & NAs
          NegIDs <- apply(Sample, MARGIN=2, FUN=function(x) any(x <= 0))
          # 
          # Next line: NonNeg is either NULL (no non-negative requirement)
          # or a vector picking out elements that should be non-negative.
          #
          if (is.null(NonNeg) | !isTRUE(any(NegIDs[NonNeg]))) { # 
            NegOK <- TRUE
          } else {
            if (!NegWarned) warning("Negative states have been generated: retrying ...")
            NegWarned <- TRUE
            if (Attempt==10) {
              NegOK=TRUE # This isn't going to work: give up
              res[i,,] <- NA
              FailedIDs <- c(FailedIDs, i)
            }
            Attempt <- Attempt + 1
          }
        } else {
          FailedIDs <- c(FailedIDs, i)
          NegOK <- TRUE # No point trying again if dlmBSample failed
        }
      }
    } else FailedIDs <- c(FailedIDs, i)
  }
  if (as.numeric(messages)>0) cat("\n")
  attr(res, "FailedIDs") <- NULL
  NFail <- length(FailedIDs)
  if (NFail > 0) {
    warning(paste("Sampling failed for",NFail,
                  "cases: their IDs are in the 'FailedIDs' attribute\n",
                  " of the result, and their values have been set to NA"))
    attr(res,"FailedIDs") <- FailedIDs
  }
  res
}
######################################################################
ExamineFailures <- function(Thetas, FailedIDs, alpha=0.2) {
  #
  # To examine the parameter values for which SampleStates()
  # failed due to problems with the Kalman Filter (typically
  # associated with near-singular matrix calculations).
  # Arguments:
  #
  # Thetas      A matrix containing sampled parameter values,
  #             as produced by (e.g.) dlm.ThetaSample()
  # FailedIDs   A vector containing the row numbers of
  #             ThetaSamples for which SampleStates() failed. 
  # alpha       Transparency to use for plotting non-failed
  #             samples. Values close to 1 are opaque, those
  #             close to 0 are more transparent and create
  #             the impression of a hazy background data cloud
  # 
  # The function produces a pairs plot, with the majority of 
  # the data points plotted as a background "data cloud" in
  # grey; and with the points for the failed samples plotted
  # in the foreground in red. This is primarily intended for
  # diagnostic purposes, hence few user-customisable options. 
  #
  N <- nrow(Thetas); Nf <- length(FailedIDs) 
  plot.chars <- c(rep(1,N-Nf), rep(16, Nf))
  plot.cols <- c(rep(grey(0.8, alpha=0.2), N-Nf), rep("red", Nf))
  if (Nf > 0) {
    pairs(rbind(Thetas[FailedIDs,], Thetas[-FailedIDs,]),
          col=plot.cols, pch=plot.chars)
  } else {
    pairs(Thetas, col=plot.cols[1], pch=plot.chars[1])
    warning("You can't see failed points on this plot, because there aren't any!")
  }
}
######################################################################
SampleObs <- function(Thetas, States, build, Y, WhichEls=1:ncol(Y), 
                      NonNeg=FALSE, ReplaceAll=FALSE, ...){
  #
  # To sample from the posterior predictive distribution of 
  # (elements of) the observation vector for a dlm. Arguments:
  #
  #   Thetas      A matrix containing sampled parameter values, 
  #               as produced by (e.g.) dlm.ThetaSample()
  #   States      An array containing sampled time series of 
  #               the state vector for each row of ThetaSamples,
  #               produced by SampleStates()
  #   build       Function to construct the dlm object being used.
  #               Typically one of the .modeldef functions in 
  #               this script. 
  #   Y           The data: a matrix containing a column for each
  #               time series. 
  #   WhichEls    Numeric vector indicating which column(s) of
  #               Y to simulate.
  #   NonNeg      Indicates whether the values to be returned 
  #               should be non-negative. If not (the default), 
  #               sampling is done from a normal distribution.
  #               Otherwise the samples are then transformed to
  #               the corresponding quantiles of a gamma 
  #               distribution with the same mean and variance. 
  #   ReplaceAll  If TRUE, the returned object will contain 
  #               random samples for all time points. If FALSE
  #               (the default), only the values missing from 
  #               Y are simulated (and the non-missing values
  #               are included in the output)
  #   ...         Additional arguments to build().
  #
  # The function returns an array of dimension c(T, p, N)
  # where N is the number of samples, T the number of time
  # points in Y, and p the number of series being simulated
  # (i.e. the length of WhichEls). The ordering of dimensions
  # ensures that if a single series is being simulated then
  # the result will be a matrix with the series in columns. 
  #
  if (!is.matrix(Y)) stop("Y must be a matrix")
  if(nrow(Thetas) != dim(States)[1]) {
    stop("Numbers of samples in Thetas and States don't match")
  }
  #
  #   As with SampleStates(), it's easiest to do this via
  #   a loop - which, again, is probably not too inefficient
  #   here. Start by defining an array to store the samples.
  #
  NSamp <- nrow(Thetas)
  Nt <- nrow(Y)
  p <- length(WhichEls)
  res <- array(dim=c(Nt, p, NSamp),
               dimnames=list(Time=1:Nt, Series=1:p, Sample=1:NSamp))
  NonMiss <- !is.na(Y[,WhichEls])
  #
  #   And, as with SampleStates, need to check for unnecessary 
  #   user-supplied arguments.
  #
  DotArgs <- list(...)
  if (!is.null(attr(Y, "RunsPerTS"))) {
    if (!is.null(DotArgs$RunsPerTS)) {
      warning("Ignoring argument RunsPerTS and using data attribute instead")
      DotArgs <- DotArgs[-which(names(DotArgs)=="RunsPerTS")]
      DotArgs$RunsPerTS <- attr(Y, "RunsPerTS")
    }
  }
  for (i in 1:NSamp) {
    Model <- do.call(build, c(list(theta=Thetas[i,]), DotArgs))
    FF <- Model$FF[WhichEls,]
    V.SpD <- eigen(Model$V[WhichEls, WhichEls])
    #
    #   State vectors need to be "vertical", hence transpose
    #   of States matrix
    #
    Means <- FF %*% t(States[i,,])
    z <- matrix(rnorm(Nt*p), ncol=Nt) # Matrix filled with standard normals
    z <- V.SpD$vectors %*% (z*sqrt(V.SpD$values)) # Uses columnwise recycling
    tmp <- Means + z
    if (NonNeg) { # For non-negative quantities, map to quantiles of gamma
      SDs <- sqrt(Model$V[WhichEls, WhichEls])
      z <- pnorm(tmp, mean=Means, sd=SDs)
      Shapes <- (Means/SDs)^2
      Rates <- Means/(SDs^2)
      tmp <- qgamma(z, shape=Shapes, rate=Rates)
    }
    if (!ReplaceAll) tmp[NonMiss] <- (Y[,WhichEls])[NonMiss]
    res[,,i] <- t(tmp)
  }
  res[,,,drop=TRUE] # Drop redundant dimensions
}
######################################################################
PostPredSample <- function(ModelBundle, Y=NULL, Build, N,  
                           Random=TRUE, Quantile=TRUE, 
                           Level=ifelse(Quantile, 1.96, 0.99),
                           df=NULL, Antithetic=c("Mean", "Scale"),
                           Importance=FALSE, CheckMax=FALSE, 
						   ReplaceOnFail=Random, PlotFails=FALSE, 
						   WhichEls=1, NonNeg=FALSE, ReplaceAll=FALSE, 
						   debug=FALSE, messages=FALSE, ...) {
  #
  #   To draw samples from the posterior predictive 
  #   distribution of observable time series (observations
  #   or ensemble members) using a fitted dlm. Arguments: 
  #
  #   ModelBundle   A list containing at least the named components
  #                 "Model" (an object of class dlm) and "Theta"
  #                 (a list containing a maximimum likelihood or 
  #                 MAP fit of a state space model). NB it's assumed
  #                 that the model is set up in such a way that the 
  #                 "observed"  time series is the first element of each
  #                 observation vector. The *Smooth() routines
  #                 in this script produce objects with the 
  #                 the required structure. 
  #   Y             Optional matrix of data, in which the first
  #                 column contains an observed time series
  #                 and the remainder contain ensemble members. If
  #                 NULL, Y will be taken from the Data$Y component
  #                 of ModelBundle if this exists.
  #   Build         A function that will construct the dlm object 
  #                 representing the model being fitted. Typically
  #                 one of the .modeldef functions in this script.
  #   N             Number of posterior samples to draw. Ignored 
  #                 if Random is FALSE (see below).
  #   Random        If TRUE (the default) then samples are drawn 
  #                 randomly from the (approximate) posterior. If
  #                 FALSE, a systematic sample is taken. See the
  #                 header to dlm.ThetaSample() for more details.  
  #   Quantile      Used only if Random is FALSE. See dlm.ThetaSample(). 
  #   Level         Used if Random is FALSE. See dlm.ThetaSample().
  #   df            If non-NULL and if Random is TRUE, sampling will be 
  #                 from a shifted t-distribution with df degrees of 
  #                 freedom. Otherwise a multivariate Gaussian is
  #                 used. 
  #   Antithetic    Controls the use of antithetic variables to
  #                 reduce sampling variation if Random is TRUE. 
  #                 The argument is ignored if Random is FALSE. 
  #                 See dlm.ThetaSample() for full details. 
  #   Importance    Logical scalar, indicating whether or not
  #                 to use importance sampling. 
  #   CheckMax		If TRUE and if Importance is TRUE, the routine
  #                 will check that none of the sampled parameter 
  #					values has a higher (penalised) log-likelihood
  #                 than ModelBundle$Theta
  #   ReplaceOnFail Logical scalar, used to control the behaviour
  #                 of the routine if the posterior sampling of 
  #                 states fails for some parameter sets. If TRUE,
  #                 these sets are replaced with alternative samples;
  #                 otherwise they are retained and the corresponding
  #                 samples of states and observable are set to NA. 
  #   PlotFails     If TRUE and if the posterior sampling of states
  #                 failed for any parameter sets, produce a plot
  #                 on the current graphics device to examine
  #                 the failed sets.
  #   WhichEls      Numeric vector choosing the series (i.e. 
  #                 columns of Data) for which sampling is
  #                 required. The default value of 1 selects
  #                 just the observations. 
  #   NonNeg        If TRUE, the sampled "observations" will be
  #                 drawn from gamma distributions with mean
  #                 and variance matching the correct values. 
  #                 If FALSE (the default), they will be drawn
  #                 from normal distributions. 
  #   ReplaceAll    If TRUE, the returned object will contain 
  #                 random samples for all time points. If FALSE
  #                 (the default), only the values missing from 
  #                 Data are simulated (and the non-missing values
  #                 are included in the output)
  #   debug         As in dlmLL.
  #   messages      Controls whether to print progress to screen. 
  #   ...     Additional arguments to build().
  #
  #   The routine is essentially a wrapper to dlm.ThetaSample(),
  #   SampleStates() and SampleObs(). It returns a list
  #   containing elements Thetas, States and Obs which are
  #   the results of the calls to these three functions; and
  #   Model (the dlm model object used to generate them).
  #   In addition, it contains elements Weights (the result 
  #   of calling dlm.ImportanceWts() if Importance is TRUE, 
  #   and NULL otherwise) and Replacements (a logical 
  #   scalar indicating whether some of the initial samples
  #   were replaced). 
  #
  #   Start by extracting the data. 
  #
  if (!is.null(Y)) {
    if (!is.matrix(Y)) stop("Y must be a matrix")
    Data <- Y
  } else {
  #
  #   Where data come from ModelBundle, need to ensure no
  #   conflicts with user-supplied arguments
  #
    Data <- ModelBundle$Data$Y
    if (is.null(Data)) {
      stop("ModelBundle contains no Data$Y component")
    }
    Data <- as.matrix(Data)
    attr(Data, "RunsPerTS") <- attr(ModelBundle$Data$Y, "RunsPerTS")
  }
  DotArgs <- list(...)
  if (!is.null(DotArgs$RunsPerTS) & !is.null(attr(Data, "RunsPerTS"))) {
    warning("Ignoring argument RunsPerTS and using data attribute instead")
    DotArgs <- DotArgs[-which(names(DotArgs)=="RunsPerTS")]
  }
  if ("Groups" %in% names(ModelBundle$Data)) {
    if (!is.null(DotArgs$Groups)) {
      warning("Ignoring argument Groups and using data attribute instead")
    }
    DotArgs$Groups <- ModelBundle$Data$Groups
  }
  if ("Groups" %in% names(DotArgs) & # Incredibly annoying, this one!
    !("Groups" %in% names(formals(Build)))) {
    DotArgs <- DotArgs[-which(names(DotArgs)=="Groups")]
  }
  if (ReplaceOnFail & !Random) {
    stop("ReplaceOnFail can't be TRUE unless Random is also TRUE")
  }
  if (CheckMax & !Importance) {
    warning("CheckMax is ignored if Importance is FALSE")
  }
  z <- list(Thetas=NULL, Samples=NULL, Obs=NULL, Weights=NULL,
            Replacements=FALSE)
  if (as.numeric(messages)>0) cat("Sampling parameter sets ...\n")
  z$Thetas <- 
    dlm.ThetaSample(ModelBundle$Theta, N=N, Random=Random,
                    Quantile=Quantile, Level=Level, df=df,
                    Antithetic=Antithetic)
  if (Importance) {
    if (as.numeric(messages)>0) cat("Calculating importance weights ...\n")
    IWArgs <- c(list(samples=z$Thetas, build=Build, Y=Data, 
                        prior.pars=ModelBundle$Theta$prior.pars,
                        debug=debug), DotArgs)
    z$Weights <- do.call(dlm.ImportanceWts, IWArgs)
	  if (CheckMax & any(z$Weights$log.h>0)) {
	    warning(paste("See warning from dlm.ImportanceWts: quitting",
	                  "PostPredSample without\n  generating 'States' or 'Obs'",
	                  "components of result. 'Thetas' component has",
	                  "\n  'BestID' attribute giving number of row",
	                  "with highest log-posterior."))
	    attr(z$Thetas, "BestID") <- which.max(z$Weights$log.h)			
	    return(z)
	  }
  }
  if (as.numeric(messages)>0) cat("Sampling state vectors ...\n")
  SSArgs <- c(list(Thetas=z$Theta, build=Build, Y=Data, 
                     debug=debug, messages=messages, 
                     NonNeg=if (NonNeg) WhichEls else NULL), DotArgs)
  z$States <- do.call(SampleStates, SSArgs)
  FailedIDs <- attr(z$States, "FailedIDs")
  if (PlotFails) ExamineFailures(z$Thetas, FailedIDs)
  if (!is.null(FailedIDs)) {
    if (ReplaceOnFail) {
      while (!is.null(FailedIDs)) {
        NFails <- length(FailedIDs)
        NewThetas <-     
          dlm.ThetaSample(ModelBundle$Theta, N=NFails, Random=Random,
                          Quantile=Quantile, Level=Level, df=df,
                          Antithetic=Antithetic)
        z$Thetas[FailedIDs,] <- NewThetas
        if (Importance) {
          IWArgs$samples <- NewThetas
          z$Weights[FailedIDs,] <- do.call(dlm.ImportanceWts, IWArgs)
	      if (CheckMax & any(z$Weights$log.h>0)) {
	        warning(paste("See warning from dlm.ImportanceWts: quitting",
                          "PostPredSample without\n   generating 'Obs'",
						  "component of result. 'Thetas' component has",
						  "\n  'BestID' attribute giving number of row",
						  "with highest log-posterior."))
	        attr(z$Thetas, "BestID") <- which.max(z$Weights$log.h)			
        	return(z)
          }
        }
        SSArgs$Thetas <- NewThetas
        NewStates <- do.call(SampleStates, SSArgs)
        z$States[FailedIDs,,] <- NewStates
        if (is.null(attr(NewStates, "FailedIDs"))) {
          FailedIDs <- NULL
        } else {
          FailedIDs <- FailedIDs[attr(NewStates, "FailedIDs")]
        }
      }
    }
  }
  if (as.numeric(messages)>0) cat("Sampling observable series ...\n")
  SOArgs <- c(list(Thetas=z$Thetas, States=z$States, build=Build, 
                   Y=Data, WhichEls=WhichEls, ReplaceAll=ReplaceAll,
                   NonNeg=NonNeg), DotArgs)
  z$Obs <- do.call(SampleObs, SOArgs)
  if (Importance) {
    #
	#  NB dlm.ImportanceWts returns normalised weights: if there have been 
	#  several calls e.g. due to the use of FailedIDs then different 
	#  subsets will have been normalised separately. Option here is 
	#  either to call dlm.ImportanceWts again (costly likelihood
	#  calculations) or to just recalculate the weights themselves
	#  from the values of log.h and log.g already contained in z$Weights
	#
	z$Weights$w <- exp(z$Weights$log.h - z$Weights$log.g)
	z$Weights$w <- z$Weights$w / sum(z$Weights$w)
    SortedWts <- sort(z$Weights$w, decreasing=TRUE)
    NWts <- nrow(z$Weights)
    if (sum(SortedWts[1:ceiling(NWts/10)]) > 0.8*sum(SortedWts)) {
      warning(paste("Importance sampling may be unreliable here:",
                    "more than 80% of the\n  total weight is",
                    "concentrated in fewer than 10% of samples"))
    }
  }
  z$Model <- ModelBundle$Model
  z
}
######################################################################
PostPred.VidFrames <- 
  function(Data, Samples, AntiLog=FALSE, DatColours, DatTypes=c(1,1), 
           PredColours, Nsamples=20, Nsteps=5, replace=TRUE,
           Folder, FileNameRoot=Folder, WarnExisting=TRUE, 
           PNG.args=list(), MakeGif=TRUE, DelFrames=FALSE,
           Annotation=FALSE, ...) {
  #
  #   To produce the video frames for an animation of samples
  #   from a postprocessed ensemble using the approach of
  #   Bowman (2019). Arguments:
  #
  #   Data    Data frame containing time in column 1, observed
  #           series in column 2 and ensemble members in the 
  #           remaining columns. 
  #   Samples A list containing at least named components
  #           Obs and Weights - as produced by, for example,
  #           PostPredSample(). 
  #   AntiLog If TRUE, samples will be antilogged before plotting
  #           (appropriate when models were fitted on a log scale)
  #   DatColours, Vectors of colours and line types used to 
  #   DatTypes    plot the ensemble members: for details, 
  #               see header to PlotEnsTS(). 
  #   PredColours Vector of two colours: the first is the
  #               base colour for plotting objects relating 
  #               to the observed trend mu[0], and the 
  #               second for plotting the ensemble consensus.
  #               Only the first element is used here: 
  #               the second element is retained so that
  #               (for example) the same object can be
  #               used for this routine as for calls to 
  #               SmoothPlot().
  #   Nsamples    Number of "primary" samples to use to 
  #               construct the video frames
  #   Nsteps      Number of interpolation points between
  #               each primary sample
  #   replace     Logical scalar, controlling whether or not
  #               to draw primary samples with replacement.
  #   Folder      Folder within which to store the 
  #               output graphics files
  #   FileNameRoot Root of filename for output graphics files.
  #               E.g. if FileNameRoot is "fish" then the
  #               files will be named fish001.png, fish002.png, ... 
  #               Defaults to the name of Folder.
  #   WarnExisting  Logical scalar indicating whether to 
  #               warn if Folder exists already. 
  #   PNG.args    List of arguments to the png() command,
  #               used to control (e.g.) the size of the
  #               images produced.
  #   MakeGif     If TRUE, an animated GIF will also be 
  #               created.
  #   DelFrames   If TRUE and if MakeGif is TRUE, the 
  #               files containing the individual video 
  #               frames will be deleted after the Gif is
  #               created.
  #   Annotation  If TRUE, plots / frames will be annotated
  #               with copyright and funding information
  #   ...         Other arguments to SmoothPlot(). 
  #      
  #   The function produces a collection of graphics files in 
  #   Folder, each representing a successive frame of an 
  #   animated "walk through" the posterior predictive 
  #   distribution. It also optionally produces an animated
  #   GIF image combining the frames. The return value is a 
  #   matrix containing the time series used for each frame:
  #   the column names are Frame1.0, Frame1.1, ..., Frame2.0, etc.
  #
  Nt <- nrow(Data) # Number of time points
  Ns <- ncol(Samples$Obs) # Number of samples
  Nfiles <- Nsamples*Nsteps
  #
  #   Input checks
  #
  if (Nt != nrow(Samples$Obs)) stop("Data and Samples cover different time periods")
  if (is.null(Samples$Weights)){
    WtsToUse <- rep(1, Ns) 
  } else {
    if (nrow(Samples$Weights) != Ns) stop("Length of Weights doesn't match number of Samples")
    WtsToUse <- Samples$Weights$w
  }
  if (DelFrames & !MakeGif) {
    warning("DelFrames is TRUE but MakeGif is FALSE: DelFrames will be ignored")
  }
  #
  #   Set up output folder and file names
  #
  if (file.exists(Folder)) {
    if (WarnExisting) {
      Confirm <- 
        readline(paste("Folder",Folder,"already exists. Replace it? (Y/N) "))
      if (!(Confirm %in% c("Y", "y"))) stop("OK, processing terminated.", call.=FALSE)
      unlink(Folder, recursive=TRUE)
    }
  } else dir.create(Folder)
  FileNames <- paste(Folder, "/", FileNameRoot, 
                     formatC(1:Nfiles, width=1+floor(log10(Nfiles)),
                             flag="0"), ".png", sep="")
  #
  #   Select primary samples and create matrix with 
  #   interpolations. In case of NAs in the arguments, 
  #   discard these first with a warning. 
  #
  OKSamples <- !is.na(colMeans(Samples$Obs))
  ObsNonMiss <- Samples$Obs[, OKSamples]
  if (ncol(ObsNonMiss) != ncol(Samples$Obs)) {
    warning(paste(ncol(Samples$Obs)-sum(OKSamples), "of", ncol(Samples$Obs),
                  "posterior samples contained missing values: these have been discarded."))
    Ns <- ncol(ObsNonMiss)
    WtsToUse <- WtsToUse[OKSamples]
  }
  Skeleton <- ObsNonMiss[,sample(1:Ns, size=Nsamples, 
                                 replace=replace, prob=WtsToUse)]
  mu <- rowMeans(Skeleton)
  res <- matrix(nrow=Nt, ncol=Nsamples*Nsteps)
  colnames(res) <- paste("Frame", rep(1:Nsamples,each=Nsteps),
                         ".",rep(1:Nsteps, Nsamples))
  for (i in 1:Nsamples) {
    #
    # z1 and z2 are the current and next values in the "skeleton".
    # Definition of z2 ensures that it wraps back round at the
    # end.
    #
    z1 <- Skeleton[,i] ; z2 <- Skeleton[,1 + (i %% Nsamples)]
    for (j in 0:(Nsteps-1)) {
      cur.col <- (i-1)*Nsteps + j + 1
      a <- j/Nsteps
      #
      #   Next line is Adrian Bowman's beautiful trick
      #
      res[, cur.col] <- 
        mu + ( ((1-a)*(z1-mu)) + (a*(z2-mu)) ) / sqrt(a^2 + (1-a)^2)
    }
  }
  #
  #   Now write the output files
  #
  MissObs <- is.na(Data[,2])
  for (i in 1:Nfiles) {
    do.call(png, args=c(PNG.args, file=FileNames[i]))
    par(mar=c(3,3,3,1),mgp=c(2,0.75,0),lwd=2)
    SmoothPlot(Data, Samples=Samples, AntiLog=AntiLog,
               DatColours=DatColours, PredColours=PredColours, 
               DatTypes=DatTypes, PlotMu0=FALSE, ...)
    Y0 <- Data[,2]; Y0[MissObs] <- res[MissObs,i]
    if (AntiLog) Y0 <- exp(Y0)
    lines(Data[,1], Y0, lwd=5, col=DatColours[1])
    if (Annotation) GraphAnnot()
    dev.off()
  }
  #
  #   And make the GIF animation if requested
  #
  if (MakeGif) {
    GifFrames <- image_read(FileNames)
    GifFileName <- paste(Folder, "/", FileNameRoot, ".gif", sep="")
    image_write(image_animate(GifFrames, optimize=TRUE, fps=5), 
                path=GifFileName)
    if (DelFrames) file.remove(FileNames)
  }
  res
}
######################################################################
SLLT.modeldef <- function(theta, m0=NULL, C0=NULL, kappa=1e6, Shift=0) {
  #
  #  Sets up the structure of a "smooth local linear trend" state space
  #  time series model. Arguments: 
  #
  #  theta  Vector of length 2: theta[1] is the log of the "measurement 
  #         error" variance, and theta[2] is the log of the innovation 
  #         variance on the slope process. 
  #  kappa  Initialisation variance for diffuse elements of the
  #         state vector. Ignored if C0 is non-NULL.
  #  m0     Optional vector of initial values for the state vector. 
  #         If NULL (the default) this is determined automatically
  #         from the model structure. 
  #  C0     Optional covariance matrix for initialising the state
  #         vector. If NULL, like m0 this is determined from the
  #         model structure and from kappa. 
  #  Shift  Fits and smooths can be done on data Y-Shift to
  #         improve numerical stability: often, Shift will be
  #         the sample mean of the observations. The model 
  #         structure returned is for the shifted data, but
  #         with a "Shift" attribute containing the vectors
  #         that must be added to the m0 component and any
  #         state estimates to get back to the original scale. 
  #
  #  The function returns a list defining a dynamical linear model, 
  #  in the format required by routines in the dlm library. 
  #
  FF <- matrix(c(1,0),nrow=1)	# Needs to be a row vector
  GG <- matrix(c(1,0,1,1),nrow=2)
  V <- as.numeric(exp(theta[1])) # Strip away names, which cause dlm to crash. 
  W <- diag(c(0,exp(theta[2]))) # Underlying mean has no innovation
  init <- make.SSinits(GG,W,kappa=kappa)
  if (!is.null(m0)) {
    if (!isTRUE(length(m0)==nrow(GG))) {
      stop(paste("m0 should be a vector containing",nrow(GG),"values"))
    }
    init$m0 <- m0 - c(Shift, 0)
  }
  if (!is.null(C0)) { # Overwrite if the user supplied a value
    if (!isTRUE(all(dim(C0)==nrow(GG)))) {
      stop(paste("C0 should be a ",nrow(GG),"x",nrow(GG)," matrix",sep=""))
    }
    init$C0 <- C0
  }
  z <- dlm(c(list(FF=FF,GG=GG,V=V,W=W),init))
  attr(z, "Shift") <- c(Shift,0)
  z
}
######################################################################
SLLT.IniPar <- function(Y, method="arima", collapse=TRUE, Tiny=1e-6) {
  #
  #   To find initial values for the two variance parameters in a 
  #   smooth local linear trend model. Arguments: 
  #   
  #   Y       Either a vector or a matrix containing observed series. 
  #           If a matrix, the columns are considered to contain
  #           the series. 
  #   method  Either "arima" (the default) or "moments". "arima"
  #           exploits the fact that the SLLT model is stochastically
  #           equivalent to an ARIMA(0,2,2) process and obtains the
  #           required values using the conditional sum of squares - 
  #           ML algorithm inplemented by the arima() command. 
  #           "moments" uses a method of moments based on the variance
  #           and lag-1 covariance of the second differences. Note,
  #           however, that the method of moments is very inefficient
  #           for models containing a moving average component:
  #           its only advantage is that it's a bit quicker, although
  #           experience suggests that it can indeed produce wildly 
  #           inaccurate initial values.
  #   collapse  If Y is a matrix then, if collapse is TRUE, a single 
  #           pair of initial values is produced using the data from 
  #           all columns; if collapse is FALSE then values are
  #           produced for each column separately. 
  #   Tiny    A small value used to prevent negative variance
  #           estimates. 
  #
  #   Value: if Y is a vector or if collapse is TRUE, a vector of 
  #   length 2 containing initial values for the logarithms of the 
  #   measurement and innovation variances in the SLLT model. If Y 
  #   is a matrix and collapse is FALSE, a 2-row matrix containing
  #   separate pairs for each column of Y.  
  #
  if (method=="moments") {
    D2Y <- # as.matrix() creates 1-column matrix if Y was a vector
      diff(as.matrix(Y), differences = 2) 
    Moments <- apply(D2Y, MARGIN=2, 
                     FUN=function(x) {
                       acf(x, plot=FALSE, type="covariance", 
                           na.action=na.pass, lag.max=1)$acf
                     }
    )
    if (collapse) Moments <- as.matrix(rowMeans(Moments))
    theta.init <- matrix(nrow=2, ncol=ncol(Moments))
    theta.init[1,] <- pmax(-Moments[2,]/4, Tiny)
    theta.init[2,] <- pmax(Moments[1,]-(6*theta.init[1,]), Tiny)
  } else if (method=="arima") {
    Ymat <- as.matrix(Y) # creates 1-column matrix if Y was a vector
    ArimaPars <- apply(Ymat, MARGIN=2, 
                       FUN=function(x) {
                         z <- arima(x, order=c(0,2,2))
                         c(coef(z), z$sigma2)
                       }
                       )
    theta.init <- matrix(nrow=2, ncol=ncol(ArimaPars))
    theta.init[1,] <- # Equate lag-1 autocovariance of fitted models
      pmax(Tiny, -ArimaPars[1,]*(1+ArimaPars[2,])*ArimaPars[3,]/4)
    theta.init[2,] <- # And variances
      pmax(Tiny, 
           ((1+(ArimaPars[1,]^2)+(ArimaPars[2,]^2))*ArimaPars[3,]) -
             (6*theta.init[1,]))
    if (collapse) theta.init <- rowMeans(theta.init)
  } else stop("'method' must be either 'arima' or 'moments'")
  log(theta.init) # Log scale for optimisation
}
######################################################################
SLLTSmooth <- function(Y, m0=NULL, C0=NULL, kappa=1e6, prior.pars=NULL,  
                       theta=NULL, messages=TRUE, Use.dlm=FALSE, 
					   debug=FALSE, ...) {
  #
  #   Fits and applies a "smooth local linear trend" model to a univariate
  #   time series. Arguments:
  #
  #   Y           Vector containing time series data
  #   m0          Optional vector of initial values for the state vector. 
  #               If NULL (the default) this is determined automatically
  #               from the model structure. 
  #   kappa       Initialisation variance for diffuse elements of the 
  #               state vector
  #   prior.pars  Optional 2*2 matrix containing means and standard
  #               deviations of Gaussian prior distributions for the log 
  #               variance parameters: in this case, maximum a posterior
  #               (MAP) estimation is done. 
  #	  theta	      Optional initial value for parameter vector. If NULL,
  #               the routine will auto-initialise the optimisation
  #   messages    Controls printing of progress details to screen. See
  #               dlm.SafeMLE
  #   Use.dlm     If TRUE, fitting is done using dlmMLE to maximise the 
  #               log-likelihood. Otherwise nlm is used. In the latter 
  #               case (the default), MAP estimation can also be done 
  #               by specifying prior.pars (see above). 
  #   debug       Controls whether to use R or C code in the dlm library 
  #               (debug=FALSE means use C code, which is faster but
  #               doesn't allow tracing of problems)
  #  
  #  The function returns a list with three components: the first is the 
  #  object containing results of the maximum likelihood parameter 
  #  estimation, the second is the fitted model itself, and the third is 
  #  the result of applying the Kalman Smoother to the input series Y 
  #  using this fitted model.
  #
  #  Start by finding initial values for a numerical maximisation of the 
  #  likelihood: based on a method of moments on the second differences
  #  of the data.
  #
  par.names <- c("log(sigsq)", "log(tausq)")
  if (is.null(theta)) {
    theta.init <- as.numeric(SLLT.IniPar(Y))
  } else theta.init <- as.numeric(theta)
  names(theta.init) <- par.names
  #
  #  Now the estimation, printing information if required. For 
  #  numerical stability, the mean of the data is subtracted
  #  from both the data and trend, and then added back later.
  #
  if (as.numeric(messages)>0) {
    cat("Estimating model parameters - please be patient ...\n")
  }
  Ybar <- mean(Y, na.rm=TRUE)
  thetahat <- 
    dlm.SafeMLE(theta.init, Y-Ybar, SLLT.modeldef,
                debug=debug, Use.dlm=Use.dlm, par.names=par.names, 
                prior.pars=prior.pars, m0=m0, C0=C0, kappa=kappa,
                Shift=Ybar, messages=messages, hessian=TRUE, ...)
  Model <- SLLT.modeldef(thetahat$par, m0=m0, C0=C0,     # Original scale
                         kappa=kappa)
  ShiftMod <- SLLT.modeldef(thetahat$par, m0=m0, C0=C0, # Shifted
                            kappa=kappa, Shift=Ybar)
  if (as.numeric(messages)>0) {
    cat("\nSUMMARY OF STATE SPACE MODEL:\n")
    print(summary.dlmMLE(thetahat))
  }
  #
  #  Kalman Smoother, and return result
  #
  Smooth <- tryCatch(dlmSmooth(Y-Ybar, ShiftMod, debug=debug), 
                     error=function(e) { 
                       z <- rep(NA, length=Y)
                       attr(z, "error") <- e; z
                       }
                     )
  if (!is.null(attr(Smooth, "error"))) {
    warning(paste("dlmSmooth failed - error message was\n ", 
                  attr(Smooth, "error")$message), immediate.=TRUE)
  } else {
    Smooth$s <-  # Add the shift back
      t(t(Smooth$s) + attr(ShiftMod, "Shift"))
  }
  list(Data=list(Y=Y), Theta=thetahat, Model=Model, Smooth=Smooth)
}
######################################################################
EnsSLLT.modeldef <- function(theta, m0=NULL, C0=NULL, kappa=1e6,
                             Shift=0, NEnsTS, Groups=NULL, RunsPerTS=NULL,
                             discrepancy="varying", UseAlpha=TRUE,
                             constrain=TRUE) {
  #
  #   Sets up the structure of a simple "Ensemble smooth local 
  #   linear trend" state space time series model for an observed
  #   time series and exchangeable ensemble. The notation for the model 
  #   components corresponds to that used in the dlm library, with an
  #   additional component "Consensus" giving the coefficients needed
  #   to create an estimate of the ensemble consensus. If discrepancy
  #   is "varying" and UseAlpha is TRUE (see below) then the 
  #   the parameter vector theta is defined as follows:
  #
  #  theta[1]    This is alpha, the scaling on the relationship between
  #              the ensemble consensus and observed slope processes. 
  #  theta[2]    This is log(sigsq[0]), where sigsq[0] is the "measurement
  #              error" variance in the observed series. 
  #  theta[3]    This is log(tausq[0]), where tausq[0] is the slope 
  #              innovation variance in the observed series. 
  #  theta[4]    This is log(tausq[death]), where tausq[death] is the slope
  #              innovation variance for the shared discrepancy between 
  #              the ensemble trend and the real trend. 
  #  theta[5]    This is log(sigsq[1]), where sigsq[1] is the common 
  #              "measurement error" variance in the ensemble members.
  #  theta[6]    This is log(tausq[1]), where tausq[1] is the common slope
  #              innovation variance for the ensemble members about the
  #              shared discrepancy.
  #
  #  If discrepancy is "constant" and UseAlpha is FALSE then theta is a 
  #  4-element vector containing alpha, log(sigsq[0]), log(tausq[0]) and 
  #  log(sigsq[1]): in this case, tausq[death] and tausq[1] are both 
  #  taken as zero.
  #
  #  If UseAlpha is FALSE then alpha is fixed to 1 and omitted from 
  #  theta - so that in this case theta[1] is log(sigsq[0]), theta[2]
  #  is log(tausq[0]) etc. This argument makes no difference to the
  #  results if discrepancy is "constant". 
  #
  #  The function contains the following additional arguments:
  #
  #   m0         Optional vector of initial values for the state vector. 
  #              If NULL (the default) this is determined automatically
  #              from the model structure. 
  #   C0         Optional covariance matrix for initialising the state
  #              vector. If NULL, like m0 this is determined from the
  #              model structure and from kappa. 
  #   kappa      "Large" value with which to initialise variances of
  #              diffuse components of state vector at time 0. Taking
  #              it too large can lead to numerical instability, 
  #              however. It should be several orders of magnitude 
  #              larger than any of the model variances. This is
  #              ignored if C0 is non-NULL.
  #   Shift      Fits and smooths can be done on data Y-Shift to
  #              improve numerical stability: often, Shift will be
  #              the sample mean of the observations. The model 
  #              structure returned is for the shifted data, but
  #              with a "Shift" attribute containing the vectors
  #              that must be added to the m0 component and any
  #              state estimates to get back to the original scale. 
  #   NEnsTS     Number of series in ensemble (often # of members)
  #   Groups     Vector of length NEnsTS, indicating group membership. 
  #              This should be used if an ensemble contains multiple 
  #              runs from one or more simulators: in this case the 
  #              elements of Groups should be integers between 1 and G
  #              where G is the total number of simulators contributing
  #              to the ensemble: Groups[i] is the index of the 
  #              simulator that was used to produce the ith ensemble
  #              member. If NULL, each member is assumed to come from
  #              a different simulator. 
  #   discrepancy Either "varying" (the default) or "constant". If 
  #              "constant" then the discrepancies between each 
  #              ensemble member and reality are constant through 
  #              time; otherwise they evolve as random walks. 
  #   UseAlpha   Logical scalar indicating whether or not to include
  #              a scaling factor to represent the expected value of
  #              the ensemble consensus trend slope beta[d] given the 
  #              observed one beta[0]. If FALSE, the expected ensemble
  #              consensus trend is equal to the observed one. 
  #   constrain  Logical scalar indicating whether or not to constrain
  #              the discrepancy trends and slopes for the individual
  #              ensemble members to sum to zero. Without this 
  #              constraint, the model is formally not identifiable -
  #              although linear combinations of elements of the 
  #              state vector are. 
  #
  if (is.null(Groups)) Groups <- 1:NEnsTS
  if (UseAlpha) alpha <- theta[1] else alpha <- 1
  IdxShift <- 1-as.numeric(UseAlpha) # Shifts indices in theta
  sigsq.0 <- exp(theta[2-IdxShift])
  tausq.0 <- exp(theta[3-IdxShift])
  if (discrepancy=="varying") {
    tausq.d <- exp(theta[4-IdxShift])
    sigsq.1 <- exp(theta[5-IdxShift])
    tausq.1 <- exp(theta[6-IdxShift])
  } else if (discrepancy=="constant") {
    tausq.d <- 0
    sigsq.1 <- exp(theta[4-IdxShift])
    tausq.1 <- 0
  } else {
    stop("discrepancy must be either 'constant' or 'varying'")
  }
  NGroups <- max(Groups)
  nS <- 2*(NGroups+2) # Length of state vector
  FF <- matrix(0, nrow=NEnsTS+1, ncol=nS)
  FF[1,1] <- 1; FF[-1,1] <- alpha; FF[-1,3] <- 1
  #
  # Next lines use matrix subsetting to pick out the correct columns
  # for each simulator in Groups.
  #
  FF[cbind(2:nrow(FF), (2*Groups)+3)] <- 1
  if (constrain) {
    FF[1+which(Groups==max(Groups)), seq(5, nS-3, 2)] <- -1 # Sum-to-zero constraint on final simulator
  }
  if (constrain) FF[NEnsTS+1, seq(5, nS-3, 2)] <- -1 # Sum-to-zero constraint on final row
  GG <- diag(rep(1,nS))
  GG[(col(GG)%%2==0) & (col(GG)-row(GG)==1)] <- 1
  V <- c(sigsq.0, rep(sigsq.1, NEnsTS))
  if (!is.null(RunsPerTS)) { # Scale variances for averaged series
    if (!isTRUE(length(RunsPerTS)==NEnsTS)) {
      stop("RunsPerTS should be a vector of length NEnsTS")
    }
    V[-1] <- V[-1] / RunsPerTS
  }
  V <- diag(V)
  W <- diag(rep(0,nS)) 
  if (constrain) { # W is block diagonal if constraints are used
    diag(W)[c(2,4)] <- c(tausq.0, tausq.d)
    tau.els <- seq(6, nS-2, 2)
    W[tau.els, tau.els] <- -tausq.1 / NGroups
    diag(W)[tau.els] <- diag(W)[tau.els] + tausq.1
    nS <- nS-2 # Remove redundant elements of state vector 
    FF <- FF[,1:nS]
    GG <- GG[1:nS, 1:nS]
    W <- W[1:nS, 1:nS]
  } else { # Unconstrained version: W is diagonal
    diag(W)[2*(1:(NGroups+2))] <- c(tausq.0, tausq.d, rep(tausq.1, NGroups))
  }
  init <- make.SSinits(GG,W,kappa=kappa)
  if (!is.null(m0)) {
    if (!isTRUE(length(m0)==nrow(GG))) {
      stop(paste("m0 should be a vector containing",nrow(GG),"values"))
    }
    init$m0 <- m0 - (Shift*c(1, 0, 1-alpha, rep(0, nS-3))) # For shifted data
  }
  if (constrain) { # For constrained model, set C0 to reflect sum-to-zero constraints as well
    init$C0[5:nS, 5:nS] <- -kappa / NEnsTS
    diag(init$C0)[5:nS] <- diag(init$C0)[5:nS] + kappa
  }
  if (!is.null(C0)) { # Overwrite if the user supplied a value
    if (!isTRUE(all(dim(C0)==nrow(GG)))) {
      stop(paste("C0 should be a ",nrow(GG),"x",nrow(GG)," matrix",sep=""))
    }
    init$C0 <- C0
  }
  nS <- ncol(FF) # dimension of state vector after eliminating unwanted rows
  Consensus <- rep(0,nS); Consensus[c(1,3)] <- c(alpha,1) # Coefficients yielding ensemble consensus
  #
  #   Ideally, would return a dlm object. However, the dlm() routine
  #   checks for positive definiteness of W and C0 via a call to
  #   eigen() which can return tiny negative eigenvalues for the
  #   structure of this model: this creates an error condition.
  #   Just return a list, therefore - the rest of the dlm 
  #   routines can cope with this. 
  #
  z <- list(FF=FF,GG=GG,V=V,W=W,m0=init$m0,C0=init$C0, Consensus=Consensus)
  attr(z, "Shift") <- as.numeric(Shift*c(1, 0, 1-alpha, rep(0, nS-3)))
  z
}
######################################################################
EnsSLLTSmooth <- 
  function(Y, m0=NULL, C0=NULL, kappa=1e6, discrepancy="varying", 
           Groups=NULL, compact=!is.null(Groups), UseAlpha=TRUE, 
           prior.pars=NULL, theta=NULL, constrain=TRUE, Tiny=1/kappa,
           ObsSmooth, Ens0Theta, messages=TRUE, Use.dlm=FALSE, 
           debug=FALSE, ...) {
  #
  #  Fits and applies a model of the form defined by EnsSLLT.modeldef or
  #  EnsSLLT0.modeldef. The arguments are:
  #
  #   Y:          Matrix containing an observed time series in its first 
  #               column and exchangeable series from an ensemble in the 
  #               remaining columns. 
  #   m0          Optional vector of initial values for the state vector. 
  #               If NULL (the default) this is determined automatically
  #               from the model structure. 
  #   kappa:      Value used to initialise variances for diffuse elements 
  #               of the state vector. 
  #   discrepancy If "varying" (the default), the discrepancy term between 
  #               the observed and ensemble trends is time-varying; for a 
  #               constant discrepancy, set this to "constant". 
  #   Groups      Vector indicating group membership for each 
  #               ensemble member. This should be used if an ensemble
  #               contains multiple runs from one or more simulators: 
  #               in this case the elements of Groups should be integers
  #               between 1 and G where G is the total number of 
  #               simulators contributing to the ensemble: Groups[i] is 
  #               the index of the simulator that was used to produce 
  #               the ith ensemble member. If NULL, each member is 
  #               assumed to come from a different simulator. 
  #   compact     Logical scalar indicating whether to compact the
  #               ensemble using CompactEns() before fitting. For 
  #               ensembles with many members per simulator, this
  #               may speed up the fitting without sacrificing 
  #               information because the simulator-specific means
  #               are sufficient statistics. Defaults to TRUE if
  #               Groups is non-NULL and FALSE otherwise (because 
  #               CompactEns() just returns the original ensemble
  #               in that case). 
  #   UseAlpha    Logical scalar indicating whether or not to include
  #               a scaling factor to represent the expected value of
  #               the ensemble consensus trend slope beta[d] given the 
  #               observed one beta[0]. If FALSE, the expected ensemble
  #               consensus trend is equal to the observed one. 
  #   prior.pars  An optional 3*2 (if discrepancy is "constant") or
  #               5*2 (if discrepancy is "varying") matrix containing
  #               means and standard deviations of independent Gaussian
  #               priors for the log variance parameters in the model. If
  #               this is provided, maximum a posteriori (or penalised
  #               maximum likelihood) estimation is used; otherwise
  #               just normal maximum likelihood. 
  #	  theta	      Optional initial value for parameter vector. If NULL,
  #               the routine will auto-initialise the optimisation
  #   Tiny        A small positive value, used to replace negative initial 
  #               estimates for variances
  #   ObsSmooth   Optional: result of a previous call to SLLTSmooth() for
  #               the first column of Y. If present, this is used to help 
  #               initialise the numerical search for an MLE if a "constant" 
  #               model needs to be fitted, or if the resulting smooth is 
  #               needed. Otherwise the relevant model will be fitted as 
  #               part of this routine. 
  #   Ens0Theta   Optional parameter vector for a previous fit of
  #               an "ensemble local linear trend" model with a constant
  #               discrepancy term. If present, this will be used to 
  #               help initialise the numerical search for an MLE if
  #               a "varying" model needs to be fitted here. Otherwise
  #               the relevant constant-discrepancy model will be fitted
  #               as part of this routine. 
  #   messages    Logical value: if TRUE then progress is printed to 
  #               screen. 
  #   Use.dlm   Logical: if TRUE then fitting is done using dlmMLE
  #               to optimise the likelihood, otherwise nlm is used.
  #   debug     Logical: argument passed to dlm routines to control
  #             whether to use faster C code (the default), or 
  #             slower R code which allows problem tracing. 
  #	  ...       Additional arguments to dlm.SafeMLE
  #
  #  The function returns a list with three components: the first is the 
  #  object containing results of the maximum likelihood parameter estimation, 
  #  the second is the fitted model itself, and the third is the result 
  #  of applying the Kalman Smoother to the input series using this fitted model. 
  #
  #  Start by checking the value of the "discrepancy" argument.
  #
  if (!(discrepancy %in% c("varying", "constant"))) {
    stop("discrepancy must be either 'constant' or 'varying'")
  }
  #
  #   Check for partially missing rows of the ensemble data
  #
  if (!NA.AllOr0(Y[,-1])) {
    warning(paste("Some but not all ensemble members contain missing",
    "value at some time points.\nThis may cause problems when fitting",
    "because the dlm library routines do not cope well\nwith missing",
    "values in the observation vector unless the corresponding",
    "measurement errors\nare uncorrelated."))
  }
  #
  #  Compact the ensemble if requested (NB Y and Groups here are
  #  local copies in the active environment, hence re-using them
  #  won't affect anything else).
  #
  if (!is.null(Groups)) {
    if (compact) {
      Y <- CompactEns(Y, Groups)
      Groups <- NULL # No need for groups now
      message("######\n",
              "######  NOTE: compacting Y and Groups. Do not use these objects",
              "\n######  when working with the results of this call: use the ",
              "\n######  versions in the `Data` component of the result.",
              "\n######")
    }
  } else {
    if (compact) warning("'compact' argument does nothing when Groups is NULL")
  }
  RunsPerTS <- attr(Y, "RunsPerTS") # Could be NULL
  NEnsTS <- ncol(Y) - 1 # Number of series in ensemble
  EnsembleMean <- EnsMean(Y, Groups)
  Ybar <- mean(Y[,1], na.rm=TRUE)
  #
  #  Produce a smooth of the observations if it hasn't been supplied.
  #  NB if m0 or C0 is NULL then the subsets are too (handy).
  #
  if (missing("ObsSmooth")) 
    ObsSmooth <- SLLTSmooth(Y[,1], messages=FALSE, Use.dlm=Use.dlm, 
                            m0=m0[1:2], C0=C0[1:2, 1:2],
                            kappa=kappa, debug=debug, ...)
  #
  #   Get approximate values of SLLT fits to each series in the 
  #   ensemble - for later use. Adjustment in second line is to 
  #   account for the possibility that Y may be compacted at 
  #   this point: the parameters are the log variances in the 
  #   measurement equation, so the adjustment gets back to the
  #   log variance before compacting by calculating group means. 
  #
  ApproxEnsPars <- SLLT.IniPar(Y[,-1], collapse=FALSE)
  if (!is.null(RunsPerTS)) {
    ApproxEnsPars[1,] <- ApproxEnsPars[1,] + log(RunsPerTS)
  }
  #
  #  Now fit a constant-discrepancy model if this is needed, which is the 
  #  case EITHER if discrepancy=="constant" OR if discrepancy=="varying"
  #  and Ens0Theta hasn't been provided. 
  #  
  if (as.numeric(messages)>0) {
    cat("Estimating model parameters - please be patient ...\n")
  }
  if (discrepancy=="constant" | 
      (discrepancy=="varying" & is.null(theta) & missing("Ens0Theta"))) {
    #
    #   For the constant-discrepancy version and if UseAlpha is TRUE,
    #   initial value for theta[1] is 1 (implying that ensemble 
    #   consensus slopes are more or less proportional to those for
    #   observations). Initial values for theta[2] and theta[3] are 
    #   the MLEs from an observation-only model, which needs to be 
    #   fitted if not supplied as an argument. The smooth also needs
    #   to be computed for use with a varying-discrepancy model 
    #   later: this isn't needed for the constant-discrepancy case, 
    #   but it's not too expensive and the code is less fiddly. For
    #   theta[4], use "auto-initialise" for each ensemble member,
    #   then convert back to the scale of individual runs and average.
    #
    par.names <- c("alpha", "log(sigsq[0])", "log(tausq[0])", "log(sigsq[1])")
    if (!UseAlpha) par.names <- par.names[-1]
    if (!is.null(theta)) {
      theta.init <- theta 
    } else {
      sigsq1.init <- log(mean(exp(ApproxEnsPars[1,])))
      theta.init <- c(1, ObsSmooth$Theta$par, sigsq1.init)
      if (!UseAlpha) theta.init <- theta.init[-1]
    }
    names(theta.init) <- par.names
    CurPriors <- prior.pars # Need to amend for this fit if discrepancy is "varying"
    Npars <- nrow(CurPriors)
    if (discrepancy=="varying") {
      CurPriors <- CurPriors[-c(Npars-2, Npars),]
    }
    thetahat <- 
      dlm.SafeMLE(theta.init, Y-Ybar, EnsSLLT.modeldef, m0=m0, C0=C0, 
                  kappa=kappa, Shift=Ybar, NEnsTS=NEnsTS, Groups=Groups, 
                  discrepancy="constant", RunsPerTS=RunsPerTS, 
                  UseAlpha=UseAlpha, prior.pars=CurPriors, 
                  constrain=constrain, hessian=(discrepancy=="constant"),
                  par.names=par.names, Use.dlm=Use.dlm, 
                  messages=(messages & discrepancy=="constant"), 
                  debug=debug, ...)
  }
  if (discrepancy=="varying") {
    par.names <- c("alpha", "log(sigsq[0])", "log(tausq[0])", 
                   "log(tausq[w])", "log(sigsq[1])", "log(tausq[1])")
    if (!UseAlpha) par.names <- par.names[-1]
    if (!is.null(theta)) {
      theta.init <- theta
    } else {
      #
      #   For the varying-discrepancy version, three of the initial
      #   values can be taken from a constant-discrepancy fit - which was
      #   either provided as Ens0Theta or calculated above.
      #
      theta.init <- rep(NA, 6)
      if (UseAlpha) {
        theta.init[c(1:3,5)] <- if (missing("Ens0Theta")) thetahat$par else Ens0Theta
      } else {
        theta.init[c(2:3,5)] <- if (missing("Ens0Theta")) thetahat$par else Ens0Theta
      }
      #
      #  For theta[4], take the difference between the ensemble mean and 
      #  observation-based smooth as an estimate of the discrepancy trend 
      #  series (don't forget to remove the first element of the 
      #  observation-based smooth corresponding to time point 0!). 
      #  theta[4] is the log variance of the second differences of this
      #  discrepancy series.
      #
      EnsembleMean <- EnsMean(Y, Groups)
      if (UseAlpha) {
        Death.hat <- 
          EnsembleMean - (theta.init[1] * ObsSmooth$Smooth$s[-1,1])
      } else Death.hat <- EnsembleMean-ObsSmooth$Smooth$s[-1,1]
      IdxShift <- 1-as.numeric(UseAlpha) # Shifts indices in theta
      theta.init[4] <- 
        log(var(diff(Death.hat, differences=2), na.rm=TRUE))
      #
      #  Initialise theta[6] to something very small, on the grounds 
      #  that the existing constant-discrepancy fit corresponds to 
      #  an inter-member variance of zero. By starting with a very 
      #  small value therefore, hopefully we're already on a ridge
      #  of the log posterior which will make it easier to locate
      #  a global optimum. 
      #
      theta.init[6] <- log(1e-12)
      if (!UseAlpha) theta.init <- theta.init[-1]
    }
    names(theta.init) <- par.names
    thetahat <- 
      dlm.SafeMLE(theta.init, Y-Ybar, EnsSLLT.modeldef, m0=m0, C0=C0, 
                  kappa=kappa, Shift=Ybar, NEnsTS=NEnsTS, Groups=Groups, 
                  RunsPerTS=RunsPerTS, discrepancy="varying", 
                  UseAlpha=UseAlpha, prior.pars=prior.pars, 
                  constrain=constrain, hessian=TRUE, messages=messages, 
                  par.names=par.names, Use.dlm=Use.dlm, debug=debug, ...)
  }
  Model <- # Original scale
    EnsSLLT.modeldef(thetahat$par, m0=m0, C0=C0, kappa=kappa, 
                     NEnsTS=NEnsTS, Groups=Groups, RunsPerTS=RunsPerTS,
                     discrepancy=discrepancy, UseAlpha=UseAlpha, 
                     constrain=constrain)
  ShiftMod <- # 
    EnsSLLT.modeldef(thetahat$par, m0=m0, C0=C0, kappa=kappa, Shift=Ybar,
                     NEnsTS=NEnsTS, Groups=Groups, RunsPerTS=RunsPerTS, 
                     discrepancy=discrepancy, UseAlpha=UseAlpha, 
                     constrain=constrain)
  if (as.numeric(messages)>0) {
    cat("\nSUMMARY OF STATE SPACE MODEL:\n")
    print(summary.dlmMLE(thetahat))
  }
  #
  #  Kalman Smoother, and return result
  #
  Smooth <- tryCatch(dlmSmooth(Y-Ybar, ShiftMod, debug=debug), 
                     error=function(e) { 
                       z <- rep(NA, length=Y)
                       attr(z, "error") <- e; z
                     }
  )
  if (!is.null(attr(Smooth, "error"))) {
    warning(paste("dlmSmooth failed - error message was\n ", 
                  attr(Smooth, "error")$message), immediate.=TRUE)
  } else {
    Smooth$s <-  # Add the shift back
      t(t(Smooth$s) + attr(ShiftMod, "Shift"))
  }
  list(Data=list(Y=Y, Groups=Groups),
       Theta=thetahat, Model=Model, Smooth=Smooth)
}
######################################################################
EBMtrend.modeldef <- function(theta, Xt, m0=NULL, C0=NULL, 
                              kappa=1e6, Shift=0, UsePhi=TRUE) {
  #
  #  Sets up the structure of a state space model inspired by 
  #  the structure of a simple energy balance model (EBM). Arguments:
  #
  #  theta    Parameter vector. If UsePhi is TRUE (see below) then
  #           this contains three elements: the log of the "measurement
  #           "error" variance, the log of the "trend drift" innovation
  #           variance, and the logit of the dependence parameter. If
  #           UsePhi is FALSE then theta contains just the first two
  #           of these elements.
  #  Xt       Vector containing a time series of effective forcings. 
  #  m0       Optional vector of initial values for the state vector. 
  #           If NULL (the default) this is determined automatically
  #           from the model structure. 
  #  C0       Optional covariance matrix for initialising the state
  #           vector. If NULL, like m0 this is determined from the
  #           model structure and from kappa. 
  #  kappa    "Large" value giving the variance used to initialise diffuse 
  #           elements of the state vector. Ignored if C0 is non-NULL.
  #  Shift    Fits and smooths can be done on data Y-Shift to
  #           improve numerical stability: often, Shift will be
  #           the sample mean of the observations. The model 
  #           structure returned is for the shifted data, but
  #           with a "Shift" attribute containing the vectors
  #           that must be added to the m0 component and any
  #           state estimates to get back to the original scale. 
  #  UsePhi   Controls whether or not to include a thermal inertia 
  #           parameter in the model structure. 
  #
  phi <- 0
  if (UsePhi) phi <- 1 / (1+exp(-theta[3])) # Inverse logit transform
  FF <- matrix(c(1,0,0),nrow=1)	# Needs to be a row vector
  GG <- matrix(c(phi, 1, -abs(kappa), # Use -abs(kappa) to mark positions of time-varying elements
                 0, 1, 0,
                 0, 0, 1), nrow=3, byrow=TRUE)
  V <- as.matrix(as.numeric(exp(theta[1]))) # Strip away names, which cause dlm to crash. 
  W <- diag(c(0,exp(theta[2]),0))# Only innovation is in offset process
  JGG <- diag(rep(0,3)); JGG[1,3] <- 1 # (1,3) element of GG is from column 1 of X
  init <- make.SSinits(GG,W,kappa=kappa) # Time-invariant initialisation should be OK
  if (!is.null(m0)) {
    if (!isTRUE(length(m0)==nrow(GG))) {
      stop(paste("m0 should be a vector containing",nrow(GG),"values"))
    }
    init$m0 <- m0 - (Shift * c(1, 1-phi, 0))
  }
  if (!is.null(C0)) { # Overwrite if the user supplied a value
    if (!isTRUE(all(dim(C0)==nrow(GG)))) {
      stop(paste("C0 should be a ",nrow(GG),"x",nrow(GG)," matrix",sep=""))
    }
    init$C0 <- C0
  }
  z <- list(FF=FF, GG=GG, V=V, W=W, JGG=JGG, X=as.matrix(Xt),m0=init$m0,C0=init$C0)
  attr(z, "Shift") <- as.numeric(Shift*c(1, 1-phi, 0))
  z
}
######################################################################
EBMtrendSmooth <- function(Y, Xt, m0=NULL, C0=NULL, kappa=1e6, UsePhi=TRUE, 
                           prior.pars=NULL, theta=NULL, messages=TRUE, 
                           Use.dlm=FALSE, debug=FALSE, ...) {
  #
  #   Fits and applies a model of the form defined by EBMtrend.modeldef.
  #   Arguments:
  #
  #   Y           Vector containing a univariate time series.
  #   Xt          Vector containing the corresponding forcings. 
  #   m0          Optional vector of initial values for the state vector. 
  #               If NULL (the default) this is determined automatically
  #               from the model structure. 
  #   kappa       Initialisation variance for diffuse elements of the
  #               state vector. 
  #   UsePhi      Controls whether or not to use the full model including 
  #               a thermal inertia parameter. 
  #   prior.pars  Optional 3*2 matrix. If present, it is taken to 
  #               contain the means and standard deviations of 
  #               independent Gaussian priors for the transformed 
  #               model parameters (log variances, and logit 
  #               thermal inertia). 
  #	  theta	      Optional initial value for parameter vector. If NULL,
  #               the routine will auto-initialise the optimisation
  #   messages    Controls whether to print progress to screen. 
  #   ...		  Additional arguments to dlm.SafeMLE
  #
  #   The function returns a list with three components: the first is 
  #   the object containing parameter estimation results, the second is 
  #   the fitted model itself, and the third is the result of applying 
  #   the Kalman Smoother to the input series using this fitted model. 
  #
  #  Start by finding initial values for a numerical maximisation of the 
  #  likelihood, by fitting an initial local linear trend model to Y. Use 
  #  this to get an initial log-likelihood so as to scale the optimisation.
  #
  par.names=c("log(sigsq)", "log(tausq)", "logit(phi)")
  if (!UsePhi) par.names <- par.names[-3]
  if (!is.null(theta)) {
    theta.init <- as.numeric(theta)
  } else {
    if (as.numeric(messages)>0) {
      cat("Carrying out initial local linear trend estimation ...\n")
    }
    SLLTModel0 <- 
      SLLTSmooth(Y, messages=FALSE, Use.dlm=Use.dlm, m0=m0[1:2], 
                 C0=C0[1:2, 1:2], kappa=kappa, debug=debug, ...)
    theta.init <- SLLTModel0$Theta$par
    if (UsePhi) {
      theta.init <- c(theta.init, 5) # the "5" is logit(0.99), roughly: thus initialising phi at ~1
    }
  }
  names(theta.init) <- par.names
  #
  #  Now the estimation, printing information if required. Shift the
  #  data and trend for numerical stability in the fitting. 
  #
  if (as.numeric(messages)>0) {
    cat("Estimating model parameters - please be patient ...\n")
  }
  Ybar <- mean(Y, na.rm=TRUE)
  thetahat <- 
    dlm.SafeMLE(theta.init, Y-Ybar, EBMtrend.modeldef, Xt=Xt, m0=m0, 
                C0=C0, kappa=kappa, Shift=Ybar, prior.pars=prior.pars, 
                UsePhi=UsePhi, hessian=TRUE, par.names=par.names, 
                Use.dlm=Use.dlm, messages=messages, debug=debug, ...)
  Model <- EBMtrend.modeldef(thetahat$par, Xt=Xt, m0=m0, C0=C0, 
                             kappa=kappa, UsePhi=UsePhi) # Original scale
  ShiftMod <- EBMtrend.modeldef(thetahat$par, Xt=Xt, m0=m0, C0=C0, 
                                kappa=kappa, Shift=Ybar, UsePhi=UsePhi) # Shifted
  if (as.numeric(messages)>0) {
    cat("\nSUMMARY OF STATE SPACE MODEL:\n")
    print(summary.dlmMLE(thetahat))
  }
  #
  #  Kalman Smoother, and return result
  #
  Smooth <- tryCatch(dlmSmooth(Y-Ybar, ShiftMod, debug=debug), 
                     error=function(e) { 
                       z <- rep(NA, length=Y)
                       attr(z, "error") <- e; z
                     }
  )
  if (!is.null(attr(Smooth, "error"))) {
    warning(paste("dlmSmooth failed - error message was\n ", 
                  attr(Smooth, "error")$message), immediate.=TRUE)
  } else {
    Smooth$s <-  # Add the shift back
      t(t(Smooth$s) + attr(ShiftMod, "Shift"))
  }
  list(Data=list(Y=Y), Theta=thetahat, Model=Model, Smooth=Smooth)
}
######################################################################
EnsEBMtrend.modeldef <- 
  function(theta, Xt, m0=NULL, C0=NULL, kappa=1e6, Shift=0, 
           NEnsTS, Groups=NULL, RunsPerTS=NULL, UseAlpha=TRUE, 
           UsePhi=TRUE, constrain=TRUE) {
  #
  #  Sets up the structure of a state space model for an observed
  #  time series and ensemble from exchangeable simulators, with trend 
  #  formulation inspired by the form of a simple energy balance model. 
  #  (EBM). The vector Xt contains a time series of effective forcings.
  #  If UseAlpha and UsePhi are both TRUE (see below) then the
  #  parameter vector theta is defined as follows:
  #
  #  theta[1]    This is alpha, the scaling on the relationship between
  #              the ensemble consensus and observed drift processes. 
  #  theta[2]    This is log(sigsq[0]), where sigsq[0] is the "measurement
  #              error" variance in the observed series. 
  #  theta[3]    This is log(tausq[0]), where tausq[0] is the slope 
  #              innovation variance in the observed series. 
  #  theta[4]    This is log(tausq[death]), where tausq[death] is the slope
  #              innovation variance for the shared discrepancy between 
  #              the ensemble drift and the real drift. 
  #  theta[5]    This is log(sigsq[1]), where sigsq[1] is the common 
  #              "measurement error" variance in the ensemble members.
  #  theta[6]    This is log(tausq[1]), where tausq[1] is the common drift
  #              innovation variance for the ensemble members about the
  #              shared discrepancy.
  #  theta[7]    This is logit(phi[0]), where phi[0] is the temporal dependence
  #              parameter in the mean component of the trend in the
  #              observations. 
  #  theta[8]    This is logit(phi[1]), where phi[1] is the common temporal
  #              dependence parameter in the mean component of the trends 
  #              in both the observations and the ensemble members.
  #
  #  If UseAlpha is FALSE then the value of alpha is fixed at 1 and is 
  #  omitted from theta so that theta[1] is log(sigsq[0]), theta[2] is
  #  log(tausq[0]) etc. If UsePhi is FALSE, then both of the final 
  #  elements of theta are omitted and the values of both phi[0] and
  #  phi[1] are set to 0. 
  #
  #  The function contains the following additional arguments:
  #
  #   m0         Optional vector of initial values for the state vector. 
  #              If NULL (the default) this is determined automatically
  #              from the model structure. 
  #   C0         Optional covariance matrix for initialising the state
  #              vector. If NULL, like m0 this is determined from the
  #              model structure and from kappa. 
  #   kappa      "Large" value with which to initialise variances of
  #              diffuse components of state vector at time 0. Taking
  #              it too large can lead to numerical instability, 
  #              however. It should be several orders of magnitude 
  #              larger than any of the model variances. This is 
  #              ignored if C0 is non-NULL.
  #   Shift      Fits and smooths can be done on data Y-Shift to
  #              improve numerical stability: often, Shift will be
  #              the sample mean of the observations. The model 
  #              structure returned is for the shifted data, but
  #              with a "Shift" attribute containing the vectors
  #              that must be added to the m0 component and any
  #              state estimates to get back to the original scale. 
  #   NEnsTS     Number of series in ensemble (often # of members)
  #   Groups     Vector of length NEnsTS, indicating group membership. 
  #              This should be used if an ensemble contains multiple 
  #              runs from one or more simulators: in this case the 
  #              elements of Groups should be integers between 1 and G
  #              where G is the total number of simulators contributing
  #              to the ensemble: Groups[i] is the index of the 
  #              simulator that was used to produce the ith ensemble
  #              member. If NULL, each member is assumed to come from
  #              a different simulator. 
  #   UseAlpha   Logical scalar indicating whether or not to include
  #              a scaling factor to represent the expected value of
  #              the ensemble consensus trend slope beta[d] given the 
  #              observed one beta[0]. If FALSE, the expected ensemble
  #              consensus trend is equal to the observed one. 
  #   UsePhi     Controls whether or not to include thermal inertia 
  #              parameters in the model structure. 
  #   constrain  Logical scalar indicating whether or not to constrain
  #              the discrepancy trends and slopes for the individual
  #              ensemble members to sum to zero. Without this 
  #              constraint, the model is formally not identifiable -
  #              although linear combinations of elements of the 
  #              state vector are. 
  #
  if (is.null(Groups)) Groups <- 1:NEnsTS
  alpha <- if (UseAlpha) theta[1] else 1
  IdxShift <- 1-as.numeric(UseAlpha) # Shifts indices in theta
  sigsq.0 <- exp(theta[2-IdxShift])
  tausq.0 <- exp(theta[3-IdxShift])
  tausq.d <- exp(theta[4-IdxShift])
  sigsq.1 <- exp(theta[5-IdxShift])
  tausq.1 <- exp(theta[6-IdxShift])
  phi.0 <- phi.1 <- 0
  if (UsePhi) {
    phi.0 <- 1 / (1+exp(-theta[7-IdxShift])) # Inverse logit transform
    phi.1 <- 1 / (1+exp(-theta[8-IdxShift]))
  }
  NGroups <- max(Groups)
  nS <- 3*(NGroups+2) # Length of state vector
  FF <- matrix(0, nrow=NEnsTS+1, ncol=nS)
  FF[1,1] <- 1; FF[-1,1] <- alpha; FF[-1,4] <- 1 
  #
  # Next lines use matrix subsetting to pick out the correct columns
  # for each RCM and GCM combination in Groups.
  #
  FF[cbind(2:nrow(FF), (3*Groups)+4)] <- 1
  if (constrain) {
    FF[1+which(Groups==max(Groups)), seq(7, nS-5, 3)] <- -1 # Sum-to-zero constraint on final simulator
  }
  GG <- diag(c(phi.0,1,1, rep(c(phi.1, 1, 1), NGroups+1)))
  GG[row(GG)%%3 == 1 & col(GG)-row(GG)==1] <- 1
  GG[4,1] <- alpha*(phi.1-phi.0) # mu.death doesn't follow the overall pattern. 
  TVEls <- (row(GG)%%3 == 1 & col(GG)-row(GG)==2) # Time-varying elements
  GG[TVEls] <- -kappa # -kappa marks time-varying component
  JGG <- diag(rep(0,3*(NGroups+2)))
  JGG[TVEls] <- 1 # Time-varying elements of GG are from column 1 of X
  V <- c(sigsq.0, rep(sigsq.1, NEnsTS))
  if (!is.null(RunsPerTS)) { # Scale variances for averaged series
    if (!isTRUE(length(RunsPerTS)==NEnsTS)) {
      stop("RunsPerTS should be a vector of length NEnsTS")
    }
    V[-1] <- V[-1] / RunsPerTS
  }
  V <- diag(V)
  W <- diag(rep(0,nS))
  if (constrain) { # W is block diagonal if constraints are used
    diag(W)[c(2,5)] <- c(tausq.0, tausq.d)
    tau.els <- seq(8, nS-4, 3)
    W[tau.els, tau.els] <- -tausq.1 / NGroups
    diag(W)[tau.els] <- diag(W)[tau.els] + tausq.1
    nS <- nS-3 # Remove redundant elements of state vector 
    FF <- FF[,1:nS]
    GG <- GG[1:nS, 1:nS]
    JGG <- JGG[1:nS, 1:nS]
    W <- W[1:nS, 1:nS]
  } else { # Unconstrained version: W is diagonal
    diag(W)[c(2,5,5+(3*(1:NEnsTS)))] <- c(tausq.0, tausq.d, rep(tausq.1, NEnsTS))
  }
  init <- make.SSinits(GG,W,kappa=kappa)
  if (!is.null(m0)) {
    if (!isTRUE(length(m0)==nrow(GG))) {
      stop(paste("m0 should be a vector containing",nrow(GG),"values"))
    }
    init$m0 <- m0 - 
      (Shift*c(1, 1-phi.0, 0, 1-alpha, 1-((1-phi.1)*alpha)-phi.0, 
               rep(0, nS-5)))
  }
  if (constrain) { # For constrained model, set C0 to reflect sum-to-zero constraints as well
    init$C0[7:nS, 7:nS] <- -kappa / NGroups
    diag(init$C0)[7:nS] <- diag(init$C0)[7:nS] + kappa
  }
  if (!is.null(C0)) { # Overwrite if the user supplied a value
    if (!isTRUE(all(dim(C0)==nrow(GG)))) {
      stop(paste("C0 should be a ",nrow(GG),"x",nrow(GG)," matrix",sep=""))
    }
    init$C0 <- C0
  }
  nS <- ncol(FF) # dimension of state vector after eliminating unwanted rows
  Consensus <- rep(0, nS); Consensus[c(1,4)] <- c(alpha,1) # Coefficients to extract ensemble consensus
  z <- list(m0=init$m0,C0=init$C0, FF=FF, V=V, GG=GG, W=W, 
       JGG=JGG,X=as.matrix(Xt), Consensus=Consensus)
  attr(z, "Shift") <- as.numeric(
    Shift*c(1, 1-phi.0, 0, 1-alpha, 1-((1-phi.1)*alpha)-phi.0, 
            rep(0, nS-5)))
  z
  }
######################################################################
EnsEBMtrendSmooth <- 
  function(Y, Xt, Groups=NULL, m0=NULL, C0=NULL,  kappa=1e6,
           compact=!is.null(Groups), prior.pars=NULL, theta=NULL, 
           UseAlpha=TRUE, UsePhi=TRUE, constrain=TRUE, messages=TRUE,
           Use.dlm=FALSE, debug=FALSE, ...) {
  #
  #  Fits and applies a model of the form defined by EnsEBMtrend.modeldef.
  #  Arguments: 
  #  
  #  Y        Matrix containing an observed time series in its first column 
  #           and exchangeable series from an ensemble in the remaining 
  #           columns
  #  Xt       Vector containing the corresponding forcings for each row of Y
  #  Groups   Vector indicating group membership for each ensemble member. 
  #           This should be used if an ensemble contains multiple 
  #           runs from one or more simulators: in this case the 
  #           elements of Groups should be integers between 1 and G
  #           where G is the total number of simulators contributing
  #           to the ensemble: Groups[i] is the index of the simulator
  #           that was used to produce the ith ensemble member. If NULL, 
  #           each member is assumed to come from a different simulator. 
  #  m0       Optional vector of initial values for the state vector. 
  #           If NULL (the default) this is determined automatically
  #           from the model structure. 
  #  kappa    Variance used to initialise diffuse elements of the state vector. 
  #  compact  Logical scalar indicating whether to compact the
  #           ensemble using CompactEns() before fitting. For 
  #           ensembles with many members per simulator, this
  #           may speed up the fitting without sacrificing 
  #           information because the simulator-specific means
  #           are sufficient statistics. Defaults to TRUE if
  #           Groups is non-NULL and FALSE otherwise (because 
  #           CompactEns() just returns the original ensemble
  #           in that case). 
  #  UseAlpha Logical scalar indicating whether or not to include
  #           a scaling factor to represent the expected value of
  #           the ensemble consensus trend slope beta[d] given the 
  #           observed one beta[0]. If FALSE, the expected ensemble
  #           consensus trend is equal to the observed one. 
  #  UsePhi   Controls whether or not to include thermal inertia 
  #           parameters in the model structure. 
  #  prior.pars Optional 2-column matrix containing the means and standard 
  #           deviations of independent Gaussian priors for the transformed 
  #           model parameters (log variances, and logit thermal inertias). 
  #	  theta	  Optional initial value for parameter vector. If NULL,
  #           routine will auto-initialise the optimisation
  #  constrain Logical scalar controlling whether to impose sum-to-zero 
  #           constraints where necessary to make the model identifiable. 
  #  messages Logical scalar controlling whether or not to print progress 
  #           to screen. 
  #   debug   Logical: argument passed to dlm routines to control
  #           whether to use faster C code (the default), or slower 
  #           R code which allows problem tracing. 
  #   ...     Other arguments to dlm.SafeMLE.
  #
  #  The function returns a list with three components: the first is the
  #  object containing results of the maximum likelihood parameter estimation, 
  #  the second is the fitted model itself, and the third is the result of 
  #  applying the Kalman Smoother to the input series using this fitted model. 
  #
  #   Check for partially missing rows of the ensemble data
  #
  if (!NA.AllOr0(Y[,-1])) {
    warning(paste("Some but not all ensemble members contain missing",
                  "value at some time points.\nThis may cause problems when fitting",
                   "because the dlm library routines do not cope well\nwith missing",
                   "values in the observation vector unless the corresponding",
                   "measurement errors\nare uncorrelated."))
  }
  #
  #  Now compact the ensemble if requested (NB Y and Groups here
  #  are local copies in the active environment, hence re-using 
  #  them won't affect anything else).
  #
  if (!is.null(Groups)) {
    if (compact) {
      Y <- CompactEns(Y, Groups)
      Groups <- NULL # No need for groups now
      message("######\n",
              "######  NOTE: compacting Y and Groups. Do not use these objects",
              "\n######  when working with the results of this call: use the ",
              "\n######  versions in the `Data` component of the result.",
              "\n######")
    }
  } else {
    if (compact) warning("'compact' argument does nothing when Groups is NULL")
  }
  RunsPerTS <- attr(Y, "RunsPerTS") # Could be NULL
  NEnsTS <- ncol(Y)-1
  #
  #  Start by finding initial values for a numerical maximisation of the 
  #  likelihood. For the real-world series these are obtained via a call to
  #  EBMtrendSmooth.For parameters that appear in the specification of the
  #  ensemble consensus discrepancy, from a fit of EBMtrendSmooth to the 
  #  difference between the ensemble mean and an estimate of the real-world
  #  trend. For the additional member-specific drift terms, just use 
  #  SLLT.IniPar to get a plausible initial value for sigsq.1.
  #
  par.names <- c("alpha", "log(sigsq[0])", "log(tausq[0])", 
                 "log(tausq[w])", "log(sigsq[1])", "log(tausq[1])")
  if (UsePhi) par.names <- c(par.names, "logit(phi[0])", "logit(phi[1])")
  if (!UseAlpha) par.names <- par.names[-1]
  if (!is.null(theta)) {
    theta.init <- as.numeric(theta)
  } else {
    if (as.numeric(messages)>0) {
      cat("Finding starting values for maximum likelihood estimation ...\n")
    } 
    theta.init <- rep(NA, ifelse(UsePhi,8,6))
    theta.init[1] <- 1 # Initialise alpha to 1 (will remove this later if not needed)
    cur.priors <- prior.pars[if (UsePhi) c(2:3,7) else 2:3, ]
    if (is.null(m0)) m0.0 <- NULL else m0.0 <- m0[1:3]
    if (is.null(C0)) C0.0 <- NULL else C0.0 <- C0[1:3, 1:3]
    Model0 <- 
      EBMtrendSmooth(Y[,1], Xt, m0=m0.0, C0=C0.0, 
                     kappa=kappa, UsePhi=UsePhi, prior.pars=cur.priors, 
                     messages=FALSE,  Use.dlm=Use.dlm, debug=debug, ...)
    theta.init[2:3] <- Model0$Theta$par[1:2] # sigsq.0 and tausq.0
    if (UsePhi) theta.init[7] <- Model0$Theta$par[3] # phi.0
    EnsembleMean <- EnsMean(Y, Groups)
    #
    #  Here's the smooth of difference between ensemble mean and 
    #  real-world trend
    #
    Death.hat <- EnsembleMean-Model0$Smooth$s[-1,1] # Remove first element (time 0)
    IdxShift <- 1-as.numeric(UseAlpha) # Shifts indices in theta
    cur.priors <- prior.pars[(if (UsePhi) c(2,4,8) else c(2,4))-IdxShift,]
    cur.priors[1,1] <- cur.priors[1,1] - log(NEnsTS) # Suppress interannual variability in trend
    if (is.null(m0)) m0.0 <- NULL else m0.0 <- m0[4:6]
    if (is.null(C0)) C0.0 <- NULL else C0.0 <- C0[4:6, 4:6]
    Model0 <- 
      EBMtrendSmooth(Death.hat, Xt, m0=m0.0, C0=C0.0, kappa=kappa, 
                     UsePhi=UsePhi, prior.pars=cur.priors, 
                     messages=FALSE, Use.dlm=Use.dlm, debug=debug, ...)
    theta.init[4] <- Model0$Theta$par[2] # Slope innovation parameter
    if (UsePhi) theta.init[8] <- Model0$Theta$par[3]    
    theta.init[5:6] <- SLLT.IniPar(Y[,-1]-EnsembleMean)
    if (!UseAlpha) theta.init <- theta.init[-1]
  }
  names(theta.init) <- par.names
  #
  #  Now the estimation, printing information if required (removing
  #  alpha first if necessary)
  #
  if (as.numeric(messages)>0) {
    cat("Estimating model parameters - please be patient ...\n")
  }
  Ybar <- mean(Y[,1], na.rm=TRUE)
  thetahat <- 
    dlm.SafeMLE(theta.init, Y-Ybar, EnsEBMtrend.modeldef, Xt=Xt, m0=m0, 
                C0=C0, kappa=kappa, Shift=Ybar, prior.pars=prior.pars,
                NEnsTS=NEnsTS, RunsPerTS=RunsPerTS, Groups=Groups,
                constrain=constrain, UseAlpha=UseAlpha, UsePhi=UsePhi,
                hessian=TRUE, par.names=par.names, Use.dlm=Use.dlm, 
                messages=messages, debug=debug, ...)
  Model <- # Original scale
    EnsEBMtrend.modeldef(thetahat$par, Xt=Xt, m0=m0, C0=C0, 
                         kappa=kappa, NEnsTS=NEnsTS, RunsPerTS=RunsPerTS, 
                         Groups=Groups, UseAlpha=UseAlpha, UsePhi=UsePhi, 
                         constrain=constrain)  
  ShiftMod <- # Shifted
    EnsEBMtrend.modeldef(thetahat$par, Xt=Xt, m0=m0, C0=C0, 
                         kappa=kappa, Shift=Ybar, NEnsTS=NEnsTS, 
                         RunsPerTS=RunsPerTS, Groups=Groups, 
                         UseAlpha=UseAlpha, UsePhi=UsePhi, 
                         constrain=constrain)  
  if (as.numeric(messages)>0) {
    cat("\nSUMMARY OF STATE SPACE MODEL:\n")
    print(summary.dlmMLE(thetahat))
  }
  #
  #  Kalman Smoother, and return result
  #
  Smooth <- tryCatch(dlmSmooth(Y-Ybar, ShiftMod, debug=debug), 
                     error=function(e) { 
                       z <- rep(NA, length=Y)
                       attr(z, "error") <- e; z
                     }
  )
  if (!is.null(attr(Smooth, "error"))) {
    warning(paste("dlmSmooth failed - error message was\n ", 
                  attr(Smooth, "error")$message), immediate.=TRUE)
  } else {
    Smooth$s <-  # Add the shift back
      t(t(Smooth$s) + attr(ShiftMod, "Shift"))
  }
  list(Data=list(Y=Y, Groups=Groups, RunsPerTS=RunsPerTS),
       Theta=thetahat, Model=Model, Smooth=Smooth)
}
######################################################################
EnsEBM2waytrend.modeldef <- 
  function(theta, Xt, m0=NULL, C0=NULL, kappa=1e6, Shift=0, 
           Groups, RunsPerTS=NULL, interactions="none", UseAlpha=TRUE, 
           UsePhi=TRUE, constrain=TRUE) {
    #
    #  Sets up the structure of a state space model for an observed
    #  time series and a regional ensemble (i.e. containing runs
    #  obtained from combinations of RCMs driven by GCMs), with 
    #  trend formulation inspired by the form of a simple energy 
    #  balance model (EBM). Arguments: 
    #
    #  theta  Parameter vector, containing eight elements if UseAlpha 
    #  and UsePhi are both TRUE (see below) in which case it is defined 
    #  as follows:
    #
    #         theta[1]    This is alpha, the scaling on the relationship 
    #                     the ensemble consensus and observed drift processes. 
    #         theta[2]    This is log(sigsq[0]), where sigsq[0] is the
    #                     "measurement error" variance in the observed 
    #                     series. 
    #         theta[3]    This is log(tausq[0]), where tausq[0] is the slope 
    #                     innovation variance in the observed series. 
    #         theta[4]    This is log(tausq[death]), where tausq[death] 
    #                     is the slope innovation variance for the shared 
    #                     discrepancy between the ensemble drift and the 
    #                     real drift. 
    #         theta[5]    This is log(sigsq[1]), where sigsq[1] is the common 
    #                     "measurement error" variance in the ensemble members.
    #         theta[6]    This is log(tausq[1]), where tausq[1] is the 
    #                     common drift innovation variance for the ensemble 
    #                     members about the shared discrepancy.
    #         theta[7]    This is logit(phi[0]), where phi[0] is the 
    #                     temporal dependence parameter in the mean 
    #                     component of the trend in the observations.
    #         theta[8]    This is logit(phi[1]), where phi[1] is the common 
    #                     temporal dependence parameter in the mean component 
    #                     of the trends in both the observations and the 
    #                     ensemble members. 
    #
    #  If UseAlpha is FALSE then the value of alpha is fixed at 1 and is 
    #  omitted from theta so that theta[1] is log(sigsq[0]), theta[2] is
    #  log(tausq[0]) etc. If UsePhi is FALSE, then both of the final 
    #  elements of theta are omitted and the values of both phi[0] and
    #  phi[1] are set to 0. 
    #
    #  Xt     Vector containing a time series of effective forcings. 
    #  Groups Two-column matrix specifying the group membership structure
    #         of the ensemble. E.g. if the ensemble contains runs from some
    #         combination of RCMs 1, 2, 3, ..., R driven by GCMs 
    #         1, 2, 3, ..., G then Groups[i,1] and Groups[i,2] should
    #         respectively contain the numbers of the RCM and GCM used 
    #         for the ith ensemble member. The values of R and G will 
    #         be taken from the largest numbers found in columns 1 and
    #         2 respectively.
    #  m0     Optional vector of initial values for the state vector. 
    #         If NULL (the default) this is determined automatically
    #         from the model structure. 
    #  C0     Optional covariance matrix for initialising the state
    #         vector. If NULL, like m0 this is determined from the
    #         model structure and from kappa. 
    #  kappa  "Large" value with which to initialise variances of diffuse
    #         components of state vector at time 0. Taking it too large 
    #         can lead to numerical instability, however. It should be 
    #         several orders of magnitude larger than any of the model 
    #         variances. It is ignored if C0 is non-NULL.
    #  interactions Logical scalar, indicating which RCM:GCM "interaction"
    #         terms to include in the state vector. Options are "none" to 
    #         exclude them completely, "all" to include all of them
    #         and "available" to include just those for which the 
    #         corresponding RCM:GCM combinations appear in Groups. Default
    #         is "none", because the interaction terms substantially 
    #         increase the dimension and are therefore likely to slow 
    #         things down considerably. 
    #   UseAlpha   Logical scalar indicating whether or not to include
    #         a scaling factor to represent the expected value of
    #         the ensemble consensus trend slope beta[d] given the 
    #         observed one beta[0]. If FALSE, the expected ensemble
    #         consensus trend is equal to the observed one. 
    #   UsePhi  Controls whether or not to include thermal inertia 
    #         parameters in the model structure. 
    #   constrain  Logical scalar indicating whether or not to constrain
    #         the discrepancy trends and slopes for the individual
    #         ensemble members to sum to zero. Without this constraint,
    #         the model is formally not identifiable - although linear
    #         combinations of elements of the state vector are. 
    #
    #  The state vector for the underlying model consists of triplets
    #  {mu, beta, gamma} where mu is the value of the EBM-inspired trend, 
    #  beta is the approximation drift and gamma is the coefficient of
    #  effective radiative forcing in the EBM approximation. The 
    #  triplets are ordered as: truth, common discrepancy, RCM main 
    #  effects, GCM main effects, interactions (ordered as 11, 12, ..., 1G, 
    #  21, 22, ..., 2G, ..., R1, R2, ..., RG in what is hopefully an 
    #  obvious notation)
    #
    if (!isTRUE(interactions %in% c("none", "all", "available"))) {
      stop("'interactions' must be one of 'none', 'all' or 'available'")
    }
    alpha <- if (UseAlpha) theta[1] else 1
    IdxShift <- 1-as.numeric(UseAlpha) # Shifts indices in theta
    sigsq.0 <- exp(theta[2-IdxShift])
    tausq.0 <- exp(theta[3-IdxShift])
    tausq.d <- exp(theta[4-IdxShift])
    sigsq.1 <- exp(theta[5-IdxShift])
    tausq.1 <- exp(theta[6-IdxShift])
    phi.0 <- phi.1 <- 0
    if (UsePhi) {
      phi.0 <- 1 / (1+exp(-theta[7-IdxShift])) # Inverse logit transform
      phi.1 <- 1 / (1+exp(-theta[8-IdxShift]))
    }
    R <- max(Groups[,1]) # Number of RCMs
    G <- max(Groups[,2]) # Number of GCMs
    NEnsTS <- nrow(Groups) # Total number of series in ensemble 
    nS <- 3*(2+R+G+(R*G)) # Length of state vector including all interactions
    R.blocks <- 3:(R+2) # Blocks of state vector corresponding to main RCM effects ...
    G.blocks <- (R+3):(R+G+2) # ... and to main GCM effects ...
    RG.blocks <- (R+G+3):(R+G+(R*G)+2) # ... and to interactions
    FF <- matrix(0, nrow=NEnsTS+1, ncol=nS)
    FF[1,1] <- 1; FF[-1, 1] <- alpha; FF[-1,4] <- 1
    #
    # Next lines use matrix subsetting to pick out the correct columns
    # for each RCM and GCM combination in Groups.
    #
    FF[cbind(2:nrow(FF), 3*(R.blocks[Groups[,1]]-1)+1)] <- 1
    FF[cbind(2:nrow(FF), 3*(G.blocks[Groups[,2]]-1)+1)] <- 1
    FF[cbind(2:nrow(FF), 3*(RG.blocks[G*(Groups[,1]-1)+Groups[,2]]-1)+1)] <- 1
    if (constrain) {
      #
      # Insert -1s in the rows corresponding to runs using the final RCM and GCM.
      # This is most easily done by noting that the relevant coding is the 
      # same as that used by a linear model with factor covariates and sum-to-
      # zero constraints on the coefficients. 
      # 
      RCMs <- as.factor(Groups[,1]); contrasts(RCMs) <- "contr.sum"
      GCMs <- as.factor(Groups[,2]); contrasts(GCMs) <- "contr.sum"
      FF[-1,3*(c(R.blocks[-R], G.blocks[-G], 
                 RG.blocks[rep(1:(R-1), each=G-1)+(1:(G-1))-1])-1)+1] <- 
        model.matrix(~RCMs*GCMs)[,-1]
    }
    GG <- diag(c(phi.0,1,1, rep(c(phi.1, 1, 1), 1 + R+G+(R*G))))
    GG[row(GG)%%3 == 1 & col(GG)-row(GG)==1] <- 1 # Contributions of {beta[t-1]} to {mu[t]}
    GG[4,1] <- alpha*(phi.1-phi.0) # mu.death doesn't follow the overall pattern. 
    TVEls <- (row(GG)%%3 == 1 & col(GG)-row(GG)==2) # Time-varying elements
    GG[TVEls] <- -kappa # -kappa marks time-varying component
    JGG <- diag(rep(0,nrow(GG)))
    JGG[TVEls] <- 1 # Time-varying elements of GG are from column 1 of X
    V <- c(sigsq.0, rep(sigsq.1, NEnsTS))
    if (!is.null(RunsPerTS)) { # Scale variances for averaged series
      if (!isTRUE(length(RunsPerTS)==NEnsTS)) {
        stop("RunsPerTS should be a vector of length NEnsTS")
      }
      V[-1] <- V[-1] / RunsPerTS
    }
    V <- diag(V)
    W <- diag(rep(0, nS))
    if (constrain) {
      #
      #   If constraints are applied then structure of W needs modifying: 
      #   deal with blocks corresponding to RCMs, GCMs and interactions
      #   separately. For simplicity, the redundant elements are also
      #   set here (with incorrect values); these are subsequently deleted.
      #
      diag(W)[c(2,5)] <- c(tausq.0, tausq.d)
      tau.els <- 3*R.blocks-1 # Second elements of triplets for each RCM block ...
      W[tau.els, tau.els] <- -tausq.1/R
      diag(W)[tau.els] <- tausq.1*(R-1)/R
      tau.els <- 3*G.blocks-1 # ... each GCM block ...
      W[tau.els, tau.els] <- -tausq.1/G
      diag(W)[tau.els] <- tausq.1*(G-1)/G
      #
      #   ... and each interaction block. This is a bit trickier, 
      #   so build up the required structure in a separate matrix
      #   and then insert into the appropriate positions of W. 
      #
      tau.els <- 3*RG.blocks-1 # ... and each interaction block
      M <- matrix(1, nrow=R*G, ncol=R*G)
      diag(M) <- rep((R-1)*(G-1), R*G)
      Wanted.Els <- # Indices of diagonal blocks of M
        cbind(rep(1:(R*G), each=G), rep(G*(0:(R-1)), each=G*G)+rep(1:G,G))
      Wanted.Els <- Wanted.Els[Wanted.Els[,1] != Wanted.Els[,2],]
      M[Wanted.Els] <- -(R-1)
      Wanted.Els <- # Diagonal elements of off-diagonal blocks
        cbind(rep(1:(R*G), each=R), G*(0:(R-1))+rep(1:G, each=R))
      Wanted.Els <- Wanted.Els[Wanted.Els[,1] != Wanted.Els[,2],]
      M[Wanted.Els] <- -(G-1)
      W[tau.els, tau.els] <- M * tausq.1/(R*G)
    } else {
      diag(W)[c(2,5,5+(3*(1:(R+G+(R*G)))))] <- 
        c(tausq.0, tausq.d, rep(tausq.1, (nS/3)-2))
    }
    #
    #   That completes the structure of a model including all interaction
    #   terms: next, remove any elements that are not needed.
    #
    Keep <- 1:ncol(FF)
    BlocksToKeep <- 1:(2+R+G+(R*G)) # Vector containing all blocks
    if (!isTRUE(interactions=="all")) {
      BlocksToKeep <- 1:(2+R+G)
      if (interactions=="available") {
        BlocksToKeep <- c(BlocksToKeep, 
                          RG.blocks[unique( (G*(Groups[,1]-1))+Groups[,2])])
      }
    }
    if (constrain) {
      BlocksToDrop <- # Anything involving RCM R or GCM G
        c(2+R, 2+R+G, 2+R+G+(G*(1:(R-1))), 2+(R*G))
      BlocksToKeep <- BlocksToKeep[!(BlocksToKeep %in% BlocksToDrop)]
    }
    Keep <- which( ((0:(ncol(FF)-1)) %/% 3) %in% (BlocksToKeep-1) )
    FF <- FF[,Keep]; W <- W[Keep, Keep]
    GG <- GG[Keep, Keep]; JGG <- JGG[Keep, Keep]
    init <- make.SSinits(GG, W, kappa=kappa)
    if (!is.null(m0)) {
      if (!isTRUE(length(m0)==nrow(GG))) {
        stop(paste("m0 should be a vector containing",nrow(GG),"values"))
      }
      init$m0 <- m0 - 
        (Shift*c(1, 1-phi.0, 0, 1-alpha, 1-((1-phi.1)*alpha)-phi.0, 
                 rep(0, nrow(GG)-5)))
    }
    if (constrain) { # For constrained model, set C0 to reflect sum-to-zero 
                     # constraints as well. This is a bit tricky: copy scaled
                     # versions of the relevant submatrix of W into each of
                     # the three blocks corresponding to mu, beta and gamma. 
                     # The scaling below preserves the pattern in W, and
                     # ensures that the median diagonal element of the
                     # relevant block of C0 is equal to kappa. 
      MinVal <- min(abs(diag(W)[abs(diag(W))>0]))
      Wanted.Els <- which(rowSums(abs(W))-abs(diag(W)) > MinVal / 1e6)
      MedVal <- median(diag(W)[Wanted.Els])
      init$C0[Wanted.Els, Wanted.Els] <- 
        kappa * W[Wanted.Els, Wanted.Els] / MedVal
      init$C0[Wanted.Els-1, Wanted.Els-1] <- 
        kappa * W[Wanted.Els, Wanted.Els] / MedVal
      init$C0[Wanted.Els+1, Wanted.Els+1] <- 
        kappa * W[Wanted.Els, Wanted.Els] / MedVal
    }
    if (!is.null(C0)) { # Overwrite if the user supplied a value
      if (!isTRUE(all(dim(C0)==nrow(GG)))) {
        stop(paste("C0 should be a ",nrow(GG),"x",nrow(GG)," matrix",sep=""))
      }
      init$C0 <- C0
    }
    nS <- ncol(FF) # dimension of state vector after eliminating unwanted rows
    Consensus <- rep(0, nS); Consensus[c(1,4)] <- c(alpha,1) # Coefficients to extract ensemble consensus
    z <- list(FF=FF, GG=GG, V=V, W=W, JGG=JGG, X=as.matrix(Xt), 
              m0=init$m0, C0=init$C0, Consensus=Consensus)
    attr(z, "Shift") <- as.numeric(
      Shift*c(1, 1-phi.0, 0, 1-alpha, 1-((1-phi.1)*alpha)-phi.0, 
              rep(0, nS-5)))
    z
  }
######################################################################
EnsEBM2waytrendSmooth <- 
  function(Y, Xt, m0=NULL, C0=NULL, kappa=1e6, Groups, compact=TRUE, 
           prior.pars=NULL, theta=NULL, interactions="none", UseAlpha=TRUE,
           UsePhi=TRUE, constrain=TRUE, messages=TRUE, Use.dlm=FALSE, 
           debug=FALSE, ...) {
  #
  #  Fits and applies a model of the form defined by EnsEBM2waytrend.modeldef.
  #  Arguments: 
  #
  #   Y       Matrix containing an observed time series in its first 
  #           column and series from a regional ensemble (i.e. an 
  #           ensemble in which each series is produced by one of R RCMS 
  #           with boundary conditions provided by one of G GCMs) in the 
  #           remaining columns.
  #   Xt      Vector containing an external forcing time series: each 
  #           element is the forcing for thhe corresponding row of Y. 
  #   m0      Optional vector of initial values for the state vector. 
  #           If NULL (the default) this is determined automatically
  #           from the model structure. 
  #   kappa   Variance used to initialise diffuse elements of the state 
  #           vector. 
  #   prior.pars  An optional 7*2 matrix containing means and standard 
  #           deviations of independent Gaussian priors for the 
  #           transformed parameters (log variances / logit inertias) 
  #           in the model. If this is provided, maximum a posteriori 
  #           (or penalised maximum likelihood) estimation is used; 
  #           otherwise just normal maximum likelihood. 
  #	  theta	  Optional initial value for parameter vector. If NULL,
  #           the routine will auto-initialise the optimisation
  #   Groups  Two-column matrix such that Groups[i,1] is the number of
  #           the RCM used to produce the ith ensemble member (i.e. 
  #           column i+1 of Y) and Groups[i,2] is the number of the
  #           GCM. 
  #   compact     Logical scalar indicating whether to compact the
  #               ensemble using CompactEns() before fitting. For 
  #               ensembles with many members per simulator, this
  #               may speed up the fitting without sacrificing 
  #               information because the simulator-specific means
  #               are sufficient statistics. Defaults to TRUE. 
  #  interactions   Logical scalar, indicating which RCM:GCM "interaction"
  #           terms to include in the state vector. Options are "none" to 
  #           exclude them completely, "all" to include all of them
  #           and "available" to include just those for which the 
  #           corresponding RCM:GCM combinations appear in Groups. Default
  #           is "none", because the interaction terms substantially 
  #           increase the dimension and are therefore likely to slow 
  #           things down considerably. 
  #   UseAlpha   Logical scalar indicating whether or not to include
  #           a scaling factor to represent the expected value of
  #           the ensemble consensus trend slope beta[d] given the 
  #           observed one beta[0]. If FALSE, the expected ensemble
  #           consensus trend is equal to the observed one. 
  #   UsePhi  Controls whether or not to include thermal inertia 
  #           parameters in the model structure. 
  #   constrain  Logical scalar indicating whether or not to constrain
  #           the discrepancy trends and slopes for the individual
  #           ensemble members to sum to zero. Without this constraint,
  #           the model is formally not identifiable - although linear
  #           combinations of elements of the state vector are. 
  #   messages  Logical scalar. If TRUE (the default), progress is
  #           printed to screen. 
  #   Use.dlm  Logical scalar. If TRUE then the dlmMLE() function will be
  #           used to maximise the log-likelihood; if FALSE (the default)
  #           then nlm will be used instead. *** NOTE*** only FALSE
  #           is currently supported for this model, and it does not 
  #           currently call dlm.SafeMLE() due to a memory leak in the C 
  #           code for dlmLL. 
  #   debug   As in dlmMLE. 
  #   ...     Additional arguments to dlm.SafeMLE
  #
  #  The function returns a list with three components: the first is the 
  #  object containing results of the maximum likelihood parameter 
  #  estimation, the second is the fitted model itself, and the third is 
  #  the result of applying the Kalman Smoother to the input series using 
  #  this fitted model. 
  #
  #
  #   Check for partially missing rows of the ensemble data
  #
  if (!NA.AllOr0(Y[,-1])) {
    warning(paste("Some but not all ensemble members contain missing",
                  "value at some time points.\nThis may cause problems when fitting",
                  "because the dlm library routines do not cope well\nwith missing",
                  "values in the observation vector unless the corresponding",
                  "measurement errors\nare uncorrelated."))
  }
  #
  #  Now compact the ensemble if requested. NB Y and Groups 
  #  here are local copies in the active environment, hence re-using 
  #  them won't affect anything else.
  #
  RunsPerTS <- NULL
  if (compact) {
    GrpStrs <- apply(Groups, MARGIN=1, FUN=paste, collapse=".&.")
    Y <- CompactEns(Y, GrpStrs)
    Groups <- Groups[!duplicated(GrpStrs),]
    GrpStrs <- GrpStrs[!duplicated(GrpStrs)]
    NewOrder <- match(colnames(Y[,-1]), GrpStrs) # Ensure that
    GrpStrs <- GrpStrs[NewOrder]                 # everything is in
    Groups <- Groups[NewOrder,]                  # the right order
    message("######\n",
            "######  NOTE: compacting Y and Groups. Do not use these objects",
            "\n######  when working with the results of this call: use the ",
            "\n######  versions in the `Data` component of the result.",
            "\n######")
    RunsPerTS <- attr(Y, "RunsPerTS")   # Could be NULL
    if (!identical(colnames(Y)[-1],GrpStrs)) {
      stop(paste("Ensemble compaction has re-ordered combinations.",
                  "This is a\nprogramming error: contact package author.",
                 sep=""))
    }
  }
  #
  #  Next, find initial values for a numerical maximisation of the 
  #  likelihood. For the observed series these are obtained as in 
  #  EBMtrendSmooth. The variance of slope innovations for the shared
  #  drift is initialised from the corresponding variance from an EBM
  #  mimic fitted to the difference between the estimated mean for the 
  #  observations and the ensemble mean. For the additional 
  #  member-specific drift terms, use SLLT.IniPar to get some plausible 
  #  initial values. 
  #
  par.names <- c("alpha", "log(sigsq[0])", "log(tausq[0])", 
                 "log(tausq[w])", "log(sigsq[1])", "log(tausq[1])")
  if (UsePhi) par.names <- c(par.names, "logit(phi[0])", "logit(phi[1])")
  if (!UseAlpha) par.names <- par.names[-1]
  if (!is.null(theta)) {
    theta.init <- as.numeric(theta)
  } else {
    if (as.numeric(messages)>0) {
      cat("Finding starting values for maximum likelihood estimation ...\n")
    }
    theta.init <- rep(NA, ifelse(UsePhi,8,6))
    theta.init[1] <- 1 # Initialise alpha to 1 (will remove this later if not needed)
    cur.priors <- prior.pars[if (UsePhi) c(2:3,7) else 2:3, ]
    Model0 <- EBMtrendSmooth(Y[,1], Xt, m0=m0[1:3], kappa=kappa, UsePhi=UsePhi, 
                             prior.pars=cur.priors, messages=FALSE, 
                             Use.dlm=Use.dlm, debug=debug, ...)
    theta.init[2:3] <- Model0$Theta$par[1:2] # sigsq.0 and tausq.0
    if (UsePhi) theta.init[7] <- Model0$Theta$par[3] # phi.0
    #
    #  Here's the smooth of difference between ensemble mean and 
    #  real-world trend
    #
    EnsembleMean <- EnsMean(Y, Groups)
    Death.hat <- EnsembleMean-Model0$Smooth$s[-1,1] # Discard first element of Model0 (time 0)
    IdxShift <- 1-as.numeric(UseAlpha) # Shifts indices in theta
    cur.priors <- prior.pars[(if (UsePhi) c(2,4,8) else c(2,4))-IdxShift,]
    cur.priors[1,1] <- cur.priors[1,1] - log(ncol(Y)-1) # Suppress interannual variability in trend
    if (is.null(m0)) m0.0 <- NULL else m0.0 <- m0[4:6]
    if (is.null(C0)) C0.0 <- NULL else C0.0 <- C0[4:6, 4:6]
    Model0 <- 
      EBMtrendSmooth(Death.hat, Xt, m0=m0.0, C0=C0.0, kappa=kappa, 
                     UsePhi=UsePhi, prior.pars=cur.priors, messages=FALSE, 
                     Use.dlm=Use.dlm, debug=debug, ...)
    theta.init[4] <- Model0$Theta$par[2]
    if (UsePhi) theta.init[8] <- Model0$Theta$par[3]    
    theta.init[5:6] <- SLLT.IniPar(Y[,-1]-EnsembleMean)
    if (!UseAlpha) theta.init <- theta.init[-1]
  }
  names(theta.init) <- par.names
  #
  #  Now the estimation, printing information if required (removing
  #  alpha first if necessary)
  #
  if (as.numeric(messages)>0) {
    cat("Estimating model parameters - please be patient ...\n")
  }
  if (!UseAlpha) { 
    theta.init <- theta.init[-1] 
    par.names <- par.names[-1] 
  }
  Ybar <- mean(Y[,1], na.rm=TRUE)
  thetahat <- 
    dlm.SafeMLE(theta.init, Y-Ybar, EnsEBM2waytrend.modeldef, Xt=Xt, 
                m0=m0, C0=C0, kappa=kappa, Shift=Ybar, 
                prior.pars=prior.pars, Groups=Groups, 
                RunsPerTS=RunsPerTS, interactions=interactions, 
                UseAlpha=UseAlpha, UsePhi=UsePhi, constrain=constrain, 
                debug=debug, hessian=TRUE, par.names=par.names, 
                Use.dlm=Use.dlm, messages=messages, ...)  
  Model <- 
    EnsEBM2waytrend.modeldef(thetahat$par, Xt=Xt, m0=m0, C0=C0, 
                             kappa=kappa, Groups=Groups, 
                             RunsPerTS=RunsPerTS,
                             interactions=interactions, 
                             UseAlpha=UseAlpha, UsePhi=UsePhi, 
                             constrain=constrain)
  ShiftMod <- 
    EnsEBM2waytrend.modeldef(thetahat$par, Xt=Xt, m0=m0, C0=C0, 
                             kappa=kappa, Shift=Ybar, Groups=Groups, 
                             RunsPerTS=RunsPerTS, 
                             interactions=interactions, 
                             UseAlpha=UseAlpha, UsePhi=UsePhi, 
                             constrain=constrain)
  if (as.numeric(messages)>0) {
    cat("\nSUMMARY OF STATE SPACE MODEL:\n")
    print(summary.dlmMLE(thetahat))
  }
  #
  #  Kalman Smoother, and return result
  #
  Smooth <- tryCatch(dlmSmooth(Y-Ybar, ShiftMod, debug=debug), 
                     error=function(e) { 
                       z <- rep(NA, length=Y)
                       attr(z, "error") <- e; z
                     }
  )
  if (!is.null(attr(Smooth, "error"))) {
    warning(paste("dlmSmooth failed - error message was\n ", 
                  attr(Smooth, "error")$message), immediate.=TRUE)
  } else {
    Smooth$s <-  # Add the shift back
      t(t(Smooth$s) + attr(ShiftMod, "Shift"))
  }
  list(Data=list(Y=Y, Groups=Groups),
       Theta=thetahat, Model=Model, Smooth=Smooth)
}
