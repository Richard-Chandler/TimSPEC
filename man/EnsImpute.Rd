\name{EnsImpute}
\alias{EnsImpute}
\alias{NApatterns}
\alias{NA.AllOr0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Utility functions for working with ensembles of time series}
\description{Utility functions for working with ensembles of time series: imputation of missing values, and missing value patterns. 
}
\usage{
EnsImpute(Y, method="lm")
NApatterns(Y)
NA.AllOr0(Y)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A matrix or data frame, with each column corresponding to an ensemble member and each row corresponding to a time point.}
  \item{method}{Either \code{"lm"} (the default) or \code{"approx"}. See 'Details' section below.}
}

\details{\code{EnsImpute} imputes missing values in \code{Y} using an additive ANOVA model with effects corresponding to ensemble member and to the time variable (treated as a factor). If \code{method="lm"} (the default), the effects are estimated using least squares: this can take a few seconds for ensembles with many members and / or time points. Alternatively, if \code{method="approx"} then an approximate fit is obtained using the formulae for a full two-way layout with no missing values: denoting by \eqn{Y_{it}}{Y[it]} the value for ensemble member \eqn{i} at time \eqn{t}, the corresponding fitted value is estimated as 

\deqn{\overline{Y}_{\cdot\cdot} + 
\left(\overline{Y}_{i\cdot}-\overline{Y}_{\cdot\cdot}\right) + 
\left(\overline{Y}_{\cdot t}-\overline{Y}_{\cdot\cdot}\right) = \overline{Y}_{i\cdot} + 
\left(\overline{Y}_{\cdot t}-\overline{Y}_{\cdot\cdot}\right)}{%
Ybar[..] + (Ybar[i.]-Ybar[..]) + (Ybar[.t]-Ybar[..]) = Ybar[i.] + 
(Ybar[.j]-Ybar[..])
}

in standard ANOVA notation. The \eqn{\{\overline{Y}_{i\cdot}\}}{{Ybar[i.]}} are estimated as the sample means of the available values for ensemble member \eqn{i}; and the temporal residuals \eqn{\{\overline{Y}_{\cdot t}-\overline{Y}_{\cdot\cdot}\}}{{Ybar[.t]-Ybar[..]}} are estimated from the set of ensemble members with complete records. 

\code{NApatterns} identifies the patterns of missingness in rows of \code{Y}.

\code{NA.AllOr0} checks whether the rows of \code{Y} are either all missing or all non-missing. This is because some of the routines in the \code{dlm} package can fail if dependent components of the data vector are partially missing.}

\value{\code{EnsImpute} returns an object like \code{Y}, with the \code{NA}s replaced as described in any row containing at least one non-missing values. \code{NA}s are preserved in rows that are completely missing, because in this case there is no basis for estimating the required time effect in the additive ANOVA.

\code{NApatterns} returns a character vector with an element for each row of \code{Y}, encoding the positions of the \code{NA} values in that row. For example, if the values in columns 1 and 3 are \code{NA} then the corresponding element of the result will be \code{"1,3"}. 

\code{NA.AllOr0} returns \code{TRUE} if all rows of \code{Y} are either all missing or all non-missing, and \code{FALSE} if there are some partially missing rows. 
}
\author{
Richard E. Chandler <r.chandler@ucl.ac.uk>
}
\note{These functions are intended to be called primarily from other routines in the \pkg{TimSPEC} package, and hence are not exported. If they are required (e.g. for debugging), use \code{TimSPEC:::EnsImpute(Y)} and \code{TimSPEC:::NApatterns(Y)}.}
