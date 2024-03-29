% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00roxygen.R, R/methods.R
\name{ADapply}
\alias{ADapply}
\alias{apply,advector-method}
\alias{sapply,ANY-method}
\title{AD apply functions}
\usage{
\S4method{apply}{advector}(X, MARGIN, FUN, ..., simplify = TRUE)

\S4method{sapply}{ANY}(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
}
\arguments{
\item{X}{As \link[base]{apply}}

\item{MARGIN}{As \link[base]{apply}}

\item{FUN}{As \link[base]{apply}}

\item{...}{As \link[base]{apply}}

\item{simplify}{As \link[base]{sapply}}

\item{USE.NAMES}{As \link[base]{sapply}}
}
\value{
Object of class \code{"advector"} with a dimension attribute.
}
\description{
These \strong{base} apply methods have been modified to keep the AD class attribute (which would otherwise be lost).
}
\section{Functions}{
\itemize{
\item \code{apply(advector)}: As \link[base]{apply}

\item \code{sapply(ANY)}: As \link[base]{sapply}

}}
\examples{
F <- MakeTape(function(x) apply(matrix(x,2,2), 2, sum), numeric(4))
F$jacobian(1:4)
}
