% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00roxygen.R, R/advector.R, R/osa.R, R/sim.R
\name{TMB-interface}
\alias{TMB-interface}
\alias{MakeADFun}
\alias{sdreport}
\alias{ADREPORT}
\alias{REPORT}
\alias{getAll}
\alias{OBS}
\alias{checkConsistency}
\title{Interface to TMB}
\usage{
MakeADFun(
  func,
  parameters,
  random = NULL,
  profile = NULL,
  integrate = NULL,
  intern = FALSE,
  map = list(),
  ADreport = FALSE,
  silent = FALSE,
  ridge.correct = FALSE,
  ...
)

sdreport(obj, ...)

ADREPORT(x)

REPORT(x)

getAll(..., warn = TRUE)

OBS(x)

checkConsistency(obj, fast = TRUE, ...)
}
\arguments{
\item{func}{Function taking a parameter list (or parameter vector) as input.}

\item{parameters}{Parameter list (or parameter vector) used by \code{func}.}

\item{random}{As \link[TMB]{MakeADFun}.}

\item{profile}{As \link[TMB]{MakeADFun}.}

\item{integrate}{As \link[TMB]{MakeADFun}.}

\item{intern}{As \link[TMB]{MakeADFun}.}

\item{map}{As \link[TMB]{MakeADFun}.}

\item{ADreport}{As \link[TMB]{MakeADFun}.}

\item{silent}{As \link[TMB]{MakeADFun}.}

\item{ridge.correct}{Experimental}

\item{...}{Passed to TMB}

\item{obj}{TMB model object (output from \link{MakeADFun})}

\item{x}{Observation object}

\item{warn}{Give a warning if overwriting an existing object?}

\item{fast}{Pass \code{observation.name} to \code{TMB} ?}
}
\value{
TMB model object.
}
\description{
Interface to TMB
}
\section{Functions}{
\itemize{
\item \code{MakeADFun()}: Interface to \link[TMB]{MakeADFun}.

\item \code{sdreport()}: Interface to \link[TMB]{sdreport}.

\item \code{ADREPORT()}: Can be used inside the objective function to report quantities for which uncertainties will be calculated by \link{sdreport}.

\item \code{REPORT()}: Can be used inside the objective function to report quantities via the model object using \code{obj$report()}.

\item \code{getAll()}: Can be used to assign all parameter or data objects from a list inside the objective function.

\item \code{OBS()}: Mark the observation to be used by either \code{oneStepPredict} or by \code{obj$simulate}.
If your objective function is using an observation \code{x}, you simply need
to run \code{x <- OBS(x)} \emph{inside the objective function}.
This will (1) allow \code{oneStepPredict} to change the class of \code{x} to
\code{"osa"} (\link{OSA-residuals}) or (2) allow \code{obj$simulate} to change the class of \code{x} to
\code{"simref"} (\link{Simulation}) on request.

\item \code{checkConsistency()}: Interface to \link[TMB]{checkConsistency}.

}}
\examples{
data(rivers)
f <- function(p) { -sum(dnorm(rivers, p$mu, p$sd, log=TRUE)) }
obj <- MakeADFun(f, list(mu=0, sd=1), silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)
## 'REML trick'
obj2 <- MakeADFun(f, list(mu=0, sd=1), random="mu", silent=TRUE)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
sdreport(obj2) ## Compare with sd(rivers)
## Single argument vector function with numeric 'parameters'
fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
obj <- MakeADFun(fr, numeric(2), silent=TRUE)
nlminb(c(-1.2, 1), obj$fn, obj$gr, obj$he)
}
