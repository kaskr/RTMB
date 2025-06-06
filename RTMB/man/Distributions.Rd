% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00roxygen.R, R/distributions.R, R/methods.R
\name{Distributions}
\alias{Distributions}
\alias{dexp,ad,ad.,logical.-method}
\alias{dexp,num,num.,logical.-method}
\alias{dexp,osa,ANY,ANY-method}
\alias{dexp,simref,ANY,ANY-method}
\alias{dweibull,ad,ad,ad.,logical.-method}
\alias{dweibull,num,num,num.,logical.-method}
\alias{dweibull,osa,ANY,ANY,ANY-method}
\alias{dweibull,simref,ANY,ANY,ANY-method}
\alias{dbinom,ad,ad,ad,logical.-method}
\alias{dbinom,num,num,num,logical.-method}
\alias{dbinom,osa,ANY,ANY,ANY-method}
\alias{dbinom,simref,ANY,ANY,ANY-method}
\alias{dbeta,ad,ad,ad,missing,logical.-method}
\alias{dbeta,num,num,num,missing,logical.-method}
\alias{dbeta,osa,ANY,ANY,ANY,ANY-method}
\alias{dbeta,simref,ANY,ANY,ANY,ANY-method}
\alias{df,ad,ad,ad,missing,logical.-method}
\alias{df,num,num,num,missing,logical.-method}
\alias{df,osa,ANY,ANY,ANY,ANY-method}
\alias{df,simref,ANY,ANY,ANY,ANY-method}
\alias{dlogis,ad,ad.,ad.,logical.-method}
\alias{dlogis,num,num.,num.,logical.-method}
\alias{dlogis,osa,ANY,ANY,ANY-method}
\alias{dlogis,simref,ANY,ANY,ANY-method}
\alias{dt,ad,ad,missing,logical.-method}
\alias{dt,num,num,missing,logical.-method}
\alias{dt,osa,ANY,ANY,ANY-method}
\alias{dt,simref,ANY,ANY,ANY-method}
\alias{dnbinom,ad,ad,ad,missing,logical.-method}
\alias{dnbinom,num,num,num,missing,logical.-method}
\alias{dnbinom,osa,ANY,ANY,ANY,ANY-method}
\alias{dnbinom,simref,ANY,ANY,ANY,ANY-method}
\alias{dpois,ad,ad,logical.-method}
\alias{dpois,num,num,logical.-method}
\alias{dpois,osa,ANY,ANY-method}
\alias{dpois,simref,ANY,ANY-method}
\alias{dgamma,ad,ad,missing,ad.,logical.-method}
\alias{dgamma,num,num,missing,num.,logical.-method}
\alias{dgamma,osa,ANY,ANY,ANY,ANY-method}
\alias{dgamma,simref,ANY,ANY,ANY,ANY-method}
\alias{pnorm,ad,ad.,ad.,missing,missing-method}
\alias{pnorm,num,num.,num.,missing,missing-method}
\alias{pgamma,ad,ad,missing,ad.,missing,missing-method}
\alias{pgamma,num,num,missing,num.,missing,missing-method}
\alias{ppois,ad,ad,missing,missing-method}
\alias{ppois,num,num,missing,missing-method}
\alias{pexp,ad,ad.,missing,missing-method}
\alias{pexp,num,num.,missing,missing-method}
\alias{pweibull,ad,ad,ad.,missing,missing-method}
\alias{pweibull,num,num,num.,missing,missing-method}
\alias{pbeta,ad,ad,ad,missing,missing,missing-method}
\alias{pbeta,num,num,num,missing,missing,missing-method}
\alias{qnorm,ad,ad.,ad.,missing,missing-method}
\alias{qnorm,num,num.,num.,missing,missing-method}
\alias{qgamma,ad,ad,missing,ad.,missing,missing-method}
\alias{qgamma,num,num,missing,num.,missing,missing-method}
\alias{qexp,ad,ad.,missing,missing-method}
\alias{qexp,num,num.,missing,missing-method}
\alias{qweibull,ad,ad,ad.,missing,missing-method}
\alias{qweibull,num,num,num.,missing,missing-method}
\alias{qbeta,ad,ad,ad,missing,missing,missing-method}
\alias{qbeta,num,num,num,missing,missing,missing-method}
\alias{besselK,ad,ad,missing-method}
\alias{besselK,num,num,missing-method}
\alias{besselI,ad,ad,missing-method}
\alias{besselI,num,num,missing-method}
\alias{besselJ,ad,ad-method}
\alias{besselJ,num,num-method}
\alias{besselY,ad,ad-method}
\alias{besselY,num,num-method}
\alias{dbinom_robust}
\alias{dsn}
\alias{dSHASHo}
\alias{dtweedie}
\alias{dnbinom_robust}
\alias{dnbinom2}
\alias{dlgamma}
\alias{logspace_add}
\alias{logspace_sub}
\alias{dnorm,ad,ad.,ad.,logical.-method}
\alias{dnorm,num,num.,num.,logical.-method}
\alias{dnorm,osa,ANY,ANY,ANY-method}
\alias{dnorm,simref,ANY,ANY,ANY-method}
\alias{dlnorm,ANY,ANY,ANY,ANY-method}
\alias{dlnorm,osa,ANY,ANY,ANY-method}
\alias{dlnorm,num,num.,num.,logical.-method}
\alias{plogis,advector,missing,missing,missing,missing-method}
\alias{qlogis,advector,missing,missing,missing,missing-method}
\alias{dcompois}
\alias{dcompois2}
\alias{pbinom,ad,ad,ad,missing,missing-method}
\alias{pbinom,num,num,num,missing,missing-method}
\alias{dmultinom,ad,ad.,ad,logical.-method}
\alias{dmultinom,num,num.,num,logical.-method}
\alias{dmultinom,osa,ANY,ANY,ANY-method}
\alias{dmultinom,simref,ANY,ANY,ANY-method}
\alias{dmultinom,ANY,ANY,ANY,ANY-method}
\alias{dcauchy,ad,ad.,ad.,logical.-method}
\alias{dcauchy,num,num.,num.,logical.-method}
\alias{dcauchy,osa,ANY,ANY,ANY-method}
\alias{dcauchy,simref,ANY,ANY,ANY-method}
\title{Distributions and special functions for which AD is implemented}
\usage{
\S4method{dexp}{ad,ad.,logical.}(x, rate = 1, log = FALSE)

\S4method{dexp}{num,num.,logical.}(x, rate = 1, log = FALSE)

\S4method{dexp}{osa,ANY,ANY}(x, rate = 1, log = FALSE)

\S4method{dexp}{simref,ANY,ANY}(x, rate = 1, log = FALSE)

\S4method{dweibull}{ad,ad,ad.,logical.}(x, shape, scale = 1, log = FALSE)

\S4method{dweibull}{num,num,num.,logical.}(x, shape, scale = 1, log = FALSE)

\S4method{dweibull}{osa,ANY,ANY,ANY}(x, shape, scale = 1, log = FALSE)

\S4method{dweibull}{simref,ANY,ANY,ANY}(x, shape, scale = 1, log = FALSE)

\S4method{dbinom}{ad,ad,ad,logical.}(x, size, prob, log = FALSE)

\S4method{dbinom}{num,num,num,logical.}(x, size, prob, log = FALSE)

\S4method{dbinom}{osa,ANY,ANY,ANY}(x, size, prob, log = FALSE)

\S4method{dbinom}{simref,ANY,ANY,ANY}(x, size, prob, log = FALSE)

\S4method{dbeta}{ad,ad,ad,missing,logical.}(x, shape1, shape2, log)

\S4method{dbeta}{num,num,num,missing,logical.}(x, shape1, shape2, log)

\S4method{dbeta}{osa,ANY,ANY,ANY,ANY}(x, shape1, shape2, log)

\S4method{dbeta}{simref,ANY,ANY,ANY,ANY}(x, shape1, shape2, log)

\S4method{df}{ad,ad,ad,missing,logical.}(x, df1, df2, log)

\S4method{df}{num,num,num,missing,logical.}(x, df1, df2, log)

\S4method{df}{osa,ANY,ANY,ANY,ANY}(x, df1, df2, log)

\S4method{df}{simref,ANY,ANY,ANY,ANY}(x, df1, df2, log)

\S4method{dlogis}{ad,ad.,ad.,logical.}(x, location = 0, scale = 1, log = FALSE)

\S4method{dlogis}{num,num.,num.,logical.}(x, location = 0, scale = 1, log = FALSE)

\S4method{dlogis}{osa,ANY,ANY,ANY}(x, location = 0, scale = 1, log = FALSE)

\S4method{dlogis}{simref,ANY,ANY,ANY}(x, location = 0, scale = 1, log = FALSE)

\S4method{dt}{ad,ad,missing,logical.}(x, df, log)

\S4method{dt}{num,num,missing,logical.}(x, df, log)

\S4method{dt}{osa,ANY,ANY,ANY}(x, df, log)

\S4method{dt}{simref,ANY,ANY,ANY}(x, df, log)

\S4method{dnbinom}{ad,ad,ad,missing,logical.}(x, size, prob, log)

\S4method{dnbinom}{num,num,num,missing,logical.}(x, size, prob, log)

\S4method{dnbinom}{osa,ANY,ANY,ANY,ANY}(x, size, prob, log)

\S4method{dnbinom}{simref,ANY,ANY,ANY,ANY}(x, size, prob, log)

\S4method{dpois}{ad,ad,logical.}(x, lambda, log = FALSE)

\S4method{dpois}{num,num,logical.}(x, lambda, log = FALSE)

\S4method{dpois}{osa,ANY,ANY}(x, lambda, log = FALSE)

\S4method{dpois}{simref,ANY,ANY}(x, lambda, log = FALSE)

\S4method{dgamma}{ad,ad,missing,ad.,logical.}(x, shape, scale, log)

\S4method{dgamma}{num,num,missing,num.,logical.}(x, shape, scale, log)

\S4method{dgamma}{osa,ANY,ANY,ANY,ANY}(x, shape, scale, log)

\S4method{dgamma}{simref,ANY,ANY,ANY,ANY}(x, shape, scale, log)

\S4method{pnorm}{ad,ad.,ad.,missing,missing}(q, mean, sd)

\S4method{pnorm}{num,num.,num.,missing,missing}(q, mean, sd)

\S4method{pgamma}{ad,ad,missing,ad.,missing,missing}(q, shape, scale)

\S4method{pgamma}{num,num,missing,num.,missing,missing}(q, shape, scale)

\S4method{ppois}{ad,ad,missing,missing}(q, lambda)

\S4method{ppois}{num,num,missing,missing}(q, lambda)

\S4method{pexp}{ad,ad.,missing,missing}(q, rate)

\S4method{pexp}{num,num.,missing,missing}(q, rate)

\S4method{pweibull}{ad,ad,ad.,missing,missing}(q, shape, scale)

\S4method{pweibull}{num,num,num.,missing,missing}(q, shape, scale)

\S4method{pbeta}{ad,ad,ad,missing,missing,missing}(q, shape1, shape2)

\S4method{pbeta}{num,num,num,missing,missing,missing}(q, shape1, shape2)

\S4method{qnorm}{ad,ad.,ad.,missing,missing}(p, mean, sd)

\S4method{qnorm}{num,num.,num.,missing,missing}(p, mean, sd)

\S4method{qgamma}{ad,ad,missing,ad.,missing,missing}(p, shape, scale)

\S4method{qgamma}{num,num,missing,num.,missing,missing}(p, shape, scale)

\S4method{qexp}{ad,ad.,missing,missing}(p, rate)

\S4method{qexp}{num,num.,missing,missing}(p, rate)

\S4method{qweibull}{ad,ad,ad.,missing,missing}(p, shape, scale)

\S4method{qweibull}{num,num,num.,missing,missing}(p, shape, scale)

\S4method{qbeta}{ad,ad,ad,missing,missing,missing}(p, shape1, shape2)

\S4method{qbeta}{num,num,num,missing,missing,missing}(p, shape1, shape2)

\S4method{besselK}{ad,ad,missing}(x, nu)

\S4method{besselK}{num,num,missing}(x, nu)

\S4method{besselI}{ad,ad,missing}(x, nu)

\S4method{besselI}{num,num,missing}(x, nu)

\S4method{besselJ}{ad,ad}(x, nu)

\S4method{besselJ}{num,num}(x, nu)

\S4method{besselY}{ad,ad}(x, nu)

\S4method{besselY}{num,num}(x, nu)

dbinom_robust(x, size, logit_p, log = FALSE)

dsn(x, alpha, log = FALSE)

dSHASHo(x, mu, sigma, nu, tau, log = FALSE)

dtweedie(x, mu, phi, p, log = FALSE)

dnbinom_robust(x, log_mu, log_var_minus_mu, log = FALSE)

dnbinom2(x, mu, var, log = FALSE)

dlgamma(x, shape, scale, log = FALSE)

logspace_add(logx, logy)

logspace_sub(logx, logy)

\S4method{dnorm}{ad,ad.,ad.,logical.}(x, mean = 0, sd = 1, log = FALSE)

\S4method{dnorm}{num,num.,num.,logical.}(x, mean = 0, sd = 1, log = FALSE)

\S4method{dnorm}{osa,ANY,ANY,ANY}(x, mean = 0, sd = 1, log = FALSE)

\S4method{dnorm}{simref,ANY,ANY,ANY}(x, mean = 0, sd = 1, log = FALSE)

\S4method{dlnorm}{ANY,ANY,ANY,ANY}(x, meanlog = 0, sdlog = 1, log = FALSE)

\S4method{dlnorm}{osa,ANY,ANY,ANY}(x, meanlog = 0, sdlog = 1, log = FALSE)

\S4method{dlnorm}{num,num.,num.,logical.}(x, meanlog = 0, sdlog = 1, log = FALSE)

\S4method{plogis}{advector,missing,missing,missing,missing}(q)

\S4method{qlogis}{advector,missing,missing,missing,missing}(p)

dcompois(x, mode, nu, log = FALSE)

dcompois2(x, mean, nu, log = FALSE)

\S4method{pbinom}{ad,ad,ad,missing,missing}(q, size, prob)

\S4method{pbinom}{num,num,num,missing,missing}(q, size, prob)

\S4method{dmultinom}{ad,ad.,ad,logical.}(x, size = NULL, prob, log = FALSE)

\S4method{dmultinom}{num,num.,num,logical.}(x, size = NULL, prob, log = FALSE)

\S4method{dmultinom}{osa,ANY,ANY,ANY}(x, size = NULL, prob, log = FALSE)

\S4method{dmultinom}{simref,ANY,ANY,ANY}(x, size = NULL, prob, log = FALSE)

\S4method{dmultinom}{ANY,ANY,ANY,ANY}(x, size = NULL, prob, log = FALSE)

\S4method{dcauchy}{ad,ad.,ad.,logical.}(x, location = 0, scale = 1, log = FALSE)

\S4method{dcauchy}{num,num.,num.,logical.}(x, location = 0, scale = 1, log = FALSE)

\S4method{dcauchy}{osa,ANY,ANY,ANY}(x, location = 0, scale = 1, log = FALSE)

\S4method{dcauchy}{simref,ANY,ANY,ANY}(x, location = 0, scale = 1, log = FALSE)
}
\arguments{
\item{x}{observation vector}

\item{rate}{parameter}

\item{log}{Logical; Return log density/probability?}

\item{shape}{parameter}

\item{scale}{parameter}

\item{size}{parameter}

\item{prob}{parameter}

\item{shape1}{parameter}

\item{shape2}{parameter}

\item{df1}{parameter}

\item{df2}{parameter}

\item{location}{parameter}

\item{df}{parameter}

\item{lambda}{parameter}

\item{q}{vector of quantiles}

\item{mean}{parameter}

\item{sd}{parameter}

\item{p}{parameter}

\item{nu}{parameter}

\item{logit_p}{parameter}

\item{alpha}{parameter}

\item{mu}{parameter}

\item{sigma}{parameter}

\item{tau}{parameter}

\item{phi}{parameter}

\item{log_mu}{parameter}

\item{log_var_minus_mu}{parameter}

\item{var}{parameter}

\item{logx}{Log-space input}

\item{logy}{Log-space input}

\item{meanlog}{Parameter; Mean on log scale.}

\item{sdlog}{Parameter; SD on log scale.}

\item{mode}{parameter}
}
\value{
In autodiff contexts an object of class \code{"advector"} is returned; otherwise a standard numeric vector.
}
\description{
The functions listed in this help page are all applicable for AD types.
Method dispatching follows a simple rule:
\emph{If at least one argument is an AD type then a special AD
implementation is selected. In all other cases a default
implementation is used} (typically that of the \bold{stats}
package).
Argument recycling follows the R standard (although wihout any warnings).
}
\details{
Specific documentation of the functions and arguments should be looked up elsewhere:
\itemize{
\item All S4 methods behave as the corresponding functions in the
\bold{stats} package. However, some arguments may not be
implemented in the AD case (e.g. \code{lower-tail}).
\item Other funtions behave as the corresponding TMB versions for
which documentation should be looked up online.
}
}
\section{Functions}{
\itemize{
\item \code{dexp(x = ad, rate = ad., log = logical.)}: AD implementation of \link[stats]{dexp}

\item \code{dexp(x = num, rate = num., log = logical.)}: Default method

\item \code{dexp(x = osa, rate = ANY, log = ANY)}: OSA implementation

\item \code{dexp(x = simref, rate = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dweibull(x = ad, shape = ad, scale = ad., log = logical.)}: AD implementation of \link[stats]{dweibull}

\item \code{dweibull(x = num, shape = num, scale = num., log = logical.)}: Default method

\item \code{dweibull(x = osa, shape = ANY, scale = ANY, log = ANY)}: OSA implementation

\item \code{dweibull(x = simref, shape = ANY, scale = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dbinom(x = ad, size = ad, prob = ad, log = logical.)}: AD implementation of \link[stats]{dbinom}

\item \code{dbinom(x = num, size = num, prob = num, log = logical.)}: Default method

\item \code{dbinom(x = osa, size = ANY, prob = ANY, log = ANY)}: OSA implementation

\item \code{dbinom(x = simref, size = ANY, prob = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dbeta(x = ad, shape1 = ad, shape2 = ad, ncp = missing, log = logical.)}: AD implementation of \link[stats]{dbeta}

\item \code{dbeta(x = num, shape1 = num, shape2 = num, ncp = missing, log = logical.)}: Default method

\item \code{dbeta(x = osa, shape1 = ANY, shape2 = ANY, ncp = ANY, log = ANY)}: OSA implementation

\item \code{dbeta(x = simref, shape1 = ANY, shape2 = ANY, ncp = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{df(x = ad, df1 = ad, df2 = ad, ncp = missing, log = logical.)}: AD implementation of \link[stats]{df}

\item \code{df(x = num, df1 = num, df2 = num, ncp = missing, log = logical.)}: Default method

\item \code{df(x = osa, df1 = ANY, df2 = ANY, ncp = ANY, log = ANY)}: OSA implementation

\item \code{df(x = simref, df1 = ANY, df2 = ANY, ncp = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dlogis(x = ad, location = ad., scale = ad., log = logical.)}: AD implementation of \link[stats]{dlogis}

\item \code{dlogis(x = num, location = num., scale = num., log = logical.)}: Default method

\item \code{dlogis(x = osa, location = ANY, scale = ANY, log = ANY)}: OSA implementation

\item \code{dlogis(x = simref, location = ANY, scale = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dt(x = ad, df = ad, ncp = missing, log = logical.)}: AD implementation of \link[stats]{dt}

\item \code{dt(x = num, df = num, ncp = missing, log = logical.)}: Default method

\item \code{dt(x = osa, df = ANY, ncp = ANY, log = ANY)}: OSA implementation

\item \code{dt(x = simref, df = ANY, ncp = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dnbinom(x = ad, size = ad, prob = ad, mu = missing, log = logical.)}: AD implementation of \link[stats]{dnbinom}

\item \code{dnbinom(x = num, size = num, prob = num, mu = missing, log = logical.)}: Default method

\item \code{dnbinom(x = osa, size = ANY, prob = ANY, mu = ANY, log = ANY)}: OSA implementation

\item \code{dnbinom(x = simref, size = ANY, prob = ANY, mu = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dpois(x = ad, lambda = ad, log = logical.)}: AD implementation of \link[stats]{dpois}

\item \code{dpois(x = num, lambda = num, log = logical.)}: Default method

\item \code{dpois(x = osa, lambda = ANY, log = ANY)}: OSA implementation

\item \code{dpois(x = simref, lambda = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dgamma(x = ad, shape = ad, rate = missing, scale = ad., log = logical.)}: AD implementation of \link[stats]{dgamma}

\item \code{dgamma(x = num, shape = num, rate = missing, scale = num., log = logical.)}: Default method

\item \code{dgamma(x = osa, shape = ANY, rate = ANY, scale = ANY, log = ANY)}: OSA implementation

\item \code{dgamma(x = simref, shape = ANY, rate = ANY, scale = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{pnorm(q = ad, mean = ad., sd = ad., lower.tail = missing, log.p = missing)}: AD implementation of \link[stats]{pnorm}

\item \code{pnorm(q = num, mean = num., sd = num., lower.tail = missing, log.p = missing)}: Default method

\item \code{pgamma(
  q = ad,
  shape = ad,
  rate = missing,
  scale = ad.,
  lower.tail = missing,
  log.p = missing
)}: AD implementation of \link[stats]{pgamma}

\item \code{pgamma(
  q = num,
  shape = num,
  rate = missing,
  scale = num.,
  lower.tail = missing,
  log.p = missing
)}: Default method

\item \code{ppois(q = ad, lambda = ad, lower.tail = missing, log.p = missing)}: AD implementation of \link[stats]{ppois}

\item \code{ppois(q = num, lambda = num, lower.tail = missing, log.p = missing)}: Default method

\item \code{pexp(q = ad, rate = ad., lower.tail = missing, log.p = missing)}: AD implementation of \link[stats]{pexp}

\item \code{pexp(q = num, rate = num., lower.tail = missing, log.p = missing)}: Default method

\item \code{pweibull(
  q = ad,
  shape = ad,
  scale = ad.,
  lower.tail = missing,
  log.p = missing
)}: AD implementation of \link[stats]{pweibull}

\item \code{pweibull(
  q = num,
  shape = num,
  scale = num.,
  lower.tail = missing,
  log.p = missing
)}: Default method

\item \code{pbeta(
  q = ad,
  shape1 = ad,
  shape2 = ad,
  ncp = missing,
  lower.tail = missing,
  log.p = missing
)}: AD implementation of \link[stats]{pbeta}

\item \code{pbeta(
  q = num,
  shape1 = num,
  shape2 = num,
  ncp = missing,
  lower.tail = missing,
  log.p = missing
)}: Default method

\item \code{qnorm(p = ad, mean = ad., sd = ad., lower.tail = missing, log.p = missing)}: AD implementation of \link[stats]{qnorm}

\item \code{qnorm(p = num, mean = num., sd = num., lower.tail = missing, log.p = missing)}: Default method

\item \code{qgamma(
  p = ad,
  shape = ad,
  rate = missing,
  scale = ad.,
  lower.tail = missing,
  log.p = missing
)}: AD implementation of \link[stats]{qgamma}

\item \code{qgamma(
  p = num,
  shape = num,
  rate = missing,
  scale = num.,
  lower.tail = missing,
  log.p = missing
)}: Default method

\item \code{qexp(p = ad, rate = ad., lower.tail = missing, log.p = missing)}: AD implementation of \link[stats]{qexp}

\item \code{qexp(p = num, rate = num., lower.tail = missing, log.p = missing)}: Default method

\item \code{qweibull(
  p = ad,
  shape = ad,
  scale = ad.,
  lower.tail = missing,
  log.p = missing
)}: AD implementation of \link[stats]{qweibull}

\item \code{qweibull(
  p = num,
  shape = num,
  scale = num.,
  lower.tail = missing,
  log.p = missing
)}: Default method

\item \code{qbeta(
  p = ad,
  shape1 = ad,
  shape2 = ad,
  ncp = missing,
  lower.tail = missing,
  log.p = missing
)}: AD implementation of \link[stats]{qbeta}

\item \code{qbeta(
  p = num,
  shape1 = num,
  shape2 = num,
  ncp = missing,
  lower.tail = missing,
  log.p = missing
)}: Default method

\item \code{besselK(x = ad, nu = ad, expon.scaled = missing)}: AD implementation of \link[base]{besselK}

\item \code{besselK(x = num, nu = num, expon.scaled = missing)}: Default method

\item \code{besselI(x = ad, nu = ad, expon.scaled = missing)}: AD implementation of \link[base]{besselI}

\item \code{besselI(x = num, nu = num, expon.scaled = missing)}: Default method

\item \code{besselJ(x = ad, nu = ad)}: AD implementation of \link[base]{besselJ}

\item \code{besselJ(x = num, nu = num)}: Default method

\item \code{besselY(x = ad, nu = ad)}: AD implementation of \link[base]{besselY}

\item \code{besselY(x = num, nu = num)}: Default method

\item \code{dbinom_robust()}: AD implementation

\item \code{dsn()}: AD implementation

\item \code{dSHASHo()}: AD implementation

\item \code{dtweedie()}: AD implementation

\item \code{dnbinom_robust()}: AD implementation

\item \code{dnbinom2()}: AD implementation

\item \code{dlgamma()}: AD implementation

\item \code{logspace_add()}: AD implementation

\item \code{logspace_sub()}: AD implementation

\item \code{dnorm(x = ad, mean = ad., sd = ad., log = logical.)}: AD implementation of \link[stats]{dnorm}

\item \code{dnorm(x = num, mean = num., sd = num., log = logical.)}: Default method

\item \code{dnorm(x = osa, mean = ANY, sd = ANY, log = ANY)}: OSA implementation

\item \code{dnorm(x = simref, mean = ANY, sd = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dlnorm(x = ANY, meanlog = ANY, sdlog = ANY, log = ANY)}: AD implementation of \link[stats]{dlnorm}.

\item \code{dlnorm(x = osa, meanlog = ANY, sdlog = ANY, log = ANY)}: OSA implementation.

\item \code{dlnorm(x = num, meanlog = num., sdlog = num., log = logical.)}: Default method.

\item \code{plogis(
  q = advector,
  location = missing,
  scale = missing,
  lower.tail = missing,
  log.p = missing
)}: Minimal AD implementation of \link[stats]{plogis}

\item \code{qlogis(
  p = advector,
  location = missing,
  scale = missing,
  lower.tail = missing,
  log.p = missing
)}: Minimal AD implementation of \link[stats]{qlogis}

\item \code{dcompois()}: Conway-Maxwell-Poisson. Calculate density.

\item \code{dcompois2()}: Conway-Maxwell-Poisson. Calculate density parameterized via the mean.

\item \code{pbinom(q = ad, size = ad, prob = ad, lower.tail = missing, log.p = missing)}: AD implementation of \link[stats]{pbinom}

\item \code{pbinom(q = num, size = num, prob = num, lower.tail = missing, log.p = missing)}: Default method

\item \code{dmultinom(x = ad, size = ad., prob = ad, log = logical.)}: AD implementation of \link[stats]{dmultinom}

\item \code{dmultinom(x = num, size = num., prob = num, log = logical.)}: Default method

\item \code{dmultinom(x = osa, size = ANY, prob = ANY, log = ANY)}: OSA implementation

\item \code{dmultinom(x = simref, size = ANY, prob = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

\item \code{dmultinom(x = ANY, size = ANY, prob = ANY, log = ANY)}: Default implementation that checks for invalid usage.

\item \code{dcauchy(x = ad, location = ad., scale = ad., log = logical.)}: AD implementation of \link[stats]{dcauchy}

\item \code{dcauchy(x = num, location = num., scale = num., log = logical.)}: Default method

\item \code{dcauchy(x = osa, location = ANY, scale = ANY, log = ANY)}: OSA implementation

\item \code{dcauchy(x = simref, location = ANY, scale = ANY, log = ANY)}: Simulation implementation. Modifies \code{x} and returns zero.

}}
\examples{
MakeTape( function(x) pnorm(x), x=numeric(5))$jacobian(1:5)
}
