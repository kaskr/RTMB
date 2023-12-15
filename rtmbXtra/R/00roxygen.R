##' Log space gamma function
##'
##' Calculates `log(gamma(exp(x)))` with impoved stability for small x.
##' @title Log space gamma function
##' @param x AD vector
##' @return AD vector
##' @name logspace_gamma
##' @examples
##' F <- RTMB::MakeTape(logspace_gamma, 1)
##' identical(F(-1), lgamma(exp(-1)))
##' lgamma(exp(-1000))
##' F(-1000)
##' F$jacobian(-1000)
NULL

##' Logit transformed inverse cloglog
##'
##' Calculates `log(gamma(exp(x)))` with impoved stability for small x.
##' @title Logit transformed inverse cloglog
##' @param x AD vector
##' @return AD vector
##' @name logit_invcloglog
##' @examples
##' F <- RTMB::MakeTape(logit_invcloglog, 1)
NULL

##' Logit transformed pnorm
##'
##' Calculates `qlogis(pnorm(x)))` with impoved accuracy for small and large x.
##' @title Logit transformed pnorm
##' @param x AD vector
##' @return AD vector
##' @name logit_pnorm
##' @examples
##' F <- RTMB::MakeTape(logit_pnorm, 1)
##' abs(F(-1) - qlogis(pnorm(-1)))
##' F(-100)
##' F(100)
NULL
