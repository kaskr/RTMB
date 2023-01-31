##' Distributions and special functions for which AD is implemented
##'
##' The functions listed in this help page are all applicable for AD types.
##' Method dispatching follows a simple rule:
##' \emph{If at least one argument an AD type then a special AD
##' implementation is selected. In all other cases a default
##' implementation is used} (typically that of the \bold{stats}
##' package).
##' Argument recycling follows the R standard (although wihout any warnings).
##'
##' @param x observation vector
##' @param q vector of quantiles
##' @param rate parameter
##' @param shape parameter
##' @param scale parameter
##' @param size parameter
##' @param prob parameter
##' @param logit_p parameter
##' @param shape1 parameter
##' @param shape2 parameter
##' @param df1 parameter
##' @param df2 parameter
##' @param location parameter
##' @param alpha parameter
##' @param df parameter
##' @param mu parameter
##' @param sigma parameter
##' @param nu parameter
##' @param tau parameter
##' @param phi parameter
##' @param p parameter
##' @param var parameter
##' @param log_mu parameter
##' @param log_var_minus_mu parameter
##' @param lambda parameter
##' @param mean parameter
##' @param sd parameter
##' @param scale parameter
##' @param log Logical; Return log density/probability?
##' @rdname Distributions
##' @name Distributions
##' @examples
##' MakeTape( function(x) pnorm(x), x=numeric(5))$jacobian(1:5)
NULL

##' RTMB: R bindings for 'TMB'
##'
##' The package 'RTMB' provides a native R interface for *a subset of*
##' 'TMB' so you can avoid coding in C++.  'RTMB' only affects the
##' 'TMB' function 'MakeADFun' that builds the objective function. Once
##' 'MakeADFun' has been invoked, everything else is \emph{exactly the same}
##' and \emph{models run as fast} as if coded in C++.
##'
##' @rdname RTMB-package
##' @name RTMB-package
##' @aliases RTMB-package RTMB
##' @docType package
##' @author Kasper Kristensen
##'
##' Maintainer: <kaskr@@dtu.dk>
NULL
