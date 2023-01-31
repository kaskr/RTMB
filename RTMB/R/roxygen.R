##' Distributions and special functions for which AD is implemeted
##'
##' The functions listed in this help page are all applicable to AD types.
##' Method dispatching follows a simple rule:
##' - If at least one argument is of an AD type then a special AD implementation is selected. In all other cases we fall back on a default implementation (typically that of the stats package).
##' Argument recycling follows the R standard (although wihout any warnings).
##'
##' @name Distributions
##' @example MakeTape(function(x) dnorm(x), x=1:5)$print()
NULL
