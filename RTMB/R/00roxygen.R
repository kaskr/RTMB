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

##' Interface to TMB
##'
##' @rdname TMB-interface
##' @name TMB-interface
##' @examples
##' data(rivers)
##' f <- function(p) { -sum(dnorm(rivers, p$mu, p$sd, log=TRUE)) }
##' obj <- MakeADFun(f, list(mu=0, sd=1), silent=TRUE)
##' opt <- nlminb(obj$par, obj$fn, obj$gr)
##' sdreport(obj)
##' ## 'REML trick'
##' obj2 <- MakeADFun(f, list(mu=0, sd=1), random="mu", silent=TRUE)
##' opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
##' sdreport(obj2) ## Compare with sd(rivers)
NULL

##' The AD vector and its methods
##'
##' An \code{advector} is a class used behind the scenes to replace
##' normal R numeric objects during automatic differentiation. An
##' \code{advector} has a 'temporary lifetime' and therefore you do not
##' \emph{see} / \emph{need to know} it as a normal user.
##'
##' An AD vector (class='advector') is an atomic R vector of 'codes'
##' that are internally interpretable as 'AD scalars'. A substantial
##' part of R's existing S3 matrix and array functionality can be
##' re-used for AD vectors.
##'
##' @param x numeric or advector
##' @param a advector with dimension attribute
##' @param e1 advector
##' @param e2 advector
##' @param perm Permutation as in \code{aperm}
##' @param mode FIXME might not be handled correctly by \code{as.vector}
##' @param na.rm Must be FALSE (default)
##' @param value Replacement value implicitly converted to AD
##' @param ... Additional arguments
##' @rdname ADvector
##' @name ADvector
##' @examples
##' x <- advector(1:9)
##' a <- array(x, c(3,3))  ## as an array
##' outer(x, x, "+") ## Implicit via 'rep'
##' rev(x)           ## Implicit via '['
NULL

##' AD matrix methods (sparse and dense)
##'
##' Matrices (**base** package) and sparse matrices (**Matrix** package) can be used inside the \code{RTMB} objective function as part of the calculations. Behind the scenes these R objects are converted to AD representations when needed. AD objects have a temporary lifetime, so you probably won't see them / need to know them. The only important thing is which *methods* work for the objects.
##'
##' @param x matrix (sparse or dense)
##' @param y matrix (sparse or dense)
##' @rdname ADmatrix
##' @name ADmatrix
##' @examples
##' F <- MakeTape(function(x) matrix(1:9,3,3) %*% x, numeric(3))
##' F$jacobian(1:3)
##' F <- MakeTape(function(x) Matrix::expm(matrix(x,2,2)), numeric(4))
##' F$jacobian(1:4)
NULL

##' AD aware numeric constructors
##'
##' These base functions, that construct numeric like objects, have been modified to return AD objects when used inside an active AD context. This is necessary due to limitations in R's method dispatch machinery (\code{[<-} cannot dispatch on the \code{value} argument).
##'
##' @rdname ADconstruct
##' @name ADconstruct
##' @examples
##' func <- function(t) {
##'   x <- numeric(5)
##'   print(class(x))
##'   t
##' }
##' invisible(func(1))            ## 'numeric'
##' invisible(MakeTape(func, 1))  ## 'advector'
NULL

##' AD apply functions
##'
##' These **base** apply methods have been modified to keep the AD class attribute (which would otherwise be lost).
##'
##' @rdname ADapply
##' @name ADapply
##' @examples
##' F <- MakeTape(function(x) apply(matrix(x,2,2), 2, sum), numeric(4))
##' F$jacobian(1:4)
NULL

##' The AD tape
##'
##' The AD tape as an R function
##'
##' A 'Tape' is a representation of a function that accepts \emph{fixed size} numeric input and returns \emph{fixed size} numeric output.
##' The tape can be constructed using \code{F <- MakeTape(f, x)} where \code{f} is a standard \emph{differentiable} R function (or more precisely: One using only functions that are documented to work for AD types).
##' Having constructed a tape F, a number of methods are available:
##'
##' Evaluation:
##' - Normal function evaluation 'F(x)' for numeric input.
##' - AD evaluation 'F(x)' as part of other tapes.
##' - Jacobian calculations using 'F$jacobian(x)'.
##'
##' Transformation:
##' - Jacobian function transformation 'F$jacfun()'
##' - Laplace approximation transformation 'F$laplace()'. Also allows Saddle Point Approximation (SPA).
##'
##' @param f R function
##' @param x numeric vector
##' @param name Name of a tape method
##' @rdname Tape
##' @name Tape
##' @examples
##' F <- MakeTape(prod, numeric(3))
##' show(F)
##' F$print()
##' H <- F$jacfun()$jacfun() ## Hessian tape
##' show(H)
##' \dontrun{
##' #### Handy way to plot the graph of F
##' if (require(igraph)) {
##'    G <- graph_from_adjacency_matrix(F$graph())
##'    plot(G, vertex.size=17, layout=layout_as_tree)
##' }
##' }
NULL

##' Distributions and special functions for which AD is implemented
##'
##' The functions listed in this help page are all applicable for AD types.
##' Method dispatching follows a simple rule:
##' \emph{If at least one argument is an AD type then a special AD
##' implementation is selected. In all other cases a default
##' implementation is used} (typically that of the \bold{stats}
##' package).
##' Argument recycling follows the R standard (although wihout any warnings).
##'
##' Specific documentation of the functions and arguments should be looked up elsewhere:
##' - All S4 methods behave as the corresponding functions in the
##'   \bold{stats} package. However, some arguements may not be
##'   implemented in the AD case (e.g. \code{lower-tail}).
##' - Other funtions behave as the corresponding TMB versions for
##'   which documentation should be looked up online.
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
##' @param mode parameter
##' @param sd parameter
##' @param scale parameter
##' @param log Logical; Return log density/probability?
##' @rdname Distributions
##' @name Distributions
##' @examples
##' MakeTape( function(x) pnorm(x), x=numeric(5))$jacobian(1:5)
NULL

##' Multivariate Gaussian densities
##'
##' Multivariate Gaussian densities
##'
##' @rdname MVgauss
##' @name MVgauss
NULL

##' Recursive quantile residuals
##'
##' OSA residuals are computed using the function
##' \code{oneStepPredict}. For this to work, you need to mark the
##' observation inside the objective function using the \link{OBS}
##' function. Thereafter, residual calculation is as simple as
##' \code{oneStepPredict(obj)}. However, you probably want specify a
##' \code{method} to use.
##'
##' @rdname OSA-residuals
##' @name OSA-residuals
##' @examples
##' set.seed(1)
##' rw <- cumsum(.5*rnorm(20))
##' obs <- rpois(20, lambda=exp(rw))
##' func <- function(p) {
##'   obs <- OBS(obs) ## Mark 'obs' for OSA calculation on request
##'   ans <- 0
##'   jump <- c(p$rw[1], diff(p$rw))
##'   ans <- ans - sum(dnorm(jump, sd=p$sd, log=TRUE))
##'   ans <- ans - sum(dpois(obs, lambda=exp(p$rw), log=TRUE))
##'   ans
##' }
##' obj <- MakeADFun(func,
##'                  parameters=list(rw=rep(0,20), sd=1),
##'                  random="rw")
##' nlminb(obj$par, obj$fn, obj$gr)
##' res <- oneStepPredict(obj,
##'                       method="oneStepGeneric",
##'                       discrete=TRUE,
##'                       range=c(0,Inf))$residual
NULL

##' Simulation
##'
##' An RTMB objective function can be run in 'simulation mode' where standard likelihood evaluation is replaced by corresponding random number generation. This facilitates automatic simulation under some restrictions. Simulations can be obtained directly from the model object by \code{obj$simulate()} or used indirectly via \link{checkConsistency}.
##'
##' In simulation mode all log density evaluation, involving either random effects or observations, is interpreted as probability assignment.
##'
##' \bold{direct vs indirect} Assignments can be 'direct' as for example
##'
##' \code{dnorm(u, log=TRUE)      ## u ~ N(0, 1)}
##'
##' or 'indirect' as in
##'
##' \code{dnorm(2*u+1, log=TRUE)  ## u ~ N(-1, .25)}
##'
##' Indirect assignment works for a limited set of easily invertible functions - see \code{methods(class="simref")}.
##'
##' \bold{Simulation order} Note that probability assignments are sequential: All information required to draw a new variable must already be simulated.
##' Vectorized assignment implicitly occurs elementwise from left to right.
##' For example the assignment
##'
##' \code{dnorm(diff(u), log=TRUE)}
##'
##' is not valid without a prior assignment of \code{u[1]}, e.g.
##'
##' \code{dnorm(u[1], log=TRUE)}
##'
##' \bold{Supported distributions} Assignment must use supported density functions. I.e.
##'
##' \code{dpois(N, exp(u), log=TRUE)}
##'
##' cannot be replaced by
##'
##' \code{N * u - exp(u)}
##'
##' The latter will have no effect in simulation mode (the simulation will be \code{NA}).
##'
##' \bold{Return value} Note that when in simulation mode, the density functions all return zero. The actual simulation is written to the input argument by reference. This is very unlike standard R semantics.
##' @rdname Simulation
##' @name Simulation
##' @examples
##' s <- simref(4)
##' s2 <- 2 * s[1:2] + 1
##' s2[] <- 7
##' s ## 3 3 NA NA
##' ## Random walk
##' func <- function(p) {
##'   u <- p$u
##'   ans <- -dnorm(u[1], log=TRUE) ## u[1] ~ N(0,1)
##'   ans <- ans - sum(dnorm(diff(u), log=TRUE)) ## u[i]-u[i-1] ~ N(0,1)
##' }
##' obj <- MakeADFun(func, list(u=numeric(20)), random="u")
##' obj$simulate()
NULL

setClass("advector") ## Virtual class
setClass("adsparse",
         slots=c(x="advector", i="integer", p="integer", Dim="integer"))
## Helpers to setMethod
## Match numeric like objects that are 'AD castable' via advector()
setClassUnion("num", c("array", "numeric", "logical"))
setClassUnion("num.", c("num", "missing"))
setClassUnion("ad",  c("advector", "num"))
setClassUnion("ad.", c("advector", "num."))
setClassUnion("logical.", c("logical", "missing"))
## For OSA residuals
setClass("osa", list(x="ad", keep="ad"))
## For simulation
setClass("simref")
