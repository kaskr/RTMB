##' @describeIn ADvector Construct a new advector
advector <- function(x) {
    if (inherits(x, "advector"))
        return (x)
    if (is.complex(x))
        stop("Invalid argument to 'advector' (lost class attribute?)")
    ans <- advec(x)
    ## FIXME: Handling of attributes has to be carefully considered
    a <- attributes(x)
    if (!is.null(a)) {
        a$class <- c("advector", a$class)
        attributes(ans) <- a
    }
    ans
}
## Make this object AD interpretable if an AD context is active
## NOTE: Should be needed rarely!
## - Intended for objects that prevent simple S3 method dispatch
## - By defatult does *nothing* if *not* in an active AD context
## - To see what the 'magic' object would look like pass condition=TRUE
## Example 1: Starting out with a sparse matrix requires a little magic
## D <- Matrix::.symDiagonal(10)
## magic(D, TRUE) ## advector with attributes
magic <- function(x, condition = ad_context()) {
    if (!condition) return (x)
    if (is(x, "advector")) return (x)
    if (is(x, "sparseMatrix")) {
        x <- as(x, "adsparse")
        return (x)
    } else if (is.numeric(x)) {
        x <- advector(x)
        return (x)
    } else
        stop("'magic' does not know this object")
}
.Compare <- getGroupMembers("Compare")
##' @describeIn ADvector Binary operations
"Ops.advector" <- function(e1, e2) {
    if (compare_allow() && (.Generic %in% .Compare)) {
        return (NextMethod(getValues(e1), getValues(e2)))
    }
    if (missing(e2)) {
        if (.Generic=="-" || .Generic=="+") {
            e2 <- e1; e1 <- 0
        }
    }
    ans <- Arith2(advector(e1),
                  advector(e2),
                  .Generic)
    ## Object determining attrib of result
    e <- if (length(e2) > length(e1) || length(e2) == 0)
             e2
         else
             e1
    a <- attributes(e)
    if (!is.null(a)) {
        a$class <- "advector"
        attributes(ans) <- a
    }
    ans
}
##' @describeIn ADvector Unary operations
"Math.advector" <- function(x, ...) {
    x[] <- Math1(x, .Generic)
    x
}

##' @describeIn ADvector Makes \code{array(x)} work.
as.vector.advector <- function(x, mode = "any") {
    ## FIXME: Rcpp export 'as_advector' and use it
    asS4(structure(NextMethod(), class="advector"))
}
## unlist.advector <- function (x, recursive = TRUE, use.names = TRUE)  {
##     structure(NextMethod(), class="advector")
## }

##' @describeIn ADvector As \bold{base} version
aperm.advector <- function(a, perm, ...) {
    structure(NextMethod(), class="advector")
}
##' @describeIn ADvector As \bold{base} version
c.advector <- function(...) {
    structure(NextMethod(), class="advector")
}
##' @describeIn ADvector As \bold{base} version
"[.advector" <- function(x, ...) {
    structure(NextMethod(), class="advector")
}
##' @describeIn ADvector As \bold{base} version
"[<-.advector" <- function(x, ..., value) {
    value <- advector(value)
    NextMethod()
}
##' @describeIn ADvector As \bold{base} version
"[[.advector" <- function(x, ...) {
    structure(NextMethod(), class="advector")
}
##' @describeIn ADvector As \bold{base} version. Makes \code{outer(x,x,...)} work.
rep.advector <- function (x, ...) {
    structure(NextMethod(), class="advector")
}
##' @describeIn ADvector As \bold{base} version except \code{na.rm} not allowed.
sum.advector <- function(x, ..., na.rm) {
  if (na.rm) stop("'na.rm=TRUE' not implemented for AD sum")
  Reduce1(x, "+") + sum(...)
}
##' @describeIn ADvector As \bold{base} version except \code{na.rm} not allowed.
prod.advector <- function(x, ..., na.rm) {
  if (na.rm) stop("'na.rm=TRUE' not implemented for AD prod")
  Reduce1(x, "*") * prod(...)
}
## Make cov2cor() work. FIXME: Any unwanted side-effects with this?
##' @describeIn ADvector Makes \code{cov2cor()} work. FIXME: Any unwanted side-effects with this?
is.numeric.advector <- function(x) TRUE
## If an overload has issues we can patch it:
diff_patch <- base::diff.default
environment(diff_patch) <- local({unclass <- function(x)x; environment()})
diff.advector <- function (x, lag = 1L, differences = 1L, ...) {
    diff_patch(x, lag = 1L, differences = 1L, ...)
}

print.advector <- function (x, ...)  {
    cat("class='advector'\n")
    y <- .adv2num(x)
    print(y, ...)
}

## Helpers to autogenerate a numeric version of a function that is only available for AD types.
.adv2adv <- function(x) x
.adv2num <- function(x) {
    a <- attributes(x)
    a$class <- NULL
    ans <- getValues(x)
    attributes(ans) <- a
    ans
}
.anstype <- function(...) {
    if (any(unlist(lapply(list(...), inherits, "advector"))))
        .adv2adv
    else
        .adv2num
}

dmvnorm <- function(x, mu, Sigma, log=FALSE) {
    ## R convention is to have samples by row
    x <- t(as.matrix(x))
    mu <- t(as.matrix(mu))
    Sigma <- as.matrix(Sigma)
    d <- nrow(Sigma)
    x0 <- as.vector(x) - as.vector(mu)
    dim(x0) <- c(d, length(x0) / d)
    anstype <- .anstype(x0, Sigma)
    anstype( dmvnorm0(advector(x0), advector(Sigma), log) )
}

dgmrf <- function(x, mu, Q, log=FALSE) {
    if (!ad_context()) { ## Workaround: see C++ code 'gmrf0'
        F <- .MakeTape(function(...)advector(dgmrf(x,mu,Q,log)),numeric(0))
        return (F$eval(numeric(0)))
    }
    ## R convention is to have samples by row
    x <- t(as.matrix(x))
    mu <- t(as.matrix(mu))
    d <- attr(Q, "Dim")[1]
    x0 <- as.vector(x) - as.vector(mu)
    dim(x0) <- c(d, length(x0) / d)
    anstype <- .anstype(x0, Q)
    anstype( dgmrf0(advector(x0), as(Q, "adsparse"), log) )
}

## Low level version: Everything available
.MakeTape <- function(f, x) {
    F <- new(adfun)
    F$start()
    ## Make sure to stop even in case of failure
    on.exit(F$stop())
    activate <- function(x) {
        x <- advector(x)
        x[] <- independent(x)
        x
    }
    if (is.list(x)) {
        for (i in seq_along(x))
            x[[i]] <- activate(x[[i]])
    } else {
        x <- activate(x)
    }
    y <- f(x)
    dependent(y)
    F
}
## High level version: Not everything available
##' @describeIn Tape Generate a 'Tape' of an R function.
MakeTape <- function(f, x) {
    mod <- .MakeTape(f, x)
    .expose(mod)
}
.expose <- function(mod) {
    structure(
        function(x) {
            if (is.list(x))
                x <- do.call("c", x)
            if (inherits(x, "advector") && ad_context())
                mod$evalAD(x)
            else
                mod$eval(x)
        },
        methods = list(
            jacobian = mod$jacobian,
            optimize = mod$optimize,
            print = mod$print,
            jacfun = function() {
                .jacfun(mod)
            },
            laplace = function(random, sparse=TRUE, SPA=FALSE, ...) {
                .laplace(mod, random, sparse=sparse, SPA=SPA, ...)
            }
        ),
        class="Tape")
}
##' @describeIn Tape Get a tape method.
"$.Tape" <- function(x, name) attr(x, "methods")[[name]]
print.Tape <- function(x,...){
    cat("Object of class='Tape'\n")
    ptr <- environment(x)$mod$ptrTMB()$ptr
    info <- (.Call)(InfoADFunObject, ptr)
    txt <- paste0(" : ","R^",info$Domain, " -> " , "R^", info$Range, "\n")
    cat(txt)
    cat( c( "Methods:\n", paste0("$", names(attr(x,"methods")), "()\n")) )
}
.pointer <- function(mod) { ## FIXME: Is this safe Rcpp?
    env <- as.environment(mod)
    get(".pointer", envir = env)
}
.copy <- function(mod) {
    ans <- new(adfun)
    ans$copy(.pointer(mod))
    ans
}
.jacfun <- function(mod) {
    mod <- .copy(mod)
    mod$jacfun()
    .expose(mod)
}
.laplace <- function(mod, random, ...) {
    mod <- .copy(mod)
    random <- as.integer(random)
    cfg <- lapply(list(...), as.double)
    .transform(mod, "laplace", config=cfg,
               random_order=random, mustWork=1L)
    .transform(mod, "remove_random_parameters",
               random_order=random, mustWork=1L)
    .expose(mod)
}
.transform <- function(mod, method, ...) {
    ptr <- mod$ptrTMB()$ptr
    (.Call)(TransformADFunObject,
            f = ptr,
            control = list(method = as.character(method), ...))
}

##' @describeIn Tape Global configuration parameters of the tape (experts only!)
##' \bold{comparision} By default, AD comparision gives an error (\code{comparison="forbid"}).
##' This is the safe and recommended behaviour, because comparison is a non-differentiable operation. If you are building a tape that requires indicator functions e.g. \code{f(x)*(x<0)+g(x)*(x>=0)} then use \code{comparision="tape"} to add the indicators to the tape. A final option \code{comparision="allow"} exists for testing/illustration purposes. Do not use.
##' @param comparison Set behaviour of AD comparision (\code{">"},\code{"=="}, etc).
##' @param atomic Set behaviour of AD BLAS opererations (\code{"%*%"},...).
##' @param vectorize Enable/disable AD vectorized 'Ops' and 'Math'.
TapeConfig <- function(comparison = c("forbid", "tape", "allow"),
                       atomic = c("enable", "disable"),
                       vectorize = c("disable", "enable")) {
    comparison <- c(forbid=0L, tape=1L, allow=2L)[match.arg(comparison)]
    atomic <- c(enable=1L, disable=0L)[match.arg(atomic)]
    vectorize <- c(enable=1L, disable=0L)[match.arg(vectorize)]
    ans <- unlist(set_tape_config(comparison, atomic, vectorize))
    invisible(ans)
}

## FIXME: Add data argument?
MakeADFun <- function(func, parameters, random=NULL, map=list(), ADreport=FALSE, ...) {
    if (is.list(func))
        func <- attr(func, "func")
    ## Make empty object
    obj <- TMB::MakeADFun(data=list(),
                          parameters=parameters,
                          random=random,
                          map=map,
                          ADreport=FALSE,
                          checkParameterOrder=FALSE,
                          DLL="RTMB")
    ## Handling maps (copied and modified parList)
    parList <- function (parameters, par) {
        ans <- parameters
        ans[] <- lapply(ans, advector) ## (!)
        nonemp <- lengths(ans) > 0
        nonempindex <- which(nonemp)
        skeleton <- as.relistable(ans[nonemp])
        li <- relist(par, skeleton)
        reshape <- function(x) {
            if (is.null(attr(x, "map")))
                return(x)
            y <- attr(x, "shape")
            y <- advector(y) ## (!)
            f <- attr(x, "map")
            i <- which(f >= 0)
            y[i] <- x[f[i] + 1]
            y
        }
        for (i in seq(skeleton)) {
            ans[[nonempindex[i]]][] <- as.vector(li[[i]])
        }
        for (i in seq(ans)) {
            ans[[i]] <- reshape(ans[[i]])
        }
        ans
    }
    ## Overload and retape
    obj$env$MakeADFunObject <- function(data,parameters, reportenv, ADreport = FALSE,...) {
        mapfunc <- function(par) {
            ADREPORT_ENV$clear()
            pl <- parList(parameters, par)
            bias.correct <- ("TMB_epsilon_" %in% names(pl))
            if (bias.correct) {
                eps <- pl$TMB_epsilon_
                pl$TMB_epsilon_ <- NULL
            }
            ans <- func(pl)
            if (ADreport || bias.correct) {
                adrep <- do.call("c", lapply(ADREPORT_ENV$result(), advector) )
                if (length(adrep) == 0) adrep <- advector(numeric(0))
                if (bias.correct)
                    ans <- ans + sum(adrep * eps)
                else
                    ans <- adrep
            }
            ans
        }
        ## In TMB the 'par' is created on the C++ side
        obj$env$par <- unlist(parameters, use.names=FALSE)
        lgt <- lengths(parameters)
        names(obj$env$par) <- rep(names(lgt), lgt)
        rcpp <- .MakeTape(mapfunc, obj$env$par)
        ans <- rcpp$ptrTMB()
        ans$DLL <- obj$env$DLL
        attr(ans$ptr, "par") <- obj$env$par
        if (ADreport) {
            attr(ans$ptr, "range.names") <- ADREPORT_ENV$namevec()
            obj$env$ADreportDims <- ADREPORT_ENV$dims()
            ADREPORT_ENV$clear()
        }
        attr(ans, "rcpp") <- rcpp ## rcpp manages this ptr (no need for finalizer)
        ans
    }
    ## FIXME: Skip for now
    obj$env$MakeDoubleFunObject <- function(...)NULL
    obj$env$EvalDoubleFunObject <- function(...)NULL
    attr(obj$env$data, "func") <- func
    obj$env$ADreport <- ADreport
    obj$retape()
    obj$par <- obj$env$par[obj$env$lfixed()]
    obj
}

sdreport <- function(obj, ...) {
    sdreport_patch <- TMB::sdreport
    tmb_envir <- environment(sdreport_patch)
    env <- local({ MakeADFun <- RTMB::MakeADFun; environment() })
    parent.env(env) <- tmb_envir
    environment(sdreport_patch) <- env
    sdreport_patch(obj, ...)
}

reporter <- function() {
    ans <- list()
    report <- function(x) {
        nm <- deparse(substitute(x))
        ans[[nm]] <<- x
        NULL
    }
    result <- function() ans
    namevec <- function() {
        lgts <- lengths(ans)
        rep(names(lgts), lgts)
    }
    dims <- function() {
        getd <- function(x) { d <- dim(x); if(is.null(d)) length(x) else d}
        lapply(ans, getd)
    }
    clear <- function() ans <<- list()
    environment()
}
ADREPORT_ENV <- reporter()
REPORT_ENV <- reporter()
## User exported
ADREPORT <- ADREPORT_ENV$report
REPORT <- REPORT_ENV$report
