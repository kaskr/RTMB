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
        x <- as(x, "generalMatrix")
        ipdim <- attributes(x)[c("i", "p", "Dim")]
        x <- advector(x@x)
        attributes(x) <- c(attributes(x), ipdim)
        return (x)
    } else if (is.numeric(x)) {
        x <- advector(x)
        return (x)
    } else
        stop("'magic' does not know this object")
}
"Ops.advector" <- function(e1, e2) {
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
"Math.advector" <- function(x) {
    x[] <- Math1(x, .Generic)
    x
}

## Matrix multiply is not a simple generic - overload entirely
"%*%" <- function(x, y) {
    if (inherits(x, "advector") || inherits(y, "advector")) {
        x <- as.matrix(x)
        y <- as.matrix(y)
        matmul(advector(x), advector(y), method="atomic")
    }
    else
        .Primitive("%*%")(x, y)
}

## Make array(x) work. Also affects as.list(x)
as.vector.advector <- function(x, mode = "any") {
    structure(NextMethod(), class="advector")
}
unlist.advector <- function (x, recursive = TRUE, use.names = TRUE)  {
    structure(NextMethod(), class="advector")
}
aperm.advector <- function(a, perm, ...) {
    structure(NextMethod(), class="advector")
}
c.advector <- function(...) {
    structure(NextMethod(), class="advector")
}
"[.advector" <- function(x, ...) {
    structure(NextMethod(), class="advector")
}
"[<-.advector" <- function(x, ..., value) {
    value <- advector(value)
    structure(NextMethod(), class="advector")
}
"[[.advector" <- function(x, ...) {
    structure(NextMethod(), class="advector")
}
## Make outer(x,x,'') work
rep.advector <- function (x, ...) {
    structure(NextMethod(), class="advector")
}
sum.advector <- function(x, na.rm)Reduce1(x, "+")
prod.advector <- function(x, na.rm)Reduce1(x, "*")

## If an overload has issues we can patch it:
diff_patch <- base::diff.default
environment(diff_patch) <- local({unclass <- function(x)x; environment()})
diff.advector <- function (x, lag = 1L, differences = 1L, ...) {
    diff_patch(x, lag = 1L, differences = 1L, ...)
}

## Danger?
apply <- function (X, MARGIN, FUN, ...)  {
    ans <- base::apply(X, MARGIN, FUN, ...)
    if (inherits(X, "advector"))
        class(ans) <- "advector"
    ans
}
sapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    ans <- base::sapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
    if (inherits(X, "advector"))
        class(ans) <- "advector"
    ans
}

print.advector <- function (x, ...)  {
    cat("class='advector'\n")
    y <- .adv2num(x)
    print(y, ...)
}

dnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
    r <- (x - mean) / sd
    ans <- - .5 * r * r - log(sqrt(2*pi)) - log(sd)
    if (log) ans else exp(ans)
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
        F <- MakeTape(function(...)advector(dgmrf(x,mu,Q,log)),numeric(0))
        return (F$eval(numeric(0)))
    }
    ## R convention is to have samples by row
    x <- t(as.matrix(x))
    mu <- t(as.matrix(mu))
    d <- attr(Q, "Dim")[1]
    x0 <- as.vector(x) - as.vector(mu)
    dim(x0) <- c(d, length(x0) / d)
    anstype <- .anstype(x0, Q)
    anstype( dgmrf0(advector(x0), magic(Q, TRUE), log) )
}

MakeTape <- function(f, x, optimize=TRUE) {
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
        rcpp <- MakeTape(mapfunc, obj$env$par)
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
    obj$par <- local(par[lfixed()], obj$env)
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
