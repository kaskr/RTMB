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
MakeADFun <- function(func, parameters, random=NULL, map=list(), ...) {
    ## Make empty object
    obj <- TMB::MakeADFun(data=list(),
                          parameters=parameters,
                          random=random,
                          map=map,
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
    obj$env$MakeADFunObject <- function(data,parameters,...) {
        mapfunc <- function(par) {
            pl <- parList(parameters, par)
            func(pl)
        }
        rcpp <- MakeTape(mapfunc, obj$env$par)
        ans <- rcpp$ptrTMB()
        ans$DLL <- obj$env$DLL
        attr(ans$ptr, "par") <- obj$env$par
        attr(ans, "rcpp") <- rcpp ## rcpp manages this ptr (no need for finalizer)
        ans
    }
    ## FIXME: Skip for now
    obj$env$MakeDoubleFunObject <- function(...)NULL
    obj$env$EvalDoubleFunObject <- function(...)NULL
    obj$retape()
    obj
}
