library(Rcpp)

mod <- sourceCpp("complex.cpp", cacheDir=getwd(), cleanupCacheDir=TRUE, verbose=TRUE)

advector <- function(x) {
    if (inherits(x, "advector"))
        return (x)
    if (is.complex(x))
        stop("Invalid argument to 'advector' (lost class attribute?)")
    advec(x)
}
"Ops.advector" <- function(e1, e2) {
    if (missing(e2)) {
        if (.Generic=="-" || .Generic=="+") {
            e2 <- e1; e1 <- 0
        }
    }
    Arith2(advector(e1),
           advector(e2),
           .Generic)
}
"Math.advector" <- function(x) {
    x[] <- Math1(x, .Generic)
    x
}

## Matrix multiply is not a simple generic - overload entirely
"%*%" <- function(x, y) {
    if (inherits(x, "advector") || inherits(y, "advector"))
        matmul(x, y, method="atomic")
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
    y <- getValues(x)
    dim(y) <- dim(x)
    print(y, ...)
}

dnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
    r <- (x - mean) / sd
    ans <- - .5 * r * r - log(sqrt(2*pi)) - log(sd)
    if (log) ans else exp(ans)
}

MakeTape <- function(f, x, optimize=TRUE) {
    F <- new(adfun)
    F$start()
    ## Make sure to stop even in case of failure
    on.exit(F$stop())
    activate <- function(x)independent(advector(x))
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

f <- function(x) { cumsum(log(x)) }
F <- MakeTape(f, 1:5)
F$print()
F$eval(2:6)
F$jacobian(2:6)
F$jacfun() ## Transform
F$eval(2:6)
F$jacobian(2:6)

## Linreg
x <- 1:10
y <- 2*x+3+rnorm(length(x))
f <- function(theta) {
    a <- theta[1]
    b <- theta[2]
    sd <- theta[3]
    mu <- a*x+b
    -sum(dnorm(y, mu, sd, log=TRUE))
}
start <- c(2,3,1)+1
F <- MakeTape(f, start)
nlminb(start, F$eval, F$jacobian)

library(TMB)
## FIXME: Add data argument?
MakeADFun <- function(func, parameters, random=NULL, ...) {
    ## Make empty object
    obj <- TMB::MakeADFun(data=list(),
                          parameters=list(),
                          DLL="sourceCpp_2",
                          checkParameterOrder=FALSE)
    ## Overload and retape
    obj$env$MakeADFunObject <- function(data,parameters,...) {
        rcpp <- MakeTape(func, parameters)
        ans <- rcpp$ptrTMB()
        lgt <- lengths(parameters)
        par <- unlist(parameters, use.names=FALSE)
        names(par) <- rep(names(lgt), lgt)
        ans$DLL <- obj$env$DLL
        attr(ans$ptr, "par") <- par
        attr(ans, "rcpp") <- rcpp ## rcpp manages this ptr (no need for finalizer)
        ans
    }
    ## FIXME: Skip for now
    obj$env$MakeDoubleFunObject <- function(...)NULL
    obj$env$EvalDoubleFunObject <- function(...)NULL
    ## Change storage mode to double
    obj$env$parameters <- lapply(parameters, "+", 0.)
    obj$env$.random <- random
    if (!is.null(random))
        obj$env$type <- c(obj$env$type, "ADGrad")
    obj$retape()
    obj$par <- local(par[lfixed()],obj$env)
    obj
}

x <- 1:10
y <- 2*x+3+rnorm(length(x))
f <- function(pl) {
    a <- pl$a
    b <- pl$b
    sd <- exp(pl$logsd)
    mu <- a*x+b
    -sum(dnorm(y, mu, sd, log=TRUE))
}

p <- list(a=2,b=3,logsd=1)
obj <- MakeADFun(f, p, random="a")
obj$fn()
obj$gr()
obj$he()
