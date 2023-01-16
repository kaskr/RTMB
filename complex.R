library(Rcpp)

mod <- sourceCpp("complex.cpp", verbose=TRUE)

advector <- function(x) {
    if (inherits(x, "advector"))
        return (x)
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
    x <- advector(x)
    x <- Independent(x)
    y <- f(x)
    Dependent(y)
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
