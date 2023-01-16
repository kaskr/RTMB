library(TMB)
library(Rcpp)

mod <- sourceCpp("example.cpp", verbose=TRUE)

advector <- function(x) {
    if (inherits(x, "advector"))
        return (x)
    if (inherits(x, "externalptr"))
        ptr <- x
    else
        ptr <- advec(x)
    structure(list(ptr=ptr), class="advector")
}
"Ops.advector" <- function(e1, e2) {
    if (missing(e2)) {
        if (.Generic=="-" || .Generic=="+") {
            e2 <- e1; e1 <- 0
        }
    }
    advector(Arith2(advector(e1)$ptr,
                    advector(e2)$ptr,
                    .Generic))
}
"Math.advector" <- function(x) {
    advector(Math1(x$ptr, .Generic))
}

## 'as.vector', 'length', 'dim', 'dim<-' are necessary for array(x,...)
as.vector.advector <- function (x, mode = "any") {
    if (mode=="list")
        return (structure(lapply(seq_along(x), function(i)x[i]), class="adlist"))
    x
}
length.advector <- function(x) Length(x$ptr)
dim.advector <- function(x) attr(x, "Dim")
"dim<-.advector" <- function(x, value) {
    attr(x, "Dim") <- value
    x
}
## Make lapply work
##is.list.advector <- function(x)FALSE
as.list.advector <- function(x)structure(lapply(seq_along(x), function(i)x[i]), class="adlist")

## Make as.matrix work
is.matrix.advector <- function(x) length(dim(x))==2L


## FIXME: externalptr is unique (OK) and so are the attributes (NOT OK):
##   x <- array(advector(1:25),c(5,5))
##   y <- x
##   dim(y) <- NULL ## Changes dim(x) !!!

c.advector <- function(...) {
    dotargs <- list(...)
    ans <- advector(numeric(0))
    for (x in dotargs) AppendInplace(ans$ptr, advector(x)$ptr)
    advector(ans)
}   
print.advector <- function(x) {
    invisible(Display(x$ptr))
    if (!is.null(dim(x)))
        cat("Dim = ","(",paste(dim(x),collapse=", "),")\n")
}

"[.advector" <- function(x, ...) {
    ##ArrayApply("[", x, ...)
    j <- seq_len(length(x))
    if (!is.null(dim(x)))
        j <- array(j, dim=dim(x))
    j <- j[...]
    newdim <- dim(j)
    structure(advector(Subset(x$ptr, j - 1L)), Dim=newdim)
}
DeepCopy <- function(x)c(advector(numeric(0)),x)
"[<-.advector" <- function(x, i, value) {
    ## Deep copy only if necessary
    if (.Call("sourceCpp_1_MaybeShared", x))
        x$ptr <- DeepCopy(x)$ptr
    ## Note: i can be 'missing', 'integer', '2 col matrix' or 'logical' !!!
    j <- seq_len(length(x))
    if (!is.null(dim(x)))
        j <- array(j, dim=dim(x))
    j <- j[i]
    SubsetAssign(x$ptr, j - 1L, advector(value)$ptr)
    x
}

"[[.advector" <- function(x,...) x[...]
## Generalizing previous: Apply index reordering function
ArrayApply <- function(F, x, ...) {
    F <- match.fun(F)
    ## Apply function with 'x' replaced by its indices
    j <- seq_len(length(x))
    if (!is.null(dim(x)))
        j <- array(j, dim=dim(x))
    j <- F(j, ...)
    newdim <- dim(j)
    structure(advector(Subset(x$ptr, j - 1L)), Dim=newdim)
}
t.advector <- function(x) ArrayApply("t", x)
aperm.advector <- function (a, perm, ...) { ArrayApply("aperm", a, perm, ...) }
diff.advector <- function(x) tail(x,-1)-head(x,-1)
sum.advector <- function(x, na.rm)advector(Reduce1(x$ptr, "sum"))
prod.advector <- function(x, na.rm)advector(Reduce1(x$ptr, "prod"))
MakeTape <- function(f, x, optimize=TRUE) {
    F <- new(adfun)
    F$start()
    ## Make sure to stop even in case of failure
    on.exit(F$stop())
    x <- advector(x)
    Independent(x$ptr)
    y <- f(x)
    Dependent(y$ptr)
    F
}

## =============== Distributions
## Are generally not generic - must overload
## FIXME: Dispatch must be based on all aguments x, mean, sd
dnorm <- function(x, mean = 0, sd = 1, log = FALSE) UseMethod("dnorm", mean)
dnorm.default <- stats::dnorm
dnorm.advector <- function(x, mean = 0, sd = 1, log = FALSE) {
    r <- (x - mean) / sd
    ans <- - .5 * r * r - log(sqrt(2*pi)) - log(sd)
    if (log) ans else exp(ans)
}


## =============== Example
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
