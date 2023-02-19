## Not export
"NACHECK<-" <- function(x, j, value) {
    if(!all(is.na(x[j]))) {
        stop("A simulation can only be assigned to once!")
    }
    x[j] <- value
    x
}
##' @describeIn simulation
simref <- function(n) {
    parent <- NULL
    value <- rep(NA, n)
    update <- function(val, j) {
        NACHECK(value, j) <<- val
    }
    asS4(structure(environment(), class="simref"))
}
##' @describeIn simulation
"dim<-.simref" <- function(x, value) {
    dim(x$value) <- value
    x
}
##' @describeIn simulation
"length.simref" <- function(x) length(x$value)
##' @describeIn simulation
"dim.simref" <- function(x) dim(x$value)
##' @describeIn simulation
is.array.simref <- function(x) is.array(x$value)
## FIXME: Not safe to modify x
##' @describeIn simulation
as.array.simref<-function(x,...){.x<-x;x<-x$value;.x$value<-NextMethod();.x}
##' @describeIn simulation
is.na.simref <- function(x)is.na(x$value)
##' @describeIn simulation
"[.simref" <- function(x, ...) {
    parent <- x
    x <- seq_along(parent$value)
    dim(x) <- dim(parent$value)
    ## Pointers into parent value
    i <- as.integer(NextMethod()) ## Handles '...'
    ## value of this object (repeat above subset)
    x <- parent$value
    value <- NextMethod() ## Handles '...'
    ## ----
    update <- function(val, j) {
        ## Update 'value' of this object
        NACHECK(value, j) <<- val
        ## Pass updated 'value' to parent
        parent$update( val, i[j] )
    }
    if (all(!is.na(value)))
        return (value)
    else
        asS4(structure(environment(), class="simref"))
}
##' @describeIn simulation
setMethod("show", "simref", function(object) {
    cat("class='simref'\n")
    print(object$value)
})
##' @describeIn simulation
"[<-.simref" <- function(x, ..., value) {
    cl <- match.call()
    cl <- head(cl, -1) ## Remove 'value'
    cl[[1]] <- as.name("[") ## Change generic
    cl[[2]] <- as.name("index")
    index <- seq_along(x)
    dim(index) <- dim(x)
    j <- as.integer(eval(cl))
    x$update(value, j)
    x
}
##' @describeIn simulation
Ops.simref <- function(e1, e2) {
    if (missing(e2)) stop("Not yet implemented")
    e1.unknown <- inherits(e1, "simref")
    e2.unknown <- inherits(e2, "simref")
    if (e1.unknown && e2.unknown) stop("Not yet implemented")
    getval <- function(e) if (inherits(e, "simref")) e$value else e
    value <- callGeneric(getval(e1), getval(e2))
    n <- length(value)
    s <- function(e, j) {
        e <- getval(e)
        if (length(e) == n)
            e[j]
        else if (length(e) == 1)
            e
        else
            rep(e, length.out=n)[j]
    }
    ## update(): Solve wrt unknown
    ##   ( e1 OP e2 ) [j] = val
    if (e1.unknown) {
        .Inverse <- c("+"="-", "-"="+", "*"="/", "/"="*")[.Generic]
    } else {
        ## e2.unknown
        .Inverse <- c("+"="-", "-"="-", "*"="/", "/"="/")[.Generic]
    }
    if (.Generic %in% c("-", "/")) {
        ## Swap args
        inverse <- function(e1, e2) .Primitive(.Inverse)(e2,e1)
    } else {
        inverse <- match.fun(.Inverse)
    }
    update <- function(val, j) {
        ve1 <- s(e1, j) ## Now length(e1) = length(j) = length(val)
        ve2 <- s(e2, j) ## Now length(e2) = length(j) = length(val)
        NACHECK(value, j) <<- val
        if (e1.unknown) e1$update(inverse(val, ve2), j)
        if (e2.unknown) e2$update(inverse(val, ve1), j)
    }
    asS4(structure(environment(), class="simref"))
}
##' @describeIn simulation
Math.simref <- function(x) {
    .Inverse <- c("exp"="log", "log"="exp")[.Generic]
    inverse <- match.fun(.Inverse)
    parent <- x
    value <- callGeneric(x$value)
    update <- function(val, j) {
        NACHECK(value, j) <<- val
        parent$update(inverse(val), j)
    }
    asS4(structure(environment(), class="simref"))
}
##' @describeIn simulation
t.simref <- function(x) {
    i <- seq_along(x)
    dim(i) <- dim(x)
    it <- t(i)
    itdim <- dim(it)
    y <- x[as.integer(it)]
    dim(y) <- itdim
    y
}
