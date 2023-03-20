## Not export
"NACHECK<-" <- function(x, j, value) {
    if(!all(is.na(x[j]))) {
        stop("A simulation can only be assigned to once!")
    }
    x[j] <- value
    x
}
##' @describeIn Simulation Construct \code{simref}
##' @param n Length
simref <- function(n) {
    parent <- NULL
    value <- rep(NA, n)
    num.missing <- length(value)
    complete <- function() (num.missing == 0)
    ## 'finalize' is called when last element of simulation has been written.
    ## A custom 'finalize' could be added to do:
    ## - store the simulation in a safe place and
    ## - override this 'simref' object with the simulation
    finalize <- function(value) NULL
    forward.update <- function(j) {
        value[j]
    }
    reverse.update <- function(val, j) {
        NACHECK(value, j) <<- val
        num.missing <<- num.missing - length(j)
        if (complete()) finalize(value)
    }
    asS4(structure(environment(), class="simref"))
}
##' @describeIn Simulation Equivalent of \link[base]{dim<-}
##' @param x Object of class 'simref'
##' @param value Replacement (numeric)
"dim<-.simref" <- function(x, value) {
    dim(x$value) <- value
    x
}
##' @describeIn Simulation Equivalent of \link[base]{length}
"length.simref" <- function(x) length(x$value)
##' @describeIn Simulation Equivalent of \link[base]{dim}
"dim.simref" <- function(x) dim(x$value)
##' @describeIn Simulation Equivalent of \link[base]{is.array}
is.array.simref <- function(x) is.array(x$value)
##' @describeIn Simulation Equivalent of \link[base]{is.matrix}
is.matrix.simref <- function(x) is.matrix(x$value)
## FIXME: Not safe to modify x
##' @describeIn Simulation Equivalent of \link[base]{as.array}
as.array.simref<-function(x,...){.x<-x;x<-x$value;.x$value<-NextMethod();.x}
##' @describeIn Simulation Equivalent of \link[base]{is.na}
is.na.simref <- function(x)is.na(x$value)
##' @describeIn Simulation Equivalent of \link[base]{[}
##' @param ... Extra arguments
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
    forward.update <- function(j) {
        value[j] <<- parent$forward.update( i[j] )
        ## Implicit return: value[j]
    }
    reverse.update <- function(val, j) {
        ## Update 'value' of this object
        NACHECK(value, j) <<- val
        ## Pass updated 'value' to parent
        parent$reverse.update( val, i[j] )
    }
    if (all(!is.na(value)))
        return (value)
    else
        asS4(structure(environment(), class="simref"))
}
##' @describeIn Simulation
setMethod("show", "simref", function(object) {
    cat("class='simref'\n")
    print(object$value)
})
##' @describeIn Simulation Equivalent of \link[base]{[<-}
"[<-.simref" <- function(x, ..., value) {
    cl <- match.call()
    cl <- cl[-length(cl)] ## Remove 'value'
    cl[[1]] <- as.name("[") ## Change generic
    cl[[2]] <- as.name("index")
    index <- seq_along(x)
    dim(index) <- dim(x)
    j <- as.integer(eval(cl))
    x$reverse.update(value, j)
    x
}
##' @describeIn Simulation Equivalent of \link[base]{Ops}
##' @param e1 First argument
##' @param e2 Second argument
Ops.simref <- function(e1, e2) {
    call <- sys.call()
    call[[1]] <- as.name(.Generic)
    if (missing(e2)) {
        if (.Generic == "+")
            return(e1)
        else if (.Generic == "-")
            return(0-e1)
        else
            stop("Cannot use '", .Generic, "' with single argument")
    }
    ## 'unknown' here means 'not completely known' i.e. there might be some elements known
    e1.unknown <- inherits(e1, "simref")
    e2.unknown <- inherits(e2, "simref")
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
    ## reverse.update(): Solve wrt unknown
    ##   ( e1 OP e2 ) [j] = val
    ## e1.unknown
    .Inverse1 <- c("+"="-", "-"="+", "*"="/", "/"="*")[.Generic]
    ## e2.unknown
    .Inverse2 <- c("+"="-", "-"="-", "*"="/", "/"="/")[.Generic]
    ## Known case?
    if (is.na(.Inverse1) || is.na(.Inverse2))
        stop("Class 'simref' does not know generic '", .Generic, "'")
    if (.Generic %in% c("-", "/")) {
        ## Swap args
        inverse1 <- function(e1, e2) .Primitive(.Inverse1)(e2, e1)
        inverse2 <- function(e1, e2) .Primitive(.Inverse2)(e2, e1)
    } else {
        inverse1 <- match.fun(.Inverse1)
        inverse2 <- match.fun(.Inverse2)
    }
    forward.update <- function(j) {
        if (e1.unknown) {
            if (length(e1) == n)
                e1$forward.update(j)
            else
                e1$forward.update()
        }
        if (e2.unknown) {
            if (length(e2) == n)
                e2$forward.update(j)
            else
                e2$forward.update()
        }
        value[j] <<- match.fun(.Generic)(s(e1, j), s(e2, j))
    }
    if (e1.unknown && e2.unknown) {
        reverse.update <- function(val, j) {
            if (length(e1) != length(e2)) stop("Unexpected")
            for (i in seq_along(j)) {
                e1$forward.update(j[i])
                e2$forward.update(j[i])
                NACHECK(value, j[i]) <<- val[i]
                if (is.na(e1$value[j[i]] && is.na(e2$value[j[i]]))) {
                    expr <- deparse(substitute((x)[j], list(x=call,j=j[i])))
                    stop("Implicit simulation failed.\n",
                         "Cannot update operands of '",
                         expr,
                         "' when both are missing")
                }
                if (is.na(e1$value[j[i]]))
                    e1$reverse.update(inverse1(val[i], e2$value[j[i]]), j[i])
                if (is.na(e2$value[j[i]]))
                    e2$reverse.update(inverse2(val[i], e1$value[j[i]]), j[i])
            }
        }
    } else {
        reverse.update <- function(val, j) {
            ve1 <- s(e1, j) ## Now length(e1) = length(j) = length(val)
            ve2 <- s(e2, j) ## Now length(e2) = length(j) = length(val)
            NACHECK(value, j) <<- val
            if (e1.unknown) e1$reverse.update(inverse1(val, ve2), j)
            if (e2.unknown) e2$reverse.update(inverse2(val, ve1), j)
        }
    }
    asS4(structure(environment(), class="simref"))
}
##' @describeIn Simulation Equivalent of \link[base]{Math}
Math.simref <- function(x, ...) {
    .Inverse <- c("exp"="log", "log"="exp")[.Generic]
    if (is.na(.Inverse))
        stop("Class 'simref' does not know generic '", .Generic, "'")
    inverse <- match.fun(.Inverse)
    parent <- x
    value <- callGeneric(x$value)
    forward.update <- function(j) {
        value[j] <<- match.fun(.Generic)(parent$forward.update(j))
    }
    reverse.update <- function(val, j) {
        NACHECK(value, j) <<- val
        parent$reverse.update(inverse(val), j)
    }
    asS4(structure(environment(), class="simref"))
}
##' @describeIn Simulation Equivalent of \link[base]{t}
t.simref <- function(x) {
    i <- seq_along(x)
    dim(i) <- dim(x)
    it <- t(i)
    itdim <- dim(it)
    y <- x[as.integer(it)]
    dim(y) <- itdim
    y
}
##' @describeIn Simulation Equivalent of \link[base]{diff}
##' @param lag As \link[base]{diff}
##' @param differences As \link[base]{diff}
diff.simref <- function (x, lag = 1L, differences = 1L, ...) {
    diff_patch(x, lag = 1L, differences = 1L, ...)
}
##' @describeIn Simulation \link{Summary} operations are not invertible and will throw an error.
##' @param na.rm Ignored
Summary.simref <- function(..., na.rm = FALSE)
    stop("Class 'simref' does not allow operation '", .Generic, "' (argument must be fully simulated)")

dGenericSim <- function(.Generic, x, ..., log) {
    if (!log) stop("'simref' is for *log* density evaluation only")
    substring(.Generic, 1, 1) <- "r"
    rfun <- match.fun(.Generic)
    x[] <- rfun(length(x), ...)
    rep(0, length(x))
}

## internal 'simref' constructor with a custom 'finalize':
## - Overrides 'simref' in its frame when complete
## - And stores the simulation in SIM_ENV
SIM_ENV <- reporter()
simref2 <- function(x, name) {
    s <- simref(length(x))
    dim(s) <- dim(x)
    s$name <- name
    s$finalize <- function(value) {
        for (e in sys.frames()) {
            if (inherits(e[[name]], "simref")) {
                e[[name]] <- e[[name]]$value
                SIM_ENV$set(name, value)
                break
            }
        }
    }
    s
}

## Special set/unset to handle cases where TMB moves items from data
## list to parameters list. 'RTMB::MakeADFun' needs to know the names
## of these items.
setdata <- function(obj) {
    obj$env$data[names(obj$env$obs)] <- obj$env$obs
    attr(obj$env$data, "setdata") <- names(obj$env$obs)
    NULL
}
unsetdata <- function(obj) {
    obj$env$data[] <- NULL
}

##' @describeIn TMB-interface Interface to \link[TMB]{checkConsistency}.
##' @param fast Pass `observation.name` to `TMB` ?
checkConsistency <- function(obj, fast=TRUE, ...) {
    ## ######
    checkConsistency_patch <- TMB::checkConsistency
    tmb_envir <- environment(checkConsistency_patch)
    env <- local({
        MakeADFun <- RTMB::MakeADFun
        formals(MakeADFun)$data <- list() ## checkConsistency passes data here...
        formals(MakeADFun)$func <- as.name("data") ## And we redirect it here
        environment()
    })
    parent.env(env) <- tmb_envir
    environment(checkConsistency_patch) <- env
    ## ######
    args <- list(obj, ...)
    if (fast) {
        supported <- "observation.name" %in% names(formals(checkConsistency_patch))
        if (!supported) {
            warning("'fast' selected but not supported")
        } else {
            setdata(obj)
            on.exit(unsetdata(obj))
            args$observation.name <- names(obj$env$obs)
        }
    }
    do.call("checkConsistency_patch", args)
}

## Internal: Not export
rtweedie <- function (n, mu, phi, p) {
    ok <- all( (p > 1) & (p < 2) )
    if (!ok) stop("'p' must be in the interval (1, 2)")
    lambda <- mu^(2 - p)/(phi * (2 - p))
    alpha <- (2 - p)/(1 - p)
    gam <- phi * (p - 1) * mu^(p - 1)
    N <- rpois(n, lambda = lambda)
    rgamma(n, shape = -N * alpha, scale = gam)
}
