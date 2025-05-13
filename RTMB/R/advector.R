################################################################################
## This file contains:
## - The 'advector' and its core methods
## - The AD tape (MakeTape)
## - RTMB::MakeADFun
################################################################################

##' @describeIn ADvector Construct a new advector
##' @return Object of class \code{"advector"}.
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
##' Convert R object to AD
##'
##' Signify that this object should be given an AD interpretation if evaluated in an active AD context. Otherwise, keep object as is.
##' @details \code{AD} is a generic constructor, converting plain R structures to RTMB objects if in an autodiff context. Otherwise, it does nothing (and adds virtually no computational overhead).
##'
##' \code{AD} knows the following R objects:
##'
##' - Numeric objects from \pkg{base}, such as `numeric()`, `matrix()`, `array()`, are converted to class \link{advector} with other attributes kept intact.
##' - Complex objects from \pkg{base}, such as `complex()`, are converted to class \link{adcomplex}.
##' - Sparse matrices from \pkg{Matrix}, such as `Matrix()`, `Diagonal()`, are converted to \link{adsparse}.
##'
##' \code{AD} provides a reliable way to avoid problems with method dispatch when mixing operand types. For instance, sub assigning `x[i] <- y` may be problematic when `x` is numeric and `y` is `advector`. A prior statement `x <- AD(x)` solves potential method dispatch issues and can therefore be used as a reliable alternative to \link{ADoverload}.
##' @param x Object to be converted.
##' @param force Logical; Force AD conversion even if no AD context? (for debugging)
##' @examples
##' ## numeric object to AD
##' AD(numeric(4), force=TRUE)
##' ## complex object to AD
##' AD(complex(4), force=TRUE)
##' ## Convert sparse matrices (Matrix package) to AD representation
##' F <- MakeTape(function(x) {
##'   M <- AD(Matrix::Matrix(0,4,4))
##'   M[1,] <- x
##'   D <- AD(Matrix::Diagonal(4))
##'   D@x[] <- x
##'   M + D
##' }, 0)
##' F(2)
AD <- function(x, force = FALSE) {
    if (!force && !ad_context()) return (x)
    if (inherits(x, "advector")) return (x)
    if (inherits(x, "anysparse")) {
        x <- as(x, "adsparse")
    } else if (is.double(x) || is.integer(x) || is.logical(x)) {
        x <- advector(x)
    } else if (is.complex(x)) {
        x <- adcomplex(advector(Re(x)), advector(Im(x)))
    } else
        stop("'AD' does not know this object")
    x
}
.Compare <- getGroupMembers("Compare")
##' @describeIn ADvector Binary operations
"Ops.advector" <- function(e1, e2) {
    if (compare_allow() && (.Generic %in% .Compare)) {
        e1 <- getValues(advector(e1))
        e2 <- getValues(advector(e2))
        return (NextMethod())
    }
    if (missing(e2)) {
        if (.Generic=="-" || .Generic=="+") {
            e2 <- e1; e1 <- 0
        }
    }
    Arith2(advector(e1),
           advector(e2),
           .Generic)
}
##' @describeIn ADvector Unary operations
"Math.advector" <- function(x, ...) {
    Math1(x, .Generic, ...)
}

##' @describeIn ADvector Makes \code{array(x)} work.
as.vector.advector <- function(x, mode = "any") {
    ans <- NextMethod()
    if (is.list(ans)) lapply(ans, as_advector)
    else as_advector(ans)
}

##' @describeIn ADvector Convert to \link{ADcomplex}. Note that dimensions are dropped for consistency with base R.
as.complex.advector <- function(x, ...) {
    ans <- adcomplex(x)
    dim(ans) <- NULL ## For base R consistency
    ans
}

## unlist.advector <- function (x, recursive = TRUE, use.names = TRUE)  {
##     structure(NextMethod(), class="advector")
## }

##' @describeIn ADvector Equivalent of \link[base]{aperm}
aperm.advector <- function(a, perm, ...) {
    as_advector(NextMethod())
}
##' @describeIn ADvector Equivalent of \link[base]{c}. However note the limitation for mixed types: If `x` is an AD type, `c(x,1)` works while `c(1,x)` does not!
c.advector <- function(...) {
    ans <- unlist(lapply(list(...), advector))
    as_advector(ans)
}
##' @describeIn ADvector Equivalent of \link[base]{[}
"[.advector" <- function(x, ...) {
    as_advector(NextMethod())
}

## Extra RTMB overloads
xtra <- local({
    "[<-" <- function(x, ..., value) {
        if (inherits(value, "advector")) {
            if (is.numeric(x))
                x <- advector(x)
        }
        base::"[<-" (x, ..., value=value)
    }
    "diag<-" <- function(x, value) {
        if (inherits(value, "advector")) {
            if (is.numeric(x)) {
                x <- advector(x)
            } else
                if (inherits(x, "sparseMatrix")) {
                    x <- as(x, "adsparse")
                }
        }
        base::"diag<-" (x, value)
    }
    c <- function(...) {
        args <- list(...)
        if (any(unlist(lapply(args, inherits, "advector")))) {
            args <- lapply(args, advector)
            ans <- as_advector(unlist(args))
            return(ans)
        }
        base::"c" (...)
    }
    environment()
})
##' Enable extra RTMB convenience methods
##' @details Work around limitations in R's method dispatch system by overloading some selected primitives, currently:
##'
##' - Inplace replacement, so you can do `x[i] <- y` when `x` is numeric and `y` is AD.
##' - Mixed combine, so you can do e.g. `c(x, y)` when `x` numeric and `y` is AD.
##' - Diagonal assignment, so you can do `diag(x) <- y` when `x` is a numeric matrix and `y` is AD.
##'
##' In all cases, the result should be AD.
##' The methods are automatically **temporarily** attached to the search path (`search()`) when entering \link{MakeTape} or \link{MakeADFun}.
##' Alternatively, methods can be overloaded locally inside functions using e.g. `"[<-" <- ADoverload("[<-")`. This is only needed when using RTMB from a package.
##'
##' @examples
##' MakeTape(function(x) {print(search()); x}, numeric(0))
##' MakeTape(function(x) c(1,x), 1:3)
##' MakeTape(function(x) {y <- 1:3; y[2] <- x; y}, 1)
##' MakeTape(function(x) {y <- matrix(0,3,3); diag(y) <- x; y}, 1:3)
##' @param x Name of primitive to overload
##' @return Function representing the overload.
ADoverload <- function(x = c("[<-", "c", "diag<-")) {
    x <- match.arg(x)
    if (!ad_context())
        get(x, envir=baseenv())
    else
        get(x, envir=xtra, inherits=FALSE)
}
## For internal use
attachADoverloads <- function() {
    attached <- ( "AD-overloads" %in% search() )
    if (!attached)
        base::attach(xtra, length(search()),name="AD-overloads", warn=FALSE)
    NULL
}
detachADoverloads <- function(enable=TRUE, ...) {
    attached <- ( "AD-overloads" %in% search() )
    if (attached && !ad_context())
        base::detach("AD-overloads")
    NULL
}

##' @describeIn ADvector Equivalent of \link[base]{[<-}
"[<-.advector" <- function(x, ..., value) {
    value <- advector(value)
    NextMethod()
}
##' @describeIn ADvector Equivalent of \link[base]{[[}
"[[.advector" <- function(x, ...) {
    as_advector(NextMethod())
}
##' @describeIn ADvector Equivalent of \link[base]{rep}. Makes \code{outer(x,x,...)} work.
rep.advector <- function (x, ...) {
    as_advector(NextMethod())
}
##' @describeIn ADvector Equivalent of \link[base]{is.nan}. Check NaN status of a *constant* `advector` expression. If not constant throw an error.
is.nan.advector <- function(x) {
    if (!compare_allow() && any(getVariables(x)))
        stop("Can only determine NaN status of constant expressions")
    is.nan(getValues(x))
}
##' @describeIn ADvector Equivalent of \link[base]{is.finite}. Check finite status of a *constant* `advector` expression. If not constant throw an error.
is.finite.advector <- function(x) {
    if (!compare_allow() && any(getVariables(x)))
        stop("Can only determine finite status of constant expressions")
    is.finite(getValues(x))
}
##' @describeIn ADvector Equivalent of \link[base]{is.infinite}. Check infinity status of a *constant* `advector` expression. If not constant throw an error.
is.infinite.advector <- function(x) {
    if (!compare_allow() && any(getVariables(x)))
        stop("Can only determine infinity status of constant expressions")
    is.infinite(getValues(x))
}
##' @describeIn ADvector Equivalent of \link[base]{is.na}. Check NA status of an `advector`. NAs can only occur directly (as constants) or indirectly as the result of an operation with NA operands. For a tape built with non-NA parameters the NA status of any expression is constant and can therefore safely be used as part of the calculations. (assuming correct propagation of NAs via C-level arithmetic).
is.na.advector <- function(x) {
    is.na(getValues(x))
}
##' @describeIn ADvector Equivalent of \link[base]{sum}. \code{na.rm=TRUE} is allowed, but note that this feature assumes correct propagation of NAs via C-level arithmetic.
sum.advector <- function(x, ..., na.rm = FALSE) {
    if (na.rm) {
        x <- x[!is.na(x)]
    }
    Reduce1(x, "+") + sum(..., na.rm = na.rm)
}
##' @describeIn ADvector Equivalent of \link[base]{mean} except no arguments beyond `x` are supported.
mean.advector <- function(x, ...) {
    if (length(list(...))) stop("AD mean only works for single argument")
    sum(x) / length(x)
}
##' @describeIn ADvector Equivalent of \link[base]{prod}.
prod.advector <- function(x, ..., na.rm = FALSE) {
    if (na.rm) {
        x <- x[!is.na(x)]
    }
    Reduce1(x, "*") * prod(..., na.rm = na.rm)
}
##' @describeIn ADvector Equivalent of \link[base]{min}.
min.advector <- function(..., na.rm = FALSE) {
    x <- c(...)
    if (na.rm) {
        x <- x[!is.na(x)]
    }
    if (length(x) == 0) return (advector(Inf))
    Reduce1(x, "min")
}
##' @describeIn ADvector Equivalent of \link[base]{min}.
max.advector <- function(..., na.rm = FALSE) {
    x <- c(...)
    if (na.rm) {
        x <- x[!is.na(x)]
    }
    if (length(x) == 0) return (advector(-Inf))
    Reduce1(x, "max")
}

## Make cov2cor() work. FIXME: Any unwanted side-effects with this?
##' @describeIn ADvector Makes \code{cov2cor()} work. FIXME: Any unwanted side-effects with this?
is.numeric.advector <- function(x) TRUE
##' @describeIn ADvector Makes \code{as.numeric()} work.
as.double.advector <- function(x, ...) {
    ## Clear all attributes except class and preserve S4 bit
    attributes(x) <- attributes(x)["class"]
    x
}
##' @describeIn ADvector \link{Complex} operations are redirected to \link{adcomplex}.
##' @param z Complex (not allowed)
Complex.advector <- function(z) {
    callGeneric(adcomplex(z))
}
##' @describeIn ADvector Unimplemented \link{Summary} operations (currently \code{all} \code{any} \code{range}) will throw an error.
Summary.advector <- function(..., na.rm = FALSE)
    stop("'advector' does not allow operation ", sQuote(.Generic))
## If an overload has issues we can patch it:
diff_patch <- base::diff.default
environment(diff_patch) <- local({unclass <- function(x)x; environment()})
##' @describeIn ADvector Equivalent of \link[base]{diff}
##' @param lag As \link[base]{diff}
##' @param differences As \link[base]{diff}
diff.advector <- function (x, lag = 1L, differences = 1L, ...) {
    diff_patch(x, lag = lag, differences = differences, ...)
}

##' @describeIn ADvector Print method
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

## Low level version: Everything available
.MakeTape <- function(f, x) {
    F <- new(adfun)
    ## Start and attach overloads
    F$start()
    attachADoverloads()
    ## Make sure to stop and detach overloads (even in case of failure)
    on.exit({
        F$stop()
        detachADoverloads()
    })
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
    ## Result might be sparse matrix => Store pattern as an attribute
    Pattern <- NULL
    if (inherits(y, "adsparse")) {
        Pattern <- new("ngCMatrix", i=y@i, p=y@p, Dim=y@Dim)
        y <- y@x
    }
    y <- advector(y)
    dependent(y)
    if (is.null(Pattern))
        attr(F, "Dim") <- dim(y)
    else
        attr(F, "Pattern") <- Pattern
    F
}
## High level version: Not everything available
##' @describeIn Tape Generate a 'Tape' of an R function.
##' @return Object of class \code{"Tape"}.
MakeTape <- function(f, x) {
    f <- match.fun(f)
    mod <- .MakeTape(f, x)
    .expose(mod)
}
.expose <- function(mod) {
    Dim <- attr(mod, "Dim")
    Pattern <- attr(mod, "Pattern")
    output <- function(x) {
        if (!is.null(Dim))
            dim(x) <- Dim
        if (!is.null(Pattern)) {
            if (inherits(x, "advector"))
                x <- new("adsparse", x=x, i=Pattern@i, p=Pattern@p, Dim=Pattern@Dim)
            else
                x <- new("dgCMatrix", x=x, i=Pattern@i, p=Pattern@p, Dim=Pattern@Dim)
        }
        x
    }
    eval <- mod$eval ## cache
    evalAD <- mod$evalAD ## cache
    structure(
        function(x) {
            if (is.list(x))
                x <- do.call("c", x)
            if (ad_context()) {
                ## Note: Tape might contain references to an outer
                ## context (and we have no way to know), so we must
                ## choose AD evaluation regardless of class(x).
                x <- advector(x)
                output(evalAD(x))
            } else {
                output(eval(x))
            }
        },
        methods = list(
            jacobian = mod$jacobian,
            simplify = function(method=c("optimize", "eliminate")) {
                method <- match.arg(method)
                if (method == "optimize")
                    mod$optimize()
                else if (method == "eliminate")
                    mod$eliminate()
                else
                    stop("Unknown method")
            },
            print = function(depth=0) mod$print(as.integer(depth)),
            jacfun = function(sparse=FALSE) {
                if (!sparse)
                    .jacfun(mod)
                else
                    .spjacfun(mod)
            },
            atomic = function() {
                .atomic(mod)
            },
            laplace = function(random, sparse=TRUE, SPA=FALSE, ...) {
                .laplace(mod, random, sparse=sparse, SPA=SPA, ...)
            },
            newton = function(random, sparse=TRUE, ...) {
                .newton(mod, random, sparse=sparse, ...)
            },
            graph = function() {
                G <- get_graph(.pointer(mod))
                colnames(G) <- rownames(G) <- sub("Op","",colnames(G))
                G
            },
            data.frame = function() {
                get_df(.pointer(mod))
            },
            node = function(i) {
                mod <- .copy(mod)
                get_node(.pointer(mod), i)
                .expose(mod)
            },
            par = mod$domainvec
        ),
        class="Tape")
}
##' @describeIn Tape Get a tape method.
"$.Tape" <- function(x, name) {
    if (name == "methods") return (function()names(attr(x, "methods")))
    attr(x, "methods")[[name]]
}
##' @describeIn Tape Print method
##' @param ... Ignored
print.Tape <- function(x,...){
    cat("Object of class='Tape'\n")
    mod <- environment(x)$mod
    txt <- paste0(" : ","R^",mod$domain(), " -> " , "R^", mod$range(), "\n")
    cat(txt)
    ## cat( c( "Methods:\n", paste0("$", names(attr(x,"methods")), "()\n")) )
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
    jacdim <- c(mod$range(), mod$domain())
    mod <- .copy(mod)
    mod$jacfun()
    ans <- .expose(mod)
    environment(ans)$Dim <- jacdim
    ans
}
.spjacfun <- function(mod) {
    ptr <- .pointer(mod)
    jac <- SpJacFun(ptr)
    P <- new("dgTMatrix",
             i = jac@i,
             j = jac@j,
             x = as.double(seq_along(jac@i) - 1L),
             Dim = jac@Dim )
    P <- as(P, "CsparseMatrix") ## Permutes @x and calculates @i and @p
    RangeProj(jac@tape, P@x) ## Permute tape output to match new pattern
    mod <- new(adfun)
    e <- as.environment(mod)
    e$.pointer <- jac@tape
    ans <- .expose(mod)
    environment(ans)$Pattern <- new("ngCMatrix", i=P@i, p=P@p, Dim=P@Dim)
    ans
}
.atomic <- function(mod) {
    mod <- .copy(mod)
    mod$atomic()
    .expose(mod)
}
.laplace <- function(mod, random, ...) {
    mod <- .copy(mod)
    random <- as.integer(random)
    cfg <- lapply(list(...), as.double)
    mod$laplace(random, cfg)
    .expose(mod)
}
.newton <- function(mod, random, ...) {
    mod <- .copy(mod)
    random <- as.integer(random)
    cfg <- lapply(list(...), as.double)
    mod$newton(random, cfg)
    .expose(mod)
}
##' @describeIn Tape Global configuration parameters of the tape (experts only!)
##' \bold{comparison} By default, AD comparison gives an error
##' (\code{comparison="forbid"}).
##' This is the safe and recommended behaviour, because comparison is a
##' non-differentiable operation. If you are building a tape that
##' requires indicator functions e.g. \code{f(x)*(x<0)+g(x)*(x>=0)}
##' then use \code{comparison="tape"} to add the indicators to the
##' tape. A final option \code{comparison="allow"} exists for
##' testing/illustration purposes. Do not use.
##' @param comparison Set behaviour of AD comparison (\code{">"},\code{"=="}, etc).
##' @param atomic Set behaviour of AD BLAS operations (notably matrix multiply).
##' @param vectorize Enable/disable AD vectorized 'Ops' and 'Math'.
TapeConfig <- function(...,
                       comparison = c("NA", "forbid", "tape", "allow"),
                       atomic = c("NA", "enable", "disable"),
                       vectorize = c("NA", "disable", "enable")) {
    dotargs <- list(...)
    if (length(dotargs) && is.list(dotargs[[1]])) {
        if (length(dotargs) != 1)
            stop("List argument must be the *only* argument")
        if (!(missing(comparison) && missing(atomic) && missing(vectorize)))
            stop("List argument must be the *only* argument")
        return (do.call(set_tape_config, dotargs[[1]]))
    }
    args <- set_tape_config() ## Current settings
    ## Set comparison
    comparison <- match.arg(comparison)
    comparison <- c("NA"="NA", "forbid"="forbid", "tape"="taped", "allow"="allow")[comparison]
    args$compare <- comparison
    atomic <- match.arg(atomic)
    ## Set atomic
    atomic <- match.arg(atomic)
    atomic <- c("NA"="NA", "disable"="plain", "enable"="atomic")[atomic]
    args$matmul <- atomic
    args$mvnorm <- atomic
    ## Set vectorize
    vectorize <- match.arg(vectorize)
    vectorize <- c("NA"="NA", "disable"="plain", "enable"="vectorize")[vectorize]
    args$ops <- vectorize
    args$math <- vectorize
    args$sum <- vectorize
    ## Set other flags
    args[names(dotargs)] <- dotargs
    ## return
    invisible(do.call(set_tape_config, args))
}

##' @describeIn Tape Move a chunk of data from R to the tape by evaluating a normal R function (replaces TMB functionality 'DATA_UPDATE').
##' @examples
##' ## Taped access of an element of 'rivers' dataset
##' F <- MakeTape(function(i) DataEval( function(i) rivers[i] , i), 1 )
##' F(1)
##' F(2)
DataEval <- function(f, x) {
    if (ad_context())
        TapedEval(f, x)
    else
        f(x)
}

##' @describeIn Tape Extract tapes from a model object created by `MakeADFun`.
##' @param obj Output from `MakeADFun`
##' @param warn Give warning if `obj` was created using another DLL?
GetTape <- function(obj, name = c("ADFun", "ADGrad", "ADHess"), warn=TRUE) {
    name <- match.arg(name)
    if (name == "ADHess")
        env <- environment(obj$env$spHess)
    else
        env <- obj$env
    ADFun <- get(name, env, inherits=FALSE)
    stopifnot(is(ADFun$ptr ,"externalptr"))
    if (ADFun$DLL != "RTMB") {
        ok <- "getSetGlobalPtr" %in% names(getDLLRegisteredRoutines(ADFun$DLL)$.Call)
        if (!ok) {
            message("'getSetGlobalPtr' not found in '", ADFun$DLL, "'")
            stop("Please update TMB and recompile DLL '", ADFun$DLL, "'")
        }
        if (! exists("getSetGlobalPtr", getNamespace("RTMB")) ) {
            message("'getSetGlobalPtr' not found in 'RTMB'")
            stop("Please update TMB and recompile 'RTMB'")
        }
        fwrtmb <- .Call((getFramework))
        fwdll <- .Call(("getFramework"), PACKAGE=ADFun$DLL)
        if (!identical(fwrtmb, fwdll)) {
            info <- function(x) c(framework=x, attributes(x))
            null2na <- function(x) if (is.null(x)) NA else x
            df1 <- as.data.frame(info(fwrtmb))
            names(df1) <- names(info(fwrtmb))
            df2 <- as.data.frame(lapply(info(fwdll)[names(info(fwrtmb))], null2na ))
            names(df2) <- names(info(fwrtmb))
            df <- rbind(df1, df2)
            row.names(df) <- c('RTMB', ADFun$DLL)
            message("Note: DLL '", ADFun$DLL, "' is not binary compatible with 'RTMB'")
            print(df)
            message("- Compile with framework='TMBad'")
            message("- Compile with openmp=FALSE")
            message("- Compile with '-DTMBAD_INDEX_TYPE=uint64_t'")
            stop()
        }
        getSetGlobalPtr <- get("getSetGlobalPtr", getNamespace("RTMB"))
        RTMBptr <- .Call((getSetGlobalPtr), NULL)
        DLLptr <- .Call(("getSetGlobalPtr"), NULL, PACKAGE=ADFun$DLL)
        if (!identical(RTMBptr, DLLptr)) {
            if (warn) warning("Permanently changing the global pointer of DLL '", ADFun$DLL, "'")
            .Call(("getSetGlobalPtr"), RTMBptr, PACKAGE=ADFun$DLL)
        }
    }
    ans <- new(adfun)
    ans$copy(ADFun$ptr)
    .expose(ans)
}

## Visible bindings:
observation.name <- NULL
data.term.indicator <- NULL
data <- NULL
##' @describeIn TMB-interface Interface to \link[TMB]{MakeADFun}.
##' @details \link{MakeADFun} builds a TMB model object mostly compatible with the \pkg{TMB} package and with an almost identical interface.
##' The main difference in \pkg{RTMB} is that the objective function **and** the data is now given via a single argument \code{func}. Because \code{func} can be a *closure*, there is no need for an explicit data argument to \link{MakeADFun} (see examples).
##' @param func Function taking a parameter list (or parameter vector) as input.
##' @param parameters Parameter list (or parameter vector) used by \code{func}.
##' @param random As \link[TMB]{MakeADFun}.
##' @param profile As \link[TMB]{MakeADFun}.
##' @param integrate As \link[TMB]{MakeADFun}.
##' @param intern As \link[TMB]{MakeADFun}.
##' @param map As \link[TMB]{MakeADFun}.
##' @param ADreport As \link[TMB]{MakeADFun}.
##' @param silent As \link[TMB]{MakeADFun}.
##' @param ridge.correct Experimental
##' @param ... Passed to TMB
##' @return TMB model object.
##' @examples
##' ## Single argument vector function with numeric 'parameters'
##' fr <- function(x) {   ## Rosenbrock Banana function
##'     x1 <- x[1]
##'     x2 <- x[2]
##'     100 * (x2 - x1 * x1)^2 + (1 - x1)^2
##' }
##' obj <- MakeADFun(fr, numeric(2), silent=TRUE)
##' nlminb(c(-1.2, 1), obj$fn, obj$gr, obj$he)
MakeADFun <- function(func, parameters, random=NULL, profile=NULL, integrate=NULL, intern=FALSE, map=list(), ADreport=FALSE, silent=FALSE, ridge.correct=FALSE, ...) {
    setdata <- NULL
    if (is.list(func)) {
        setdata <- attr(func, "setdata")
        func <- attr(func, "func")
    }
    if (is.numeric(parameters)) {
        parnames <- names(formals(args(func)))
        if (length(parnames) != 1)
            stop("When 'parameters' is numeric 'func' must have a single argument")
        func1 <- func
        func <- function(p) do.call(func1, p)
        parameters <- structure(list(parameters), names=parnames)
    }
    ## Make empty object
    TMBArgs <- list(data=list(),
                    parameters=parameters,
                    random=random,
                    profile=profile,
                    map=map,
                    ADreport=FALSE,
                    checkParameterOrder=FALSE,
                    silent=silent,
                    integrate=NULL,
                    intern=FALSE,
                    ...)
    TMBArgs[duplicated(names(TMBArgs))] <- NULL
    TMBArgs$DLL <- "RTMB" ## Override if included in ...
    obj <- do.call(TMB::MakeADFun, TMBArgs)
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
            ## obj$env$data is normally empty, however some TMB
            ## functions (checkConsistency, others?) modifies the data
            ## and retapes. Let's make this data visible from
            ## objective:
            if (length(obj$env$data) > 0) {
                OBS_ENV$ans <- obj$env$data
            }
            pl <- parList(parameters, par)
            ## TMB is allowed to move data items to parameter list. We
            ## move such parmeters to the exchange environment OBS_ENV
            ## so can be obtained by 'OBS()' rather than 'parameters':
            if (length(obj$env$data) == 0) {
                setnm <- attr(obj$env$data, "setdata")
                if (!is.null(setnm)) {
                    for (nm in setnm) {
                        OBS_ENV$set(nm, pl[[nm]])
                        pl[[nm]] <- NULL
                    }
                }
            }
            bias.correct <- ("TMB_epsilon_" %in% names(pl))
            do.osa <- ("_RTMB_keep_" %in% names(pl))
            if (bias.correct) {
                eps <- pl$TMB_epsilon_
                pl$TMB_epsilon_ <- NULL
            }
            if (do.osa) {
                obn <- obj$env$observation.name
                dti <- obj$env$data.term.indicator
                x <- if (!is.null(pl[[obn]]))
                         pl[[obn]]
                     else
                         obj$env$data[[obn]]
                obs <- new("osa",
                           x = x,
                           keep = pl[[dti]])
                OBS_ENV$set(obn, obs)
                pl[[obn]] <- NULL
                pl[[dti]] <- NULL
            }
            ans <- func(pl)
            if (!do.osa) {
                asnum <- function(x) {
                    if (inherits(x, "advector"))
                        structure(getValues(x), dim=dim(x))
                    else x
                }
                ## Place OSA marked observations in obj
                obj$env$obs <- lapply(OBS_ENV$result(), asnum)
            }
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
        ## 'mapfunc' uses the 'exchange environments'. Clear before and after use:
        clear_all()
        on.exit(clear_all())
        rcpp <- .MakeTape(mapfunc, obj$env$par)
        if (!ADreport) {
            if (rcpp$range() != 1) {
                stop("'func' must return a *scalar* (forgot to sum?)")
            }
        }
        if (TMB::config(DLL="RTMB")$optimize.instantly)
            rcpp$optimize()
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
    ## OSA
    obj$env$observation.name <- observation.name
    obj$env$data.term.indicator <- data.term.indicator
    if (length(observation.name) && !is.null(data[[observation.name]]))
        obj$env$data[[observation.name]] <- data[[observation.name]]
    ## Simulate
    obj$simulate <- function(par=obj$env$last.par,...) {
        clear_all()
        on.exit(clear_all())
        p <- obj$env$parList(par=par)
        for (nm in obj$env$.random) {
            p[[nm]] <- simref2(p[[nm]], nm)
        }
        for (nm in names(obj$env$obs)) {
            obs <- simref2(obj$env$obs[[nm]], nm)
            OBS_ENV$set(nm, obs)
        }
        func(p)
        c(SIM_ENV$result(),
          REPORT_ENV$result())
    }
    ## Report
    obj$report <- function(par=obj$env$last.par,...) {
        clear_all()
        on.exit(clear_all())
        p <- obj$env$parList(par=par)
        func(p)
        REPORT_ENV$result()
    }
    ## FIXME: Skip for now
    obj$env$MakeDoubleFunObject <- function(...)NULL
    obj$env$EvalDoubleFunObject <- function(...)NULL
    attr(obj$env$data, "func") <- func
    attr(obj$env$data, "setdata") <- setdata
    obj$env$ADreport <- ADreport
    obj$env$profile <- profile
    obj$env$integrate <- integrate
    obj$env$intern <- intern
    obj$retape()
    obj$par <- obj$env$par[obj$env$lfixed()]
    if (ridge.correct) {
        obj$env$ridge.correct <- ridge.correct
        p <- ridge.correct
        p <- if (is.numeric(p)) p else .5
        obj$env$altHess <- altHessFun(obj, p)
    }
    obj
}

## 'Patch' a TMB function by overloading 'TMB::MakeADFun' with
## 'RTMB::MakeADFun', possibly accounting for TMB/RTMB argument
## differences.
## FIXME (?): The current patching takes place in RTMB namespace at
## RTMB install time. It follows, that changes to TMB
## (e.g. TMB::sdreport, TMB::oneStepPredict, ...) require
## re-installation of RTMB to take place ! We could alternatively
## apply the patches from inside the '.onLoad' function.
TMB_patch <- function(fun, ...) {
    environment(fun) <- new.env( parent = environment(fun) )
    environment(fun)$MakeADFun <- MakeADFun ## RTMB::MakeADFun
    dotargs <- list(...)
    formals(environment(fun)$MakeADFun)[names(dotargs)] <- dotargs
    fun
}

## Text substitution in function body (without changing environment)
bodysub <- function(fun, pattern, replacement, ...) {
    body(fun) <- parse(text=gsub(pattern,
                                 replacement,
                                 deparse(body(fun)),
                                 ...))[[1]]
    fun
}
sdreport_patch <- TMB_patch(TMB::sdreport)
sdreport_patch <- bodysub(sdreport_patch,
                          "\\(.Call\\(.*have_tmb_symbolic.*\\)\\)",
                          "(FALSE)")

##' @describeIn TMB-interface Interface to \link[TMB]{sdreport}.
##' @param obj TMB model object (output from \link{MakeADFun})
sdreport <- function(obj, ...) {
    sdreport_patch(obj, ...)
}

reporter <- function() {
    ans <- list()
    set <- function(nm, x) {
        ans[[nm]] <<- x
        NULL
    }
    get <- function(nm) ans[[nm]]
    report <- function(x) {
        nm <- deparse(substitute(x))
        set(nm, x)
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
##' @describeIn TMB-interface Can be used inside the objective function to report quantities for which uncertainties will be calculated by \link{sdreport}.
##' @param x Object to report
ADREPORT <- ADREPORT_ENV$report
##' @describeIn TMB-interface Can be used inside the objective function to report quantities via the model object using \code{obj$report()}.
REPORT <- REPORT_ENV$report

## Clear *all* exchange environments:
clear_all <- function() {
    OBS_ENV$clear()
    SIM_ENV$clear()
    REPORT_ENV$clear()
    ADREPORT_ENV$clear()
}

##' @describeIn TMB-interface Can be used to assign all parameter or data objects from a list inside the objective function.
##' @param warn Give a warning if overwriting an existing object?
getAll <- function(..., warn=TRUE) {
    fr <- parent.frame()
    x <- c(...)
    if (!is.list(x)) stop("'getAll' is for lists only")
    nm <- names(x)
    if (is.null(nm) || any(nm==""))
        stop("'getAll' is for *named* lists only")
    anyADvars <- FALSE ## For warn=TRUE case
    for (i in seq_along(x)) {
        if (warn) {
            if (!is.null(fr[[nm[i]]]))
                warning("Object '", nm[i], "' already defined")
            anyADvars <- anyADvars || inherits(x[[i]], "advector")
        }
        fr[[nm[i]]] <- x[[i]]
    }
    if (warn) {
        if (ad_context() && !anyADvars)
            warning("No active parameters found")
    }
    invisible(NULL)
}

