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
        e1 <- getValues(advector(e1))
        e2 <- getValues(advector(e2))
        return (NextMethod())
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

##' @describeIn ADvector Equivalent of \link[base]{aperm}
aperm.advector <- function(a, perm, ...) {
    asS4(structure(NextMethod(), class="advector"))
}
##' @describeIn ADvector Equivalent of \link[base]{c}. However note the limitation for mixed types: If `x` is an AD type, `c(x,1)` works while `c(1,x)` does not!
c.advector <- function(...) {
    ans <- unlist(lapply(list(...), advector))
    asS4(structure(ans, class = "advector"))
}
##' @describeIn ADvector Equivalent of \link[base]{[}
"[.advector" <- function(x, ...) {
    asS4(structure(NextMethod(), class="advector"))
}
##' @describeIn ADvector Equivalent of \link[base]{[<-}
"[<-.advector" <- function(x, ..., value) {
    value <- advector(value)
    NextMethod()
}
##' @describeIn ADvector Equivalent of \link[base]{[[}
"[[.advector" <- function(x, ...) {
    asS4(structure(NextMethod(), class="advector"))
}
##' @describeIn ADvector Equivalent of \link[base]{rep}. Makes \code{outer(x,x,...)} work.
rep.advector <- function (x, ...) {
    structure(NextMethod(), class="advector")
}
##' @describeIn ADvector Equivalent of \link[base]{sum} except \code{na.rm} not allowed.
sum.advector <- function(x, ..., na.rm) {
  if (na.rm) stop("'na.rm=TRUE' not implemented for AD sum")
  Reduce1(x, "+") + sum(...)
}
##' @describeIn ADvector Equivalent of \link[base]{prod} except \code{na.rm} not allowed.
prod.advector <- function(x, ..., na.rm) {
  if (na.rm) stop("'na.rm=TRUE' not implemented for AD prod")
  Reduce1(x, "*") * prod(...)
}
## Make cov2cor() work. FIXME: Any unwanted side-effects with this?
##' @describeIn ADvector Makes \code{cov2cor()} work. FIXME: Any unwanted side-effects with this?
is.numeric.advector <- function(x) TRUE
##' @describeIn ADvector \link{Complex} operations are not allowed and will throw an error.
##' @param z Complex (not allowed)
Complex.advector <- function(z)
    stop("'advector' does not allow complex operations")
##' @describeIn ADvector Non differentiable \link{Summary} operations (e.g. \code{min} \code{max}) are not allowed and will throw an error.
Summary.advector <- function(..., na.rm = FALSE)
    stop("'advector' does not allow operation ", sQuote(.Generic))
## If an overload has issues we can patch it:
diff_patch <- base::diff.default
environment(diff_patch) <- local({unclass <- function(x)x; environment()})
##' @describeIn ADvector Equivalent of \link[base]{diff}
##' @param lag As \link[base]{diff}
##' @param differences As \link[base]{diff}
diff.advector <- function (x, lag = 1L, differences = 1L, ...) {
    diff_patch(x, lag = 1L, differences = 1L, ...)
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

##' @describeIn Distributions Conway-Maxwell-Poisson. Calculate density.
dcompois <- function(x, mode, nu, log = FALSE) {
    if (inherits(x,"osa"))    return (dGenericOSA( "dcompois" , x=x, mode=mode, nu=nu, log=log ))
    if (inherits(x,"simref")) return (dGenericSim( "dcompois" , x=x, mode=mode, nu=nu, log=log ))
    loglambda <- nu * log(mode)
    ans <- x * loglambda - nu * lfactorial(x)
    ans <- ans - compois_calc_logZ(loglambda, nu)
    if (log) ans else exp(ans)
}
## internal (simref)
rcompois <- function(n, mode, nu) {
    loglambda <- nu * log(mode)
    loglambda <- rep(loglambda, length.out=n)
    mapply(distr_rcompois, loglambda, nu) ## mapply re-cycles
}

##' @describeIn Distributions Conway-Maxwell-Poisson. Calculate density parameterized via the mean.
dcompois2 <- function(x, mean, nu, log = FALSE) {
    if (inherits(x,"osa"))    return (dGenericOSA( "dcompois2" , x=x, mean=mean, nu=nu, log=log ))
    if (inherits(x,"simref")) return (dGenericSim( "dcompois2" , x=x, mean=mean, nu=nu, log=log ))
    logmean <- log(mean)
    loglambda <- compois_calc_loglambda(logmean, nu)
    ans <- x * loglambda - nu * lfactorial(x)
    ans <- ans - compois_calc_logZ(loglambda, nu)
    if (log) ans else exp(ans)
}
## internal (simref)
rcompois2 <- function(n, mean, nu) {
    logmean <- log(mean)
    loglambda <- getValues(compois_calc_loglambda(advector(logmean), advector(nu)))
    loglambda <- rep(loglambda, length.out=n)
    mapply(distr_rcompois, loglambda, nu) ## mapply re-cycles
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
    y <- advector(y)
    dependent(y)
    F
}
## High level version: Not everything available
##' @describeIn Tape Generate a 'Tape' of an R function.
MakeTape <- function(f, x) {
    f <- match.fun(f)
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
            atomic = function() {
                .atomic(mod)
            },
            laplace = function(random, sparse=TRUE, SPA=FALSE, ...) {
                .laplace(mod, random, sparse=sparse, SPA=SPA, ...)
            },
            graph = function() {
                G <- get_graph(.pointer(mod))
                colnames(G) <- rownames(G) <- sub("Op","",colnames(G))
                G
            }
        ),
        class="Tape")
}
##' @describeIn Tape Get a tape method.
"$.Tape" <- function(x, name) attr(x, "methods")[[name]]
##' @describeIn Tape Print method
##' @param ... Ignored
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
.atomic <- function(mod) {
    mod <- .copy(mod)
    mod$atomic()
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
TapeConfig <- function(comparison = c("forbid", "tape", "allow"),
                       atomic = c("enable", "disable"),
                       vectorize = c("disable", "enable")) {
    comparison <- c(forbid=0L, tape=1L, allow=2L)[match.arg(comparison)]
    atomic <- c(enable=1L, disable=0L)[match.arg(atomic)]
    vectorize <- c(enable=1L, disable=0L)[match.arg(vectorize)]
    ans <- unlist(set_tape_config(comparison, atomic, vectorize))
    invisible(ans)
}

## Visible bindings:
observation.name <- NULL
data.term.indicator <- NULL
data <- NULL
##' @describeIn TMB-interface Interface to \link[TMB]{MakeADFun}.
##' @param func Function taking a parameter list (or parameter vector) as input.
##' @param parameters Parameter list (or parameter vector) used by \code{func}.
##' @param random As \link[TMB]{MakeADFun}.
##' @param map As \link[TMB]{MakeADFun}.
##' @param ADreport As \link[TMB]{MakeADFun}.
##' @param silent As \link[TMB]{MakeADFun}.
##' @param ... Passed to TMB
##' @examples
##' ## Single argument vector function with numeric 'parameters'
##' obj <- MakeADFun(function(x)-sum(dnorm(x, log=TRUE)), 1:10)
MakeADFun <- function(func, parameters, random=NULL, map=list(), ADreport=FALSE, silent=FALSE,...) {
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
    obj <- TMB::MakeADFun(data=list(),
                          parameters=parameters,
                          random=random,
                          map=map,
                          ADreport=FALSE,
                          checkParameterOrder=FALSE,
                          silent=silent,
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
                ## Place OSA marked observations in obj
                obj$env$obs <- OBS_ENV$result()
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
    obj$retape()
    obj$par <- obj$env$par[obj$env$lfixed()]
    obj
}

##' @describeIn TMB-interface Interface to \link[TMB]{sdreport}.
##' @param obj TMB model object (output from \link{MakeADFun})
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
    for (i in seq_along(x)) {
        if (warn) {
            if (!is.null(fr[[nm[i]]]))
                warning("Object '", nm[i], "' already defined")
        }
        fr[[nm[i]]] <- x[[i]]
    }
    invisible(NULL)
}
