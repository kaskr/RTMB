##' @describeIn OSA-residuals Calculate the residuals. See documentation of \code{TMB::}\link[TMB]{oneStepPredict}.
##' @param obj TMB model object (output from \code{MakeADFun})
##' @param observation.name Auto detected - use the default
##' @param data.term.indicator Auto detected - use the default
##' @param ... Passed to \code{TMB::}\link[TMB]{oneStepPredict} - \bold{please carefully read the documentation}, especially the \code{method} argument.
oneStepPredict <- function(obj,
                           observation.name=names(obj$env$obs)[1],
                           data.term.indicator="_RTMB_keep_",
                           ...) {
    ## FIXME: Check against valid OSA names
    ## names(obj$env$obs)
    obj$env$data[[observation.name]] <- obj$env$obs[[observation.name]]
    obj$env$observation.name <- observation.name
    obj$env$data.term.indicator <- data.term.indicator
    on.exit({
        obj$env$data[[observation.name]] <- NULL
        obj$env$observation.name <- NULL
    })
    oneStepPredict_patch <- TMB::oneStepPredict
    tmb_envir <- environment(oneStepPredict_patch)
    env <- local({
        MakeADFun <- RTMB::MakeADFun
        formals(MakeADFun)$data <- list() ## Dummy
        formals(MakeADFun)$observation.name <- list() ## Dummy
        formals(MakeADFun)$data.term.indicator <- list() ## Dummy
        formals(MakeADFun)$func <- as.name("data")
        environment()
    })
    parent.env(env) <- tmb_envir
    environment(oneStepPredict_patch) <- env
    oneStepPredict_patch(obj,
                         observation.name=observation.name,
                         data.term.indicator=data.term.indicator,
                         ...)
}

OBS_ENV <- reporter()

##' @describeIn TMB-interface Mark the observation
##' Mark observation to be used by either \code{oneStepPredict} or by \code{obj$simulate}.
##' If your objective function is using an observation \code{x}, you simply need
##' to run \code{x <- OBS(x)} \emph{inside the objective function}.
##' This will (1) allow \code{oneStepPredict} to change the class of \code{x} to
##' \code{"osa"} (\link{OSA-residuals}) or (2) allow \code{obj$simulate} to change the class of \code{x} to
##' \code{"simref"} (\link{Simulation}) on request.
##' @param x Observation object
OBS <- function(x) {
    ## Four evaluation modes (!)
    ## 1. MakeADFun for the first time (ad_context()==TRUE and typeof(obs)=="numeric")
    ## 2. MakeADFun for OSA (ad_context()==TRUE and OBS_ENV contains 'obs' of class=="osa")
    ## 3. Normal 'double' evaluation mode (ad_context()==FALSE and typeof(obs)="numeric")
    ## 4. Simulation mode (ad_context()==FALSE and OBS_ENV contains 'obs' of class=="simref")
    nm <- deparse(substitute(x))
    xobs <- OBS_ENV$get(nm)
    if (is.null(xobs)) {
        ## Are we running MakeADFun for the first time (i.e. taping)?
        ## ==> Store the observation inside OBS_ENV
        if (ad_context()) {
            OBS_ENV$set(nm, x)
        }
        return (x)
    }
    xobs
}

##' @describeIn OSA-residuals Subset observations marked for OSA calculation.
##' This function makes sure that when you subset an observation of class \code{"osa"} such as
##' \code{obs <- new("osa", x=advector(matrix(1:10,2)), keep = cbind(rep(TRUE,10),FALSE,FALSE))}
##' the 'keep' attribute will be adjusted accordingly
##' \code{obs[,1:2]}
##' @param x Object of class 'osa'
"[.osa" <- function(x, ...) {
    keep <- x@keep
    ord <- attr(keep, "ord")
    x <- x@x
    x. <- structure(NextMethod(), class=class(x))
    ## Should be fast due to 'ALTREP' (?)
    ind <- seq_along(x); dim(ind) <- dim(x)
    x <- ind; ind <- NextMethod()
    ## Remaining subsets uses index vector
    keep <- keep[ind, , drop=FALSE]
    if (!is.null(ord)) {
        attr(keep, "ord") <- ord[ind]
    }
    new("osa", x=x., keep=keep)
}

##' @describeIn OSA-residuals Equivalent of \link[base]{length}
"length.osa" <- function(x) length(x@x)
##' @describeIn OSA-residuals Equivalent of \link[base]{dim}
"dim.osa" <- function(x) dim(x@x)
##' @describeIn OSA-residuals Equivalent of \link[base]{is.array}
is.array.osa <- function(x) is.array(x@x)
##' @describeIn OSA-residuals Equivalent of \link[base]{is.matrix}
is.matrix.osa <- function(x) is.matrix(x@x)

dGenericOSA <- function(.Generic, x, ..., log) {
    if (!log) stop("'OSA' is for *log* density evaluation only")
    keep <- x@keep
    x <- x@x
    dfun <- match.fun(.Generic)
    ans <- dfun(x, ..., log=TRUE) * keep[,1]
    if (ncol(keep) == 3) { ## CDF method
        substring(.Generic, 1, 1) <- "p"
        pfun <- match.fun(.Generic)
        F <- pfun(x, ...) ## log=FALSE lower.tail=TRUE (default)
        ## NOTE: Very innocent fudge factor. It holds that
        ##   F + 1e-300 == F
        ## for all F in [1e-183, 1]
        ans <- ans + log(F + 1e-300) * keep[,2]   ## lower
        ans <- ans + log((1-F) + 1e-300) * keep[,3] ## upper
    }
    if (log) ans else exp(ans)
}
