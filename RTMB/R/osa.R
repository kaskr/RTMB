## Patched
oneStepPredict <- function(obj,
                           observation.name=names(obj$env$osa)[1],
                           data.term.indicator="_RTMB_keep_",
                           ...) {
    ## FIXME: Check against valid OSA names
    ## names(obj$env$osa)
    obj$env$data[[observation.name]] <- obj$env$osa[[observation.name]]
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

OSA_ENV <- reporter()
OSA <- function(x) {
    if (!ad_context()) return (x)
    nm <- deparse(substitute(x))
    xosa <- OSA_ENV$get(nm)
    if (is.null(xosa)) {
        OSA_ENV$set(nm, x)
        xosa <- x
    }
    xosa
}

setClass("osa", list(x="ad", keep="ad"))

##' @DescribeIn OSA-residuals Subset observations marked for OSA calculation
##' @examples
##' x <- new("osa", x=advector(matrix(1:10,2)), keep = cbind(rep(TRUE,10),FALSE,FALSE))
##' x[,1:2]
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
        ans <- ans + log(F) * keep[,2]   ## lower
        ans <- ans + log(1-F) * keep[,3] ## upper
    }
    if (log) ans else exp(ans)
}

## FIXME: Autogenerate using 'distr.R'
setMethod("dpois", "osa", function(x, lambda, log) {
    dGenericOSA(.Generic, x, lambda, log=log)
})
setMethod("dnorm", "osa", function(x, mean, sd, log) {
    dGenericOSA(.Generic, x, mean, sd, log=log)
})
