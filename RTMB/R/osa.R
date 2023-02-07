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
setMethod("dpois", "osa", function(x, lambda, log) {
    keep <- x@keep
    x <- x@x
    ans <- dpois(x, lambda, log=TRUE) * keep
    if (log) ans else exp(ans)
})
