## Patched
oneStepPredict <- function(obj, ...) {
    oneStepPredict_patch <- TMB::oneStepPredict
    tmb_envir <- environment(oneStepPredict_patch)
    env <- local({
        MakeADFun <- RTMB::MakeADFun
        formals(MakeADFun)$data <- list() ## Dummy
        formals(MakeADFun)$func <- as.name("data")
        environment()
    })
    parent.env(env) <- tmb_envir
    environment(oneStepPredict_patch) <- env
    oneStepPredict_patch(obj, ...)
}
