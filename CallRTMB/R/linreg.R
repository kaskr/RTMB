getObj <- function() {
    set.seed(123)
    data <- list(x = as.numeric(1:10),
                 Y = rnorm(10) + 1:10,
                 foo = function(x) { print(x); exp(x) } )
    ## Disable data sanity check (would reject function object)
    attr(data, "check.passed") <- TRUE
    parameters <- list(a=0, b=0, logSigma=0)
    obj <- TMB::MakeADFun(data, parameters, DLL="CallRTMB")
    obj
}
