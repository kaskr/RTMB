getObj <- function() {
    set.seed(123)
    data <- list(Y = rnorm(10) + 1:10, x=as.numeric(1:10), hej=function(x){print(x);exp(x)})
    attr(data, "check.passed") <- TRUE
    parameters <- list(a=0, b=0, logSigma=0)
    obj <- TMB::MakeADFun(data, parameters, DLL="CallRTMB")
    obj
}
