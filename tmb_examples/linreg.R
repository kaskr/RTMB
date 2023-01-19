library(RTMB)
set.seed(123)
data <- list(Y = rnorm(10) + 1:10, x=1:10)
parameters <- list(a=0, b=0, logSigma=0)

f <- function(parameters) {
    Y <- data$Y
    x <- data$x
    a <- parameters$a
    b <- parameters$b
    logSigma <- parameters$logSigma
    ## FIXME:
    ## ADREPORT(exp(2*logSigma));
    nll = -sum(dnorm(Y, a+b*x, exp(logSigma), TRUE))
    nll
}
obj <- MakeADFun(f, parameters)

obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)
