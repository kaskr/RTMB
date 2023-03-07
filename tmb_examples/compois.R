library(RTMB)

set.seed(123)
nu <- .1
mode <- 10
domain <- 0:100
prob <- dpois(domain, lambda=mode)^nu; prob <- prob / sum(prob)
sum(prob * domain) ## mean
x <- sample(domain, size=1e4, replace=TRUE, prob = prob)

## TMB data
data <- list( x = x )
parameters <- list( logmu = 0, lognu = 0 )

## Parameterization through the mode
parameterization <- "mode"
func <- function(p) {
    mu <- exp(p$logmu)
    nu <- exp(p$lognu)
    ADREPORT(mu)
    ADREPORT(nu)
    if (parameterization == "mode")
        -sum(dcompois(x, mu, nu, TRUE))
    else ## mean
        -sum(dcompois2(x, mu, nu, TRUE))        
}
obj <- MakeADFun(func, parameters)
system.time(fit.mode <- nlminb(obj$par, obj$fn, obj$gr, obj$he))
rep.mode <- sdreport(obj)

## Parameterization through the mean
parameterization <- "mean"
obj <- MakeADFun(func, parameters)
system.time(fit.mean <- nlminb(obj$par, obj$fn, obj$gr, obj$he))
rep.mean <- sdreport(obj)

summary(rep.mode, "report")
summary(rep.mean, "report")
