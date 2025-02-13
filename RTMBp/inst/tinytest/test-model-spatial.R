library(RTMBp)

## Read data
source("spatial_data.R")

## Euclidian distance matrix
dd <- as.matrix(dist(Z))

## Parameters
parameters <- list(
    b = c(0, 0),
    loga = 0,
    logsigma = 0,
    u = rep(0,n) )

## Objective
f <- function(params) {
    ## Get parameters
    getAll(params)
    ## Mark the response (optional)
    y <- OBS(y)
    ## Parameter transforms
    a <- exp(loga)
    sigma <- exp(logsigma)
    ## Covariance matrix of random effects
    cov <- exp( - a * dd )
    ## Draw random effects (see ?MVgauss)
    u %~% dmvnorm(mu=0, Sigma=cov)
    ## Latent log intensity
    eta <- X %*% b + sigma * u
    ## Draw data given random effects
    y %~% dpois(exp(eta))
    ## Add to sdreport
    ADREPORT(a)
    ADREPORT(sigma)
}

## Fit model
obj <- MakeADFun(f, parameters, random = "u", silent=TRUE)
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))

## Report output
rep <- sdreport(obj)

## Check
tab <- summary(rep, c("fixed", "report"))

expect_true(rep$pdHess, info="Hessian positive definite")
expect_true(max(abs(rep$gradient.fixed)) < 1e-6, info="Gradient close to zero")
expect_equal(tab, tab.expected, tol=1e-4, info="Parameter estimates and std errors")
