## Multivatiate SV model from Table 5 of Skaug and Yu "A flexible and automated likelihood based 
## framework for inference in stochastic volatility models." Computational Statistics & Data Analysis 76 (2014): 642-654.

library(RTMB)

## Read data
source("sdv_multi_data.R")
X=y

## Parameter initial guess
us <- unstructured(p) ## Unstructured pxp correlation
parameters <- list(
    phi        = rep(0.97,p),               # See eqn (12) in Skaug and Yu (2014)
    log_sigma  = rep(-1.7,p),               #       ---------||---------
    mu_x       = rep(-0.5,p),               #       ---------||---------
    off_diag_x = us$parms(),                #       ---------||---------
    h          = matrix(0.0,nrow=n,ncol=p)  #       ---------||---------
)

parameters$h <- t(parameters$h) ## Workaround: order as old TMB example for unittest

## Negative joint likelihood (nll) of data and parameters
f <- function(parms) {
    ## Optionally mark the observation object (so can 'checkConsistency' and 'oneStepPredict')
    X <- OBS(X)
    ## Make all parms visible
    getAll(parms)
    ## Parameters on natural scale
    sigma <- exp(log_sigma)
    sigma_init <- sigma/sqrt(1-phi^2)  
    h <- t(h)  # Workaround: order as old TMB example for unittest
    nll <- 0  # Start collecting contributions
    ## Process likelihood
    for (j in 1:p) {
        nll <- nll - dautoreg(h[,j], phi=phi[j], scale=sigma_init[j], log=TRUE)
    }
    ## Parameterizes correlation matrix of X in terms of Cholesky factor
    R <- us$corr(off_diag_x)
    ## Scaling parameter
    mux <- matrix(mu_x, nrow(h), ncol(h), byrow=TRUE)
    sigma_y <- exp( .5 * (h + mux) )
    ## Observation likelihood (dmvnorm is row-wise vectorized)
    nll <- nll - sum(dmvnorm(X, 0, R, scale=sigma_y, log=TRUE))
    nll
}

obj <- MakeADFun(f, parameters, random="h")
system.time(
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  lower = c(rep(-.99,3), rep(-3.0,3), rep(-3.0,3), rep(-5.0,3)),
                  upper = c(rep( .99,3), rep( 3.0,3), rep( 3.0,3), rep( 5.0,3)) )
)
rep <- sdreport(obj)
rep
