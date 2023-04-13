## Multivatiate SV model from Table 5 of Skaug and Yu "A flexible and automated likelihood based 
## framework for inference in stochastic volatility models." Computational Statistics & Data Analysis 76 (2014): 642-654.
## Negative joint likelihood (nll) of data and parameters
f <- function(parms) {
    "[<-" <- ADoverload("[<-")
    ## Optionally mark the observation object (so can 'checkConsistency' and 'oneStepPredict')
    X <- OBS(X)
    ## Random effects
    h <- parms$h
    ## Parameters on natural scale and remove "parms"
    sigma <- exp(parms$log_sigma)
    phi <- parms$phi
    sigma_init <- sigma/sqrt(1-phi^2)  
    nll <- 0  # Start collecting contributions
    ## Process likelihood
    for (j in 1:p) {
        nll <- nll - dautoreg(h[,j], phi=phi[j], scale=sigma_init[j], log=TRUE)
    }
    ## Parameterizes correlation matrix of X in terms of Cholesky factor
    L <- diag(p)
    L[lower.tri(L)] <- parms$off_diag_x   
    R <- cov2cor(L%*%t(L))  # Correlation matrix of X (guarantied positive definite)
    ## Scaling parameter
    mux <- matrix(parms$mu_x, nrow(h), ncol(h), byrow=TRUE)
    sigma_y <- exp( .5 * (h + mux) )
    ## Observation likelihood (dmvnorm is row-wise vectorized)
    nll <- nll - sum(dmvnorm(X, 0, R, scale=sigma_y, log=TRUE))
    nll
}

##' Fit sdvmulti
##' @param X Data
##' @examples
##' sdvmulti()
sdvmulti <- function(X=y) {
    p <- ncol(X)
    ## Parameter initial guess
    parameters <- list(
        phi        = rep(0.97,p),               # See eqn (12) in Skaug and Yu (2014)
        log_sigma  = rep(-1.7,p),               #       ---------||---------
        mu_x       = rep(-0.5,p),               #       ---------||---------
        off_diag_x = rep(0.0,p),                #       ---------||---------
        h          = matrix(0.0,nrow=n,ncol=p)  #       ---------||---------
    )
    data <- local({X <- X; environment()})
    environment(f) <- data
    obj <- MakeADFun(f, parameters, random="h")
    opt <- nlminb(obj$par, obj$fn, obj$gr,
           lower = c(rep(-.99,3), rep(-3.0,3), rep(-3.0,3), rep(-5.0,3)),
           upper = c(rep( .99,3), rep( 3.0,3), rep( 3.0,3), rep( 5.0,3)) )
    rep <- sdreport(obj)
    environment()
}
