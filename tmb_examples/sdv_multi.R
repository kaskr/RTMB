# Multivatiate SV model from Table 5 of Skaug and Yu "A flexible and automated likelihood based 
# framework for inference in stochastic volatility models." Computational Statistics & Data Analysis 76 (2014): 642-654.

library(RTMB)

# Read data
source("sdv_multi_data.R")
X=y

## Parameter initial guess
parameters <- list(
    phi        = rep(0.97,p),               # See eqn (12) in Skaug and Yu (2014)
    log_sigma  = rep(-1.7,p),               #       ---------||---------
    mu_x       = rep(-0.5,p),               #       ---------||---------
    off_diag_x = rep(0.0,p),                #       ---------||---------
    h          = matrix(0.0,nrow=n,ncol=p)  #       ---------||---------
)

# Negative joint likelihood (nll) of data and parameters
f <- function(parms) {
  
  # Parameters on natural scale and remove "parms"
  sigma <- exp(parms$log_sigma)
  phi <- parms$phi
  sigma_init <- sigma/sqrt(1-phi^2)  
  h <- parms$h

  nll = 0  # Start collecting contributions
  
  # Prior on state variable (y)
  nll <- nll - sum(dnorm(h[1,],0,sigma_init,log=TRUE))   # Start in stationary distribution
  for(j in 1:p) 
    nll <- nll - sum(dnorm(h[-1,j],phi[j]*h[-n,j],sigma[j],log=TRUE))  # AR(1) process

  # Parameterizes correlation matrix of X in terms of Cholesky factor
  L <- diag(p)
  L[lower.tri(L)] <- parms$off_diag_x   
  row_norms <- apply(L, 1, function(row) sqrt(sum(row^2)))
  L <- t(t(L)/row_norms)
  R <- L%*%t(L)  # Correlation matrix of X (guarantied positive definite)
  
  # Likelihood of data (X) given h
  for(i in 1:n){
    sigma_y <- exp(0.5*(parms$mu_x + h[i,]));
    Cov <- diag(sigma_y) %*% R %*% diag(sigma_y)
    nll <- nll - sum(dmvnorm(X[i,], 0, Cov, log=TRUE))
  }
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

