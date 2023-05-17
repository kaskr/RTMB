library(RTMB)

## From package 'geoR':
matern <- function (u, phi, kappa) 
{
    if (is.vector(u)) 
        names(u) <- NULL
    if (is.matrix(u)) 
        dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf, 
        gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)), 
        1)
    if (!RTMB:::ad_context()) ## FIXME
        uphi[u > 600 * phi] <- 0
    return(uphi)
}

## Test:
n <- 100
set.seed(123)
x <- matrix(runif(2*n, 0, 10), n)
D <- as.matrix(dist(x))
C <- matern(D, phi=1.3, kappa=2.7)
x <- t(chol(C)) %*% rnorm(n)
data <- list(x=x, D=D)
parameters <- list(phi=.5, kappa=.5)
map <- NULL

################################################################################

func <- function(p) {
    C <- matern(D, p$phi, p$kappa)
    -dmvnorm(x, 0, C, log=TRUE)
}

################################################################################

obj <- MakeADFun(func, parameters, map=map, DLL="matern")
system.time( fit <- nlminb(obj$par, obj$fn, obj$gr, obj$he) )
system.time( rep <- sdreport(obj, fit$par, obj$he(fit$par)) )
print(rep)
