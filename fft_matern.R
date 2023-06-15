library(RTMB)
set.seed(1)
log_dmvnorm_circ <- function(x, C) {
    sd <- sqrt(Re(fft(C)))
    x <- Mod(fft(x)) / sqrt(length(x))
    sum(dnorm(x, 0, sd, log=TRUE))
}

circDist <- function(dim) {
    cd <- function(n) pmin(1:n-1, n:1)
    g <- expand.grid(lapply(dim, cd))
    structure(sqrt(rowSums(g*g)), dim=dim)
}

matern <- function (u, phi, kappa) 
{
    uphi <- u/phi
    ifelse(u > 0, (((2^(-(kappa - 1)))/gamma(kappa)) * (uphi^kappa) * besselK(x = uphi, nu = kappa)), 1)
    ##return(uphi)
    ## dim(ans) <- dim(u)
    ## ans
}

## config:
dim <- c(200,200)
Dmax <- 10
## dim <- c(20,20)
## Dmax <- 2


##D <- circDist()
D <- circDist(dim)
C <- matern(D, 10, 2.5)
s <- Re(fft( fft(structure(rnorm(length(C)),dim=dim(C))) * sqrt(fft(C)) , inverse=TRUE)) / length(C)

## Eigen value approximation analysis (depends on smoothness)
lam <- fft(C)
##q <- quantile(Mod(lam), .99)
##lam[Mod(lam)<q] <- 0
lam[D>Dmax] <- 0
Chat <- Re(fft(lam, inverse=TRUE)) / length(C)
range(C-Chat)
plot(as.vector(C),as.vector(Chat))
abline(0,1)
sum(lam!=0)

## Observations
loc <- sample(1:length(C), 10000)
obs <- s[loc] + rnorm(length(loc), sd=1)
## objective function
parms <- list(mu=0, phi=10, kappa=2.5, sd=1, sd_obs=1, u=numeric(sum(D<=Dmax)))
f <- function(parms) {
    getAll(as.list(parms))
    C <- matern(D, phi, kappa)
    i <- order(D)[1:length(u)]
    ans <- -sum(dnorm(u, log=TRUE))
    U <- C * 0
    U[i] <- u
    ##U <- adcomplex(U)
    ## simulation
    S <- Re(fft( fft(U) * Re(sqrt(fft(C))) , inverse=TRUE)) / length(C)
    REPORT(S)
    ans <- ans - sum(dnorm(obs, sd * S[loc] + mu, sd=sd_obs, log=TRUE))
    ans
}
TMB::config(tmbad.sparse_hessian_compress=1)
obj <- MakeADFun(f,parms,random="u")
obj$fn()
obj$gr()
opt <- nlminb(obj$par,obj$fn,obj$gr)

lpb <- obj$env$last.par.best
qw <- obj$report(lpb)
