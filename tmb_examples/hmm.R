library(RTMB)
TapeConfig(vectorize="enable") ## Optional (speeds up this model)
library(Matrix) ## expm

## Simulate SDE
lambda <- 1
gamma <- 1
sigmaX <- 0.5
sigmaY <- 0.1
x0 <- 1
par.true <- c(x0     = x0,
              lambda = lambda,
              gamma  = gamma,
              logsX  = log(sigmaX),
              logsY  = log(sigmaY))
f <- function(x) lambda * x - gamma * x^3
g <- function(x) x * 0 + sigmaX
Tsim <- 0.01; T <- 50
set.seed(1)
euler <- function(x0,f,g,tvec,dB=NULL){
    X <- numeric(length(tvec))
    X[1] <- x0
    dt <- diff(tvec)
    if(is.null(dB)) dB <- rnorm(length(dt),sd=sqrt(dt))
    for(i in 1:(length(tvec)-1))
        X[i+1] <- X[i] + f(X[i])*dt[i] + g(X[i])*dB[i]
    return(X)
}
tsim <- seq(0,T,Tsim)
Xsim <- euler(x0,f,g,tsim)
# Measure every 10th simulated value
iobs <- seq(1,length(tsim),10)
## Measurements
Y <- rnorm(length(iobs), mean = Xsim[iobs], sd = sigmaY)
plot(tsim,Xsim,type="l")
points(tsim[iobs],Y,col="red")
grid <- seq(-3,3,length=101)

data <- list(
    grid = grid,
    dt = diff(tsim[iobs])[1],
    yobs = findInterval(Y, grid) ##-1
)
parameters <- list(
    lambda=0, gamma=0, logsX=0, logsY=0
)

## ================ Translation of original TMB example

## Get sde_t
sde <- function(lambda, gamma, sigmaX) {
    list(advection  = function (x) lambda * x - gamma * x^3,
         dispersion = function (x) x*0 + sigmaX )
}
## Finite volume discritize advection diffusion equation. Assuming
## equidistant grid and using central difference scheme for advection.
fvade <- function(sde, grid) {
    h <- diff(grid)[1]
    D <- .5 * sde$dispersion(grid)^2
    ## Helper
    ## Simpler (FIXME): rbind(0, cbind(diag(x), 0))
    subdiag <- function(x) {
        n <- length(x) + 1
        m <- matrix(0, n, n)
        m[col(m) == row(m) - 1] <- x
        m
    }
    ## L matrix (dim nvol = n-1)
    L <- subdiag(D[-c(1, length(D))])
    L <- L+t(L)
    diag(L) <- -colSums(L)
    L <- L / (h*h)
    ## Setup advection term
    xm <- .5 * (head(grid, -1) + tail(grid, -1))
    v <- sde$advection(xm)
    G <- -.5 * subdiag(head(v, -1)) + .5 * t(subdiag(tail(v, -1)))
    diag(G) <- -colSums(G)
    G <- G / h
    ## Generator
    A <- -G + L
    environment()
}
## hmm filter
hmm.filter <- function(A, grid, dt) {
    P <- expm(A*dt)
    nvol <- nrow(A)
    h <- diff(grid)[1]
    P0 <- diag(nvol)
    setGaussianError <- function(sd) {
        xm <- .5 * (head(grid, -1) + tail(grid, -1))
        P0 <<- dnorm(outer(xm, xm, "-"), sd=sd)
        P0 <<- P0 / colSums(P0)[col(P0)]
        NULL
    }
    ## sum(px)=1 obs=integer
    px <- NULL
    update <- function(yobs) {
        Ppx <- P %*% px
        Jslice <- Ppx * P0[yobs, ]
        Ly <- sum(Jslice)
        px <<- Jslice / Ly
        log(Ly)
    }
    ## log likelihood
    loglik <- function(obs) {
        px <<- rep(1, nvol)
        px <<- px / sum(px) ## prob
        px <<- px / h ## density
        ans <- 0
        for (o in obs)
            ans <- ans + update(o)
        ans
    }
    environment()
}

## FIXME:
colSums <- function(x) apply(x, 2, sum)
## Helper to make all objects visible inside objective:
## (Is this generally useful...?)
getAll <- function(x) {
    fr <- parent.frame()
    for (nm in names(x))
        fr[[nm]] <- x[[nm]]
    invisible(NULL)
}

func <- function(parameters) {
    getAll(parameters)
    getAll(data)
    sigmaX <- exp(logsX)
    sigmaY <- exp(logsY)
    ## Construct the SDE
    sde_object <- sde(lambda, gamma, sigmaX)
    ## Finite volume disretize and report generator
    fvol <- fvade(sde_object, grid)
    REPORT(fvol$A)
    ## Construct likelihood function
    hmm <- hmm.filter(fvol$A, fvol$grid, dt)
    hmm$setGaussianError(sigmaY)
    -hmm$loglik(yobs)
}

func(parameters)

obj <- MakeADFun(func, parameters=parameters)
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
(sdr <- sdreport(obj, opt$par))

dp <- par.true[-1] - opt$par
H <- solve(sdr$cov.fixed)
## Both less than qchisq(.95,df=4):
t(dp) %*% H %*% dp
2 * (obj$fn(par.true[-1]) - obj$fn(opt$par))
