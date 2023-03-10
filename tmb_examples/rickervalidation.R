## rickervalidation
##
## Estimate and validate a Ricker model based on data
## simulated from the logistic map
##
## Compare Thygesen et al (submitted, 2016): Validation of state space models
## fitted as mixed effects models

library(RTMB)

## Simulate the logistic map 
R <- 3.8  # With this value, we're in the chaotic regime
f <- function(x) R*x*(1-x)

N <- 200
x <- numeric(N)
x[1] <- 0.1
for(i in 2:N) x[i] <- f(x[i-1])

## Initial guess on parameters for a Ricker model
K <- 1       # Carrying capacity
Q <- 0.1     # Process noise 
r <- 0.5     # Growth rate
theta <- 1   # With theta=1, the "theta logistic" model is the Ricker model
S <- 50      # Sample volume controlling measurement uncertainty

X=numeric(N)

## Generate data with Poisson noise
set.seed(1)
Y <- rpois(length(x),S*x)

## TMB Data
data <- list(Y=Y)

## Initial guess on parameters for the fitted Ricker model
parameters0 <- list(
  X=log((data$Y+1)/S),
  logr=log(r),
  logtheta=log(theta),
  logK=log(K),
  logQ=log(Q),
  logS=log(S)
  )

## Make TMB model 
fixed <- factor(NA)

fun <- function(parms) {
    Y <- OBS(Y) ## Mark Y for OSA calculation on request
    r = exp(parms$logr)
    theta = exp(parms$logtheta)
    K = exp(parms$logK)
    Q = exp(parms$logQ)
    S = exp(parms$logS)
    X <- parms$X ## Initially I forgot this line => big problems (X exists in .GlobalEnv!)
    Xprev <- head(X, -1)
    Xnext <- tail(X, -1)
    dX <- Xnext - Xprev
    dXpred <- r * ( 1.0 - (exp(Xprev) / K)^theta )
    nll <- -sum(dnorm(dX, dXpred, sqrt(Q), TRUE))
    nll <- nll - sum(dpois(Y, S * exp(X), TRUE))
    nll
}

obj <- MakeADFun(fun,
                  parameters0,
                  random="X",
                  map=list(logtheta=fixed,logS=fixed))

## Fit the model 
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))    
sd <- sdreport(obj)

pred  <- oneStepPredict(obj,
                        method="oneStepGeneric",discrete=TRUE,range=c(0,Inf),
                        conditional=1, ## Skip first residual
                        parallel=FALSE)

layout(t(1:3))
plot(Y,xlab="Time",ylab="Y")
## Plot residual against previous observation:
res <- pred$res
obs <- head(Y, -1)
plot(obs, res,
     xlab=expression(y[i]),ylab=expression(u[i+1]))
abline(0,0)
summary(lm(res ~ obs))
## Plot data generating and fitted model
## The fitted model in natural abundance (not log)
g <- function(n,r,K) n*exp(r*(1-n/K))
nv <- seq(0,1,0.01)
## Median and +/- 1 sd in log domain for the fitted model
gn <- g(n=nv,r=exp(opt$par["logr"]),K=exp(opt$par["logK"]))
qu <- gn * exp(sqrt(exp(opt$par["logQ"])))
ql <- gn / exp(sqrt(exp(opt$par["logQ"])))
plot(c(0,1),c(0,max(qu)),type="n",xlab=expression(N[t]),ylab=expression(N[t+1]))
## Confidence regions
polygon(c(nv,rev(nv)),c(qu,rev(ql)),col="grey")
## Add median
lines(nv,gn,lwd=2)
## Add data generating logistic map
lines(nv,f(nv),lty="dashed",lwd=3)
