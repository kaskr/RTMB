library(RTMB)

set.seed(1)
library(MASS)

simdata <- function(){
  local({
    rho=0.9
    sds=seq(0.5,2,length=stateDim)
    sdObs=rep(1,stateDim);
    corrMat=matrix(0.0,stateDim,stateDim)
    for(i in 1:stateDim){
      for(j in 1:stateDim){
        corrMat[i,j] = rho^abs(i-j)
      }
    }
    Sigma=corrMat*(sds %o% sds)
    d=matrix(NA,timeSteps,stateDim)
    obs=d;
    ##init state
    d[1,] = rnorm(stateDim);
    i=1;
    obs[i,] = d[i,] + rnorm(stateDim,rep(0,stateDim),sdObs)
    for(i in 2:timeSteps){
      d[i,] = d[i-1,] + mvrnorm(1,rep(0,stateDim),Sigma=Sigma)
      obs[i,] = d[i,] + rnorm(stateDim,rep(0,stateDim),sdObs)
    }
    matplot(d,type="l")
    matpoints(obs);
  },.GlobalEnv)
}
stateDim=3
timeSteps=100
simdata()
data <- list(obs=t(obs))
parameters <- list(
  u=data$obs*0,
  transf_rho=0.1,
  logsds=sds*0,
  logsdObs=sdObs*0
)

## =============== RTMB objective function

### Parameter transform
trf <- function(x) {2/(1 + exp(-2 * x)) - 1;}
### Objective
f <- function(
              transf_rho,
              logsds,
              logsdObs,
              u) {
    rho <- trf(transf_rho)
    sds <- exp(logsds)
    sdObs <- exp(logsdObs)
    cov <- outer(1:stateDim,
                 1:stateDim,
                 function(i,j) rho^(abs(i-j)) * sds[i] * sds[j])
    ADREPORT(cov)
    du <- diff(t(u))
    ans <- 0
    ans <- ans - sum(dmvnorm(du, 0, cov, log=TRUE))
    ans <- ans - sum(dnorm(data$obs, u, sdObs, log=TRUE))
    ans
}

## Test eval:
do.call("f", parameters)

obj <- MakeADFun(function(p)do.call("f",p), parameters, random="u")
TMB::newtonOption(obj, smartsearch=FALSE)

obj$fn()
obj$gr()
system.time(opt <- do.call("optim",obj))
pl <- obj$env$parList() ## <-- List of predicted random effects
matpoints(t(pl$u),type="l",col="blue",lwd=2)
sdreport(obj)
