source("sam_data.R")

library(RTMB)

parameters <- list(
  logFpar      = c(-11, -10, -10, -9, -9),
  logQpow      = numeric(0),
  logSdLogFsta = rep(-0.693147,max(data$keyVarF)+1),
  logSdLogN    = c(0.356675, -0.356675),
  logSdLogObs  = c(-0.356675, -0.356675, -0.356675, -0.356675, -0.356675),
  rec_loga     = 1,
  rec_logb     = -12,
  rho          = 0.5,
  logScale     = numeric(data$noScaledYears),
  logScaleSSB  = if(any(data$fleetTypes %in% c(3,4))) {numeric(1)} else {numeric(0)},
  logPowSSB    = if(any(data$fleetTypes == 4))        {numeric(1)} else {numeric(0)},
  logSdSSB     = if(any(data$fleetTypes %in% c(3,4))) {numeric(1)} else {numeric(0)},
  U = matrix(0, nrow=max(data$keyLogFsta)+1 + data$maxAge-data$minAge+1, ncol=data$noYears)
)
data$nlogF = max(data$keyLogFsta)+1
data$nlogN = data$maxAge-data$minAge+1

if (FALSE) {
    attach(data)       ## discouraged but
    attach(parameters) ## VERY handy while developing!
}

### Parameter transform
f <- function(x) {2/(1 + exp(-2 * x)) - 1}
square <- function(x){x*x}
func <- function(
                 logFpar,
                 logQpow,
                 logSdLogFsta,
                 logSdLogN,
                 logSdLogObs,
                 rec_loga,
                 rec_logb,
                 rho,
                 logScale,
                 logScaleSSB,
                 logPowSSB,
                 logSdSSB,
                 U) {
    logN <- U[  1:nlogN , , drop=FALSE]
    logF <- U[-(1:nlogN), , drop=FALSE]
    timeSteps <- ncol(logF)
    stateDimF <- nrow(logF)
    stateDimN <- nrow(logN)
    sdLogFsta = exp(logSdLogFsta)
    varLogN = exp(logSdLogN*2)
    varLogObs = exp(logSdLogObs*2)
    ## First take care of F
    fcor <- outer(1:stateDimF,
                  1:stateDimF,
                  function(i,j)(i!=j)*rho + (i==j))
    fsd <- sdLogFsta[keyVarF[1,]+1L]
    fvar <- outer(1:stateDimF,
                  1:stateDimF,
                  function(i,j)fcor[cbind(i,j)]*fsd[i]*fsd[j])
    ans <- 0
    ans <- ans - sum(dmvnorm( diff(t(logF)) , 0, fvar, log=TRUE))
    ## SIMULATE {
    ##   logF.col(i) = logF.col(i-1) + neg_log_densityF.simulate();
    ## }
    calcssb <- function(i)sum(exp(logN[,i]) * exp(-exp(logF [keyLogFsta[1,]+1L ,i] ) * propF[i,]-natMor[i,]*propM[i,])*propMat[i,]*stockMeanWeight[i,])
    ## FIXME: ssb <- sapply(1:timeSteps, calcssb)
    ssb <- do.call("c",lapply(1:timeSteps, calcssb))
    logssb <- log(ssb)
    ## Now take care of N
    nvar <- outer(1:stateDimN, 1:stateDimN,
                  function(i,j) (i==j)*varLogN[ keyVarLogN[1,i]+1L ])
    predN <- numeric(stateDimN)
    for(i in 2:timeSteps) {
        if(stockRecruitmentModelCode==0){ ## straight RW
            predN[1] = logN[1, i-1]
        } else {
            if (stockRecruitmentModelCode==1){ ##ricker
                predN[1] = rec_loga+log(ssb[i-1])-exp(rec_logb)*ssb[i-1]
            }else{
                if(stockRecruitmentModelCode==2){##BH
                    predN[1]=rec_loga+log(ssb[i-1])-log(1+exp(rec_logb)*ssb[i-1])
                }else{
                    stop("SR model code not recognized");
                }
            }
        }
        for(j in 2:stateDimN) {
            predN[j]=logN[j-1,i-1]-exp(logF[(keyLogFsta[1,j-1]+1L),i-1])-natMor[i-1,j-1]
        }
        if(maxAgePlusGroup==1){
            predN[stateDimN] = log(exp(logN[stateDimN-1,i-1]-exp(logF[(keyLogFsta[1,stateDimN-1]+1L),i-1])-natMor[i-1,stateDimN-1])+
                                   exp(logN[stateDimN,i-1]-exp(logF[(keyLogFsta[1,stateDimN]+1L),i-1])-natMor[i-1,stateDimN]))
        }
        ## SIMULATE {
        ##   logN.col(i) = predN + neg_log_densityN.simulate();
        ## }
        ans <- ans - dmvnorm(logN[,i], predN, nvar, log=TRUE) ## N-Process likelihood
    }
    ## Now finally match to observations
    minYear <- obs[1,1]
    predObs <- 0
    zz <- 0
    var <- 0
    for(i in 1:nobs){
        y <- obs[i,1] - minYear + 1 ## 1 based
        f <- obs[i,2] ## 1 based
        ft <- fleetTypes[f] ## 1 based
        a <- obs[i,3] - minAge + 1 ## 1 based
        zz <- exp(logF[keyLogFsta[1,a]+1L,y])+natMor[y,a]
        if(ft==0){## residual fleet
            predObs <- logN[a,y]-log(zz)+log(1-exp(-zz))
            if((keyLogFsta[f,a])>(-1)){
                predObs <- predObs + logF[keyLogFsta[1,a]+1L,y]
            }
        }else{
            if(ft==1){## comm fleet
                stop("Not implemented yet!!!")
            }else{
                if(ft==2){## survey
                    predObs=logN[a,y]-zz*sampleTimes[f]
                    if(keyQpow[f,a]>(-1)){
                        predObs = predObs*exp(logQpow[keyQpow[f,a]+1L])
                    }
                    if(keyLogFpar[f,a]>(-1)){
                        predObs <- predObs+logFpar[keyLogFpar[f,a]+1L]
                    }
                }else{
                    if(ft==3){## SSB survey -- nevermind for now
                        stop("Not implemented yet!!!")
                    }else{
                        if(ft==4){## SSB survey -- nevermind for now
                            stop("Not implemented yet!!!")
                        }
                    }
                }
            }
        }
        var <- varLogObs[keyVarObs[f,a]+1L]
        ans <- ans - dnorm(log(obs[i,4]),predObs,sqrt(var),log=TRUE)
        ## SIMULATE {
        ##   obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
        ## }
    }
    ans
}

## Test that we can evaluate using numeric types
environment(func) <- list2env(data)
do.call(func, parameters)

obj <- MakeADFun(function(p)do.call(func,p), parameters, random=c("U"), DLL="sam")
lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower["rho"] <- 0.01
upper["rho"] <- 0.99

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper))
rep <- sdreport(obj)
rep
