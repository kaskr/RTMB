## Length based model

## Parameters
g <- 15                          ## Growth rate
z <- 1                           ## Total mortality
sdrate <- 4                      ## sd = sdrate * dt
meanw <- 2                       ## Mean of log weights
sdw <- 2                         ## SD of log weights
## Data
inputtimes <- c(-1,-2,-3,-4)     ## Recruitment times
time <- 0                        ## Now
sizegrid <- 1:100
agegrid <- time - inputtimes
## Simulate
sizeByAge <- function(age) {
    mu <- g * age
    sd <- sdrate * age
    exp(-z * age) * dnorm(sizegrid, mu, sd)
}
P <- sapply(agegrid, sizeByAge)

## Simulate hauls
set.seed(123)
nage <- ncol(P)
nhaul <- 100
W <- matrix(exp(rnorm(nhaul * nage, mean=meanw, sd=sdw)), nage, nhaul)
F <- P %*% W
N <- F; N[] <- rpois(length(F),lambda=F)

########################################################
## RTMB model

library(RTMB)

data <- list(sizegrid=sizegrid,
             agegrid=agegrid,
             N=N)

parameters <- list(g=g,
                   z=z,
                   logsdrate=log(sdrate),
                   meanw=meanw,
                   logsdw=log(sdw),
                   logW=log(W)*0
                   )

## Mark data terms (not yet exported)
## Note: 'Term' does *nothing* in normal evaluation mode
Term <- RTMB:::Term

pla <- function(sizegrid, age, g, z, sdrate) {
  mu <- g * age
  sd <- sdrate * age
  exp(-z*age) * dnorm(sizegrid,mu,sd)
}
func <- function(parms) {
    getAll(data, parms)
    N <- OBS(N)
    sdrate <- exp(logsdrate)
    sdw <- exp(logsdw)
    nage <- length(agegrid)
    nsize <- length(sizegrid)
    nhaul <- ncol(N)
    res <- 0
    res <- res - sum(dnorm(logW, meanw, sdw, TRUE))
    P <- AD(matrix(0, nsize, nage))
    for (i in 1:nage) {
        P[,i] <- pla(sizegrid, agegrid[i], g, z, sdrate)
    }
    W <- exp(logW)
    F <- P %*% W
    REPORT(P)
    REPORT(W)
    REPORT(F)
    res <- res - sum(Term(dpois(N, F, TRUE)))
    res
}

########################################################

## Fit model using plain Laplace Approximation
## obj <- MakeADFun(func,
##                  parameters=parameters,
##                  random="logW",
##                  silent=TRUE)
## start <- obj$par * 1.4
## obj$env$tracemgc <- TRUE
## system.time(opt <- nlminb(start,obj$fn,obj$gr)) ## elapsed 2.064

## ====>
## opt = list(par = c(g = 14.8552113710921, z = 0.871393836817725, logsdrate = 1.33823539748694,  meanw = 1.8246507955666, logsdw = 0.65392120940321), objective = 3172.08961587239,      convergence = 1L, iterations = 39L, evaluations = c("function" = 80L,      gradient = 40L), message = "false convergence (8)")

## Takes a while...
## prof <- TMB::tmbprofile(obj,
##                         1,
##                         adaptive=FALSE,
##                         h=5e-3,
##                         maxit=40,
##                         ytol=Inf)

## Fit model using *ridge corrected* Laplace Approximation
obj2 <- MakeADFun(func,
                  parameters=parameters,
                  random="logW",
                  silent=TRUE,
                  ridge.correct=TRUE) ## NOTE: Enables 'Term'wise correction
start <- obj2$par * 1.4
system.time(opt2 <- nlminb(start,obj2$fn,obj2$gr)) ## elapsed 1.556

opt2.expected = list(par = c(g = 14.9401312393329, z = 0.969063999683462, logsdrate = 1.3434745572236,  meanw = 2.08151774558094, logsdw = 0.573799886824663), objective = 3222.10860363538,      convergence = 0L, iterations = 25L, evaluations = c("function" = 33L,      gradient = 26L), message = "relative convergence (4)")

expect_true(all(abs(opt2$par - opt2.expected$par) < 1e-4))
expect_true(abs(opt2$objective - opt2.expected$objective) < 1e-4)
