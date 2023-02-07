## randomwalkvalidation
##
## Estimate and validate a random walk model with and without drift
##
## Compare Thygesen et al (submitted, 2016): Validation of state space models
## fitted as mixed effects models
##
## Uffe Høgsbro Thygesen and Kasper Kristensen, 2016

library(RTMB)

## For reproducible results
set.seed(123)

## Simulate data with these parameters
mu <- 0.75
sigma <- 1
s <- 1
huge <- 1e3

## Simulate random track
nT <- 100
X <- c(0,cumsum(rnorm(nT-1,mean=mu,sd=sigma)))

## Simulate random measurements
Y <- X + rnorm(nT,sd=s)

data <- list(y=Y,huge=huge)
parameters <- list(x=X,mu=0,logsigma=log(sigma),logs=log(s))

## objective
func <- function(x, mu, logsigma, logs) {
    ## Enable OSA on request
    y <- OSA(data$y)
    ## Initial condition
    nll <- -dnorm(x[1], sd=huge, log=TRUE)
    ## Increments
    dx <- diff(x)
    nll <- nll - sum(dnorm(dx, mean=mu, sd=exp(logsigma), log=TRUE))
    ## Observations
    nll <- nll - sum(dnorm(y, x, exp(logs), TRUE))
    nll
}
f <- function(p)do.call("func", p)

## Estimate states and parameters under H0: mu=0
obj0 <- MakeADFun(f,parameters,random=c("x"),map=list(mu=factor(NA)))
opt0 <- do.call("optim",obj0)
sdr0 <- sdreport(obj0)
estX0 <- summary(sdr0,"random")

## Estimate states and parameters under H1: mu != 0
obj1 <- MakeADFun(f,parameters,random=c("x"))
opt1 <- do.call("optim",obj1)
sdr1 <- sdreport(obj1)
estX1 <- summary(sdr1,"random")

## ### One-step predictions

## ## Generate one step predictions with the models fitted under H0 and H1
predict0  <- oneStepPredict(obj0, method="fullGaussian")
predict1  <- oneStepPredict(obj1, method="fullGaussian")

## ## Test if the prediction errors are unbiased
print(TestPred1 <- anova(lm(predict1$residual ~ 0),lm(predict1$residual ~ 1)))
print(TestPred0 <- anova(lm(predict0$residual ~ 0),lm(predict0$residual ~ 1)))
