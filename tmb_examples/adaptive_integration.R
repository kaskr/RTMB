## Exact Likelihood inference using adaptive integration
require(RTMB)

## Simulate data
set.seed(123)
ndat <- 1000
c1 <- rnorm(ndat, sd=1)
c2 <- rnorm(ndat, sd=1)
ngroup <- 11
c3 <- cut(runif(ndat), ngroup)
eps <- rnorm(ndat, sd=.5)
p <- plogis(c1 + c2 + seq(-1,1,length=ngroup)[c3] + eps)
n <- rep(10, ndat)
x <- rbinom(ndat, size=n, prob=p)
A <- model.matrix(~ c1 + c2 + c3 - 1)

## Data and parameters for TMB
data <- list(x=x, n=n, A=A)
parameters <- list(b = rep(0,ncol(A)), logsd=0)

## Density of Gauss-binomial mixture
GaussBinomial <- function(x, n, mu, sd) {
  integrand <- function(u) {
    loglik <- dnorm(u, log=TRUE)
    logitp <- sd * u + mu
    loglik <- loglik + dbinom_robust(x, n, logitp, log=TRUE)
    exp(loglik)
  }
  integrate(integrand, -Inf, Inf)$value
}
## Helper to avoid taping the same function many times.
## TODO: Move to ?RTMB::Tape
TapeOnFirstEval <- function(f) {
  if (!RTMB:::ad_context()) return(f)
  f <- match.fun(f)
  firstEval <- TRUE
  firstArgs <- NULL
  firstTape <- NULL
  function(...) {
    if (firstEval) {
      firstArgs <<- as.relistable(list(...))
      x0 <- RTMB:::getValues(do.call("c", firstArgs))
      vf <- function(x) do.call("f", relist(x, firstArgs))
      firstTape <<- MakeTape(vf, x0)
      firstEval <<- FALSE
    }
    firstTape(c(...))
  }
}

func <- function(p) {
  getAll(p, data)
  mu <- A %*% b
  sd <- exp(logsd)
  ans <- 0
  ## This line is to speedup MakeADFun:
  GaussBinomial <- TapeOnFirstEval(GaussBinomial)
  for(i in 1:length(x)) {
    ans <- ans - log( GaussBinomial(x[i], n[i], mu[i], sd) );
  }
  ans
}

## Fit model
model <- MakeADFun(func, parameters)
model$fn()
model$gr()
system.time( fit <- nlminb(model$par, model$fn, model$gr) )
rep <- sdreport(model)
print(rep)
