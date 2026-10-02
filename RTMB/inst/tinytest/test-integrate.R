## Testing RTMB::integrate

library(tinytest)
library(RTMB)
tol <- .Machine$double.eps^.25 ## Default used by integrate

## Gauss-weibull
f <- function(p) {
  getAll(p)
  integrand <- function(u) exp( dweibull(u,shape,scale,log=TRUE) + dnorm(x,u,sd,log=TRUE) )
  integrate(integrand, 0, Inf)$val
}
F <- MakeTape(f, list(shape=1,scale=1,sd=1,x=1))
p <- list(shape=1,scale=1,sd=1,x=7.982944)
expect_true(abs(f(p)-F(p)) < tol, info="integrate 0 to Inf")
p <- list(shape=1.1,scale=1.1,sd=1.1,x=5.779196)
expect_true(abs(f(p)-F(p)) < tol, info="integrate 0 to Inf")

## Spiky dnorm
f <- function(x, sd) dnorm(x, sd=1e-3)
C <- function(sd) integrate(dnorm, -Inf, Inf, sd=sd, rel.tol=1e-8, abs.tol=0)$value
F <- MakeTape(C, 1)
one <- sapply(10^(-3:3),F)
expect_true( all( abs(one - 1) < tol ) , info="integrate spiky dnorm")

## Singular beta
s <- 0.000001
f <- function(x) dbeta(x, s, s)
F <- MakeTape(function(x)integrate(f, 0, x, subdivisions=1e5, abs.tol=0)$val, 0)
expect_true( abs(F(.5) - .5) < tol , info="integrate singular dbeta")

## stats::integrate can't do this
if (FALSE) {
  integrate(f,0,.5) ## Wrong
  integrate(f,0,.5,subdivisions=1e5, abs.tol=0) ## Error
}

## Subdivision test
if (FALSE) {
  f <- function(x) sin(exp(x))
  F <- MakeTape(function(x) {qw<-integrate(f, 0, x,subdivisions = 1e5); do.call("c",qw)}, 0)
  F(12) ## Uses 11234 subdivisions
  integrate(f, 0, 12, subdivisions=1e5) ## Uses 11063 subdivisions
}

## Tolerance test
if (FALSE) {
  f <- function(x) sin(exp(x))*exp(x)
  Ftrue <- function(x) -cos(exp(x)) + cos(1)
  F <- MakeTape(function(x) {qw<-integrate(f, 0, x,subdivisions = 1e5); do.call("c",qw)}, 0)
  F(6)[1] - Ftrue(6) ## OK
  F(7)[1] - Ftrue(7) ## OK
  F(12)[1] - Ftrue(12) ## oops not within tolerance!
}

if (FALSE) {
  f <- function(x) sin(exp(3*sin(3*x))) + 1
  F <- MakeTape(function(x) {qw<-integrate(f, 0, x,subdivisions = 1e5); do.call("c",qw)}, 0)
  F(10)[1] ## 25 subdiv
  integrate(f,0,10) ## 31 subdiv
}

## stats::integrate can't do this
f <- function(x) integrate(dgamma, 0, Inf, shape=x, subdivisions=1e5, rel.tol=1e-10, abs.tol=0)$value
F <- MakeTape(f, 1)
x <- 10^-(1:5)
y <- sapply(x, F)
expect_true(max(abs(y-1)) < 1e-8) ## observed tol=1.184941e-11

f <- function(t) dnorm(t, 0, 0.1)
b <- 0.01
F <- MakeTape(function(x) integrate(f, -3, x, rel.tol = 1e-8, abs.tol = 0)$value, b)
ans <- F$jacobian(b)
expect_true(max(abs(ans-f(b))) < 1e-10) ## observed tol=1.332268e-15

## stats::integrate can't do this
f <- function(x) integrate(dt, -Inf, Inf, df=x, rel.tol=1e-10,abs.tol=0)$value
F <- MakeTape(f, 1)
x <- 10^-(0:5)
y <- sapply(x, F)
expect_true(max(abs(y-1)) < 1e-5) ## observed tol=2.517603e-07

f <- function(x) integrate(df, 0, x, df1=x, df2=x, abs.tol=0, rel.tol=1e-10)$val
F <- MakeTape(f, 1)
x <- 2.025203
expect_true( abs(F(x) - pf(x,x,x)) < 1e-12 ) ## observed tol=1.887379e-15

if (at_home()) {
  ## From glmmAdaptive:
  set.seed(1234)
  n <- 300 # number of subjects
  K <- 3*4 # number of measurements per subject
  t_max <- 15 # maximum follow-up time
  ## we construct a data frame with the design:
  ## everyone has a baseline measurement, and then measurements at K time points
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   time = gl(K, 1, n*K, labels = paste0("Time", 1:K)),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))
  ## design matrices for the fixed and random effects
  X <- model.matrix(~ sex * time, data = DF)
  Z <- model.matrix(~ 1, data = DF)
  betas <- c(-2.13, 1, rep(c(1.2, -1.2), K-1)) # fixed effects coefficients
  D11 <- 1 # variance of random intercepts
  ## we simulate random effects
  b <- cbind(rnorm(n, sd = sqrt(D11)))
  ## linear predictor
  eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
  ## we simulate binary longitudinal data
  DF$y <- rbinom(n * K, 1, plogis(eta_y))
  ## Density of Gauss-binomial mixture
  GaussBinomial <- function(x, mu, sd) {
    integrand <- Vectorize(function(u) {
      loglik <- dnorm(u, log=TRUE)
      logitp <- sd * u + mu
      loglik <- loglik + sum(dbinom_robust(x, 1, logitp, log=TRUE))
      exp(loglik)
    })
    integrate(integrand, -Inf, Inf, rel.tol=1e-8, abs.tol=0)$value
  }
  ## Use 'Vectorize' to speedup MakeADFun:
  GaussBinomial <- Vectorize(GaussBinomial)
  X <- model.matrix(~sex + time,data=DF)
  func <- function(p) {
    getAll(p, DF)
    mu <- X %*% beta
    sd <- exp(logsd)
    y <- split(y, id)
    mu <- split(mu, id)
    -sum(log(GaussBinomial(y, mu, sd)))
  }
  parameters <- list(beta = rep(0,ncol(X)), logsd=0)
  ## func(parameters) ## 2052.784
  obj <- MakeADFun(func, parameters, silent=TRUE)
  expect_true(abs(obj$fn() - 2052.78381412415) < 1e-6, info="glmm adaptive start value")
  expect_true(all(is.finite(obj$gr())), info="glmm adaptive finite gradient")
  system.time( fit <- nlminb(obj$par, obj$fn, obj$gr) )
  expect_true( abs(fit$objective - 1814.28886426115) < 1e-6, info="glmm adaptive final value")
}
