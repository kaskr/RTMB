## Testing RTMB::integrate

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
