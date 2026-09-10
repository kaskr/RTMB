## Test sequential reduction
library(tinytest)
library(RTMB)

if (at_home()) {

#########################################################
## Ising model
#########################################################

## States
n <- 3
s <- matrix(0, n, n)

## option 1: Lattice connectivity
e <- expand.grid(1:nrow(s), 1:ncol(s))
J0 <- as.matrix(dist(e)) == 1
## option 2: Dense connectivity
J1 <- J0 * 0 + 1

## Data and parameters
data <- list(J=J0, flag=1)
parameters <- list(s=s, beta=.1)

f <- function(pars) {
  getAll(data, pars)
  J <- J * beta
  s <- as.vector(s)
  if (flag)
    -s%*%(J%*%s)  ## few 'big' terms => slow SR
  else
    -sum(J*outer(s,s))  ## many 'small' terms => fast SR
}

## Integration grid
G <- c(-1, +1)

ArgComb <- expand.grid(
  perm=c(TRUE,FALSE), flag=c(TRUE,FALSE))

testFun <- function(perm, flag, dense, debug=FALSE) {
  data$J <<- if (dense) J1 else J0
  data$flag <<- flag
  integrate <- list(
    "s" = TMB::SR(G, discrete=TRUE, perm=perm, debug=debug)
  )
  obj <- MakeADFun(f, parameters, random="s", integrate=integrate)
  obj$fn(1) ## beta=1
}

ans0 <- Vectorize(testFun)(ArgComb$perm, ArgComb$flag, FALSE)
expect_equal(ans0, rep(-24.6945889885603, 4L), info="SR ising sparse")

ans1 <- Vectorize(testFun)(ArgComb$perm, ArgComb$flag, TRUE)
expect_equal(ans1, rep(-81.6931471805601, 4L), info="SR ising dense")

## Brute force:
## data$J <- J1
## myf <- function(s) {parameters$s[] <- s; parameters$beta <- 1; f(parameters)}
## e <- expand.grid(rep(list(G),length(parameters$s)))
## -log(sum(exp(-apply(e,1,myf))))

#########################################################
## Test random clique permutations
#########################################################
library(RTMB)
g <- function(x) sin(x[1])^2 * cos(x[2])^2 * exp(x[3])
G <- -2:2
nx <- 100
f <- function(x) {
  set.seed(1)
  s <- 0
  for (i in 1:(length(x)-2)) {
    s <- s + g(x[sample(i:(i+2), replace=TRUE)])
  }
  s
}
integrate=list("x"=TMB::SR(G, discrete=TRUE))
obj <- MakeADFun(f, numeric(nx), random="x", integrate=integrate)
## e <- expand.grid(rep(list(G),nx))
## -log(sum(exp(-apply(e,1,f))))
expect_equal(obj$fn(), -135.717352212, info = "SR random clique permutations")

#########################################################
## Test constant function
#########################################################
library(RTMB)
f <- function(x) 5
integrate=list("x"=TMB::SR(1:3, discrete=TRUE))
obj <- MakeADFun(f, numeric(4), random="x", integrate=integrate)
expect_equal(obj$fn(), -log(exp(-5) * 3^4) , info = "SR constant function")

#########################################################
## Binomial state space - density integrates to one
#########################################################
library(RTMB)
prob <- .3
f <- function(x) {
  s <- 0
  for (i in 1:length(x)) {
    xprev <- if (i==1) 0 else x[i-1]
    s <- s - dbinom(x[i] , 10 - xprev, prob, log=TRUE)
  }
  s
}
obj <- MakeADFun(f, numeric(100), random="x",
                 integrate=list("x"=TMB::SR(0:10, discrete=TRUE)))
expect_equal(obj$fn(), 0, info = "SR binomial state space")

#########################################################
## Binomial branching process - density integrates to one
#########################################################
library(RTMB)
prob <- .3
f <- function(x) {
  s <- 0
  for (i in 1:length(x)) {
    parent <- floor(i/2)
    xprev <- if (parent==0) 0 else x[parent]
    s <- s - dbinom(x[i] , 10 - xprev, prob, log=TRUE)
  }
  s
}
obj <- MakeADFun(f, numeric(100), random="x",
                 integrate=list("x"=TMB::SR(0:10, discrete=TRUE)))
expect_equal(obj$fn(), 0, info = "SR binomial branching process")

}
