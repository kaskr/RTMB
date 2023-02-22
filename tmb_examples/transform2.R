library(RTMB)

## Generate AR1 variables
set.seed(123)
n <- 1000
phi <- .6
u <- numeric(n)
u[1] = rnorm(1)
for(i in 2:n){
    u[i] = phi * u[i-1] + rnorm(1, sd = sqrt(1 - phi^2))
}

## Transform to beta
shape1 <- .5
shape2 <- 2
x <- qbeta(pnorm(u), shape1 = shape1, shape2 = shape2)

## Samples with noise. Note: we need *many* replicates per random
## effect ( or equivalently one replicate with *small* sd) for the
## Laplace approximation to work for this problem):
sd <- .005
y <- x + rnorm(n, sd=sd)

data=list(y=y)
parameters=list(phi=0, shape1=1, shape2=1, sd=1, u=u*0)

dar1 <- function(x, phi, log=FALSE) {
    ans <-
        dnorm(x[1], 0, 1, log=TRUE) +
        sum(dnorm( x[-1] - phi * x[-length(x)], 0, sqrt(1-phi*phi), log=TRUE ))
    if (log) ans else exp(ans)
}

func <- function(phi,
                 shape1,
                 shape2,
                 sd,
                 u) {
    y <- data$y
    y <- OBS(y) ## Optional
    res <- 0
    res <- res - dar1(u, phi, log=TRUE)
    unif <- pnorm(u)
    x <- qbeta(unif, shape1, shape2)
    res <- res - sum(dnorm(y, x, sd, log=TRUE))
    res
}
do.call(func,parameters)

obj <- MakeADFun(
    function(p)do.call(func,p),
    parameters=parameters,
    random="u"
)
obj$fn()
opt <- nlminb(obj$par, obj$fn, obj$gr, lower=1e-6)
sdr <- sdreport(obj)
summary(sdr,"fixed")
