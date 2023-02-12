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

## Transform to gamma
shape <- 2
scale <- 3
x <- qgamma(pnorm(u), shape=shape, scale=scale)

## Samples with noise
sd <- 2
y <- x + rnorm(n, sd=sd)

data <- list(y=y)
parameters <- list(phi=0, shape=1, scale=1, sd=1, u=u*0)

dar1 <- function(x, phi, log=FALSE) {
    ans <-
        dnorm(x[1], 0, 1, log=TRUE) +
        sum(dnorm( x[-1] - phi * x[-length(x)], 0, sqrt(1-phi*phi), log=TRUE ))
    if (log) ans else exp(ans)
}

func <- function(phi,
                 shape,
                 scale,
                 sd,
                 u) {
    y <- data$y
    res <- 0
    res <- res - dar1(u, phi, log=TRUE)
    unif <- pnorm(u)
    x <- qgamma(unif, shape, scale=scale)
    res <- res - sum(dnorm(y, x, sd, log=TRUE))
    res
}
do.call(func,parameters) ## 27638.34

obj <- MakeADFun(
    function(p)do.call("func",p),parameters=parameters,
    random="u"
    )
obj$fn()
opt <- nlminb(obj$par,obj$fn,obj$gr)
sdr <- sdreport(obj)
summary(sdr,"fixed")
