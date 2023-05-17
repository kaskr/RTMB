library(RTMB)

## Consider a case where the same objective function has to be evaluated for many different identically shaped datasets.
## We want to create a single tape and reuse it for all the data.
## DataEval is used to solve this problem. It 'abuses' a parameter for an index and allows evaluation of an arbitrary function with the index as argument. The evaluation will become part of the AD tape (with derivatives set to zero). See example below.

## For example 10 linear regressions
df <- data.frame(x=1:100,y=1:100 + rnorm(100))
dfs <- split(df, gl(10,10))

## Setup a parameter list for the first data
plist <- as.relistable(list(a=0, b=0, sd=1, i=1))

## Objective function for one chunk of data
F <- MakeTape(function(x) {
    ## Get parameter list
    parms <- relist(x)
    ## Fetch data from R (x and y)
    x <- DataEval( function(i) dfs[[i]]$x, parms$i)
    y <- DataEval( function(i) dfs[[i]]$y, parms$i)
    -sum(dnorm(y , parms$a * x + parms$b, parms$sd, log=TRUE))
}, unlist(plist))

F$simplify() ## Simplify the operation sequence if possible
F <- F$atomic() ## Turn into atomic function

## Aggregate tape by summing over 'i'
obj <- MakeADFun(function(x) {
    s <- 0
    for (i in 1:10) s <- s + F(c(x,i))
    s
}, unlist(plist)[1:3])

opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
sdreport(obj)
