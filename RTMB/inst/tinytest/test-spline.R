library(RTMB)

################################################################################
## Test RTMB::splinfun 'method' consistent with stats::splinefun
################################################################################

## stats::splinefun
x <- 1:10
f1 <- splinefun(sin(x), method="fmm")
f2 <- splinefun(sin(x), method="natural")
f3 <- splinefun(sin(c(x, x[1])), method="periodic")
xnew <- seq(-2,12,length=100)

## RTMB::splinefun
T1 <- MakeTape(function(x)splinefun(sin(x), method="fmm")(xnew), numeric(10))
T2 <- MakeTape(function(x)splinefun(sin(x), method="natural")(xnew), numeric(10))
T3 <- MakeTape(function(x)splinefun(sin(c(x,x[1])), method="periodic")(xnew), numeric(10))

expect_equal(f1(xnew), T1(x), info="fmm")
expect_equal(f2(xnew), T2(x), info="natural")
expect_equal(f3(xnew), T3(x), info="periodic")

################################################################################
## Test constant spline with parameter dependent input
################################################################################

s <- splinefun(sin(1:10))
T <- MakeTape(s, 0)
expect_equal(Vectorize(T)(xnew), s(xnew), info="constant spline")
