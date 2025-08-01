## Test that Bessel functions work as expected
library(RTMB)

################################################################################
## Test all bessel combinations
## (Skipping derivative testing for now - see #65)
################################################################################

for (FUN in c("besselK", "besselI")) {
  for (expo in c(FALSE, TRUE)) {
    x <- 1:10
    nu <- 10:1
    f <- function(p) get(FUN)(p*x, p*nu, expo)
    F <- MakeTape(f, numeric(1))
    expect_equal(F(1),
                 f(1),
                 info=paste(FUN, "expo=", expo))
    ## expect_equal(F$jacobian(1),
    ##              numDeriv::jacobian(f, 1),
    ##              info=paste(FUN, "expo=", expo, "jacobian"))
  }
}

for (FUN in c("besselJ", "besselY")) {
  x <- 1:10
  nu <- 10:1
  f <- function(p) get(FUN)(p*x, p*nu)
  F <- MakeTape(f, numeric(1))
  expect_equal(F(1),
               f(1),
               info=paste(FUN, "expo=", expo))
  ## expect_equal(F$jacobian(1),
  ##              numDeriv::jacobian(f, 1),
  ##              info=paste(FUN, "expo=", expo, "jacobian"))
}
