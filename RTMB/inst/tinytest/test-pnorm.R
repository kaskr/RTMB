## Test that pnorm works as expected
library(RTMB)

################################################################################
## Test all pnorm combinations
################################################################################

for (lower.tail in c(TRUE, FALSE)) {
  for (log.p in c(FALSE, TRUE)) {
    x <- 1:10
    mu <- 10:1
    sd <- c(2,3) ## recycling
    f <- function(p) pnorm(x, p*mu, p*sd, lower.tail, log.p)
    F <- MakeTape(f, numeric(1))
    info <- paste("pnorm",
                  "lower.tail=", lower.tail,
                  "log.p=", log.p)
    expect_equal(F(1),
                 f(1),
                 info=paste(info, "tape eval"))
    expect_equal(F$jacobian(1),
                 numDeriv::jacobian(f, 1),
                 info=paste(info, "jacobian eval"))
    expect_equal(F$jacfun()(1),
                 F$jacobian(1),
                 info=paste(info, "jacobian replay"))
  }
}
