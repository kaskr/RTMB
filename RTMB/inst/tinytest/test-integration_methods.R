## Testing arguments 'integrate' and 'intern'
library(RTMB)
formals(MakeADFun)$silent <- TRUE

dat <- list(x = c(-0.5605, -0.2302, 1.5587))
par <- list(u = c(0, 0, 0), sd=1, mu=0)
map <- list(sd=factor(NA))
what <- "Est"
f <- function(par) {
  getAll(dat, par)
  nll <- -sum(dnorm(u, mu, log=TRUE))
  nll <- nll - sum(dnorm(x, mean=u, sd=sd, log=TRUE))
  nll
}

test <- function(...) {
  obj <- MakeADFun(f, par, random="u", map=map, ...)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  summary(sdr)
}

TMB.baseline <- test()
other <- list(
  intern = test(intern=TRUE),
  GK = test(integrate=list("u"=TMB:::GK())),
  LA = test(integrate=list("u"=TMB:::LA())),
  SR = test(integrate=list("u"=TMB:::SR(seq(-5,5,length=100))))
)
res <- lapply(other, function(x) x - TMB.baseline)

tol <- 1e-6
expect_true( all(abs(res$intern) < tol) , info="TMB intern")
expect_true( all(abs(res$GK) < tol) , info="GK")
expect_true( all(abs(res$LA) < tol) , info="LA")
expect_true( all(abs(res$SR) < tol) , info="SR")
