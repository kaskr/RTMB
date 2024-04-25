## Leukemia example from Lindgren et al 2011, JRSS-B
## http://www.r-inla.org/examples/case-studies/lindgren-rue-and-lindstrom-rss-paper-2011
library(RTMBp)
library(Matrix)

## SPDE data from INLA
##  'spde' : List of matrices M0, M1, M2 (Converted to compressed symmetric)
##  'Leuk' : data.frame
##  'loc'  : Integer code that maps Leuk rows to spatial grid id
load("spde_data.RData")

## Fixed effects part of model
X <- model.matrix( ~ 1 + sex + age + wbc + tpi, data = Leuk)

data <- list(time       = Leuk$time,
             notcens    = Leuk$cens,
             meshidxloc = loc,
             X          = as.matrix(X))

## SPDE part: builds Q (precision matrix) from M0, M1, M2.
Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2 ## M0=G0, M1=G1, M2=G2
}

n_s <- nrow(spde$M0) # Number of points in mesh (including supporting points)
parameters <- list(beta      = c(-5.0,0,0,0,0),
                   log_tau   = -2.0,
                   log_kappa = 2.5,
                   log_omega = -1,
                   x         = rep(0.0, n_s) )

## Objective function
f <- function(parms) {
    getAll(parms, data)
    tau <- exp(log_tau)
    kappa <- exp(log_kappa)
    omega <- exp(log_omega)  ## Parameter of Weibull distribution
    ## GMRF prior
    Q <- Q_spde(spde, kappa)
    nll <- -dgmrf(x, 0, Q, log=TRUE)
    Xbeta <- X %*% beta
    eta <- Xbeta + x[meshidxloc] / tau
    lambda <- exp(eta)
    t_omega <- time^omega
    S <- exp(-lambda*t_omega) ## Survival function
    f <- lambda * omega * t_omega / time * S ## Weibull density
    ## Likelihood contribution depends on truncation status
    for (i in 1:length(notcens)) {
        if (notcens[i])
            nll <- nll - log(f[i])
        else
            nll <- nll - log(S[i])
    }
    nu <- 1.0 ## nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
    rho <- sqrt(8*nu)/kappa  ## Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
    ADREPORT(rho)
    nll
}

## Check tape
obj <- MakeADFun(f, parameters, silent=TRUE)
expect_equal(f(parameters), obj$fn())

## Fit model (NaNs during minimization expected - supress warnings)
obj <- MakeADFun(f, parameters, random="x", silent=TRUE)
opt <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr))

# Calculate standard deviations, and extract rho
rep <- sdreport(obj)

## Check
tab <- summary(rep, c("fixed", "report"))

tab.expected <-
structure(c(-5.68807047389909, 0.0714771144613715, 0.0326402402338425, 
0.00306919177460965, 0.0248088827166649, -2.44407701215713, 2.50660923625429, 
-0.518240883366069, 0.230642020470054, 0.225699157756241, 0.0692812590113958, 
0.00223230252873417, 0.000456785153247373, 0.00988478153762717, 
0.531820909501231, 0.507805126113752, 0.0273250883581759, 0.117121200291926
), dim = c(9L, 2L), dimnames = list(c("beta", "beta", "beta", 
"beta", "beta", "log_tau", "log_kappa", "log_omega", "rho"), 
c("Estimate", "Std. Error")))

expect_true(rep$pdHess, info="Hessian positive definite")
expect_true(max(abs(rep$gradient.fixed)) < 1e-1, info="Gradient close to zero")
expect_equal(tab, tab.expected, tol=1e-4, info="Parameter estimates and std errors")
