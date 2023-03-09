# Illustration SPDE/INLA approach to spatial modelling via Matern correlation function
# Leukemia example from Lindgren et al 2011, JRSS-B
# http://www.r-inla.org/examples/case-studies/lindgren-rue-and-lindstrom-rss-paper-2011

library(RTMB)
library(Matrix)

## get cached objects - See 'spde_mesh.R'
##  'inla_mesh'
##  'inla_spde'
##  'Leuk'
load("spde_mesh.RData")

## Fixed effects part of model
X <- model.matrix( ~ 1 + sex + age + wbc + tpi, data = Leuk)

data <- list(time       = Leuk$time,
             notcens    = Leuk$cens,
             meshidxloc = inla_mesh$idx$loc, ## RTMB: 1-based indexing !
             X          = as.matrix(X))

## SPDE part: builds 3 components of Q (precision matrix)
data$spde <- inla_spde$param.inla[c("M0","M1","M2")]  # Encapsulation of 3 matrices
n_s <- nrow(data$spde$M0)                             # Number of points in mesh (including supporting points)

parameters <- list(beta      = c(-5.0,0,0,0,0),
                   log_tau   = -2.0,
                   log_kappa = 2.5,
                   log_omega = -1,
                   x         = rep(0.0, n_s) )

## Objective function
Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2    ## M0=G0, M1=G1, M2=G2
}
f <- function(parms) {
    getAll(parms, data)
    tau <- exp(log_tau)
    kappa <- exp(log_kappa)
    omega <- exp(log_omega)  ## Parameter of Weibull distribution
    ## GMRF prior
    Q <- Q_spde(spde, kappa)
    nll <- -dgmrf(x, 0, Q, log=TRUE)        ## Negative log likelihood
    Xbeta <- X %*% beta
    eta <- Xbeta + x[meshidxloc] / tau
    lambda <- exp(eta)
    t_omega <- time^omega
    S <- exp(-lambda*t_omega)                ## Survival function
    f <- lambda * omega * t_omega / time * S ## Weibull density
    ## Likelihood contribution depends on truncation status
    for (i in 1:length(notcens)) {
        if (notcens[i])
            nll <- nll - log(f[i])
        else
            nll <- nll - log(S[i])
    }
    nu <- 1.0            ## nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
    rho <- sqrt(8*nu)/kappa  ## Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
    ADREPORT(rho)

    return(nll)
}

f(parameters)

## Phase 1: Fit non-spatial part first to get good starting values for fixed effects
not_phase1 <- list(log_tau   = as.factor(NA),
                   log_kappa = as.factor(NA),
                   x         = factor(rep(NA, n_s)) )
obj <- MakeADFun(f, parameters, map=not_phase1)
opt1 <- nlminb(obj$par, obj$fn, obj$gr)

## Modify starting values after phase 1
parameters <- list(beta      = opt1$par[1:5],
                   log_tau   = -2.0,
                   log_kappa =  2.5,
                   log_omega = opt1$par["log_omega"],
                   x         = rep(0.0, n_s))

## Phase 2: Include spatial part. Use starting values from phase 1
obj <- MakeADFun(f, parameters, random="x")
L   <- c(-7, -1, -1, -1, -1, -3.0, 2.0, log(0.1) )
U   <- c(-4,  1,  1,  1,  1, -1.0, 3.0, log(10.0))
opt <- nlminb(obj$par, obj$fn, obj$gr, lower=L, upper=U)

# Calculate standard deviations, and extract rho
Rep <- sdreport(obj)
rho_est <- summary(Rep,"report")
