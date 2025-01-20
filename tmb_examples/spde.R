# Illustration SPDE/INLA approach to spatial modelling via Matern correlation function
# Leukemia example from Lindgren et al 2011, JRSS-B

library(RTMB)
library(fmesher) # Replaces INLA. Used to create the mesh and precision matrix
library(Matrix)
library(spBayesSurv) # Contains the Leukemia dataset

data(LeukSurv)

## Spatial part part of model
loc <- cbind(LeukSurv$xcoord,LeukSurv$ycoord)  # Spatial coordinates
bnd1 <- fmesher::fm_nonconvex_hull(loc, convex=0.05)
bnd2 <- fmesher::fm_nonconvex_hull(loc, convex=0.25)

mesh <- fmesher::fm_mesh_2d(
  loc=loc,
  boundary=list(bnd1, bnd2),
  min.angle=24,
  max.edge=c(0.05, 0.2),
  cutoff=0.005,
  plot.delay=0.5
)

plot(mesh)

spde <- fmesher::fm_fem(mesh) # Calculate the sparse matrices c0,g1, g2 need for precission matrix

## Fixed effects part of model
X <- model.matrix( ~ 1 + sex + age + wbc + tpi, data = LeukSurv)

data <- list(time       = LeukSurv$time,
             notcens    = LeukSurv$cens,
             meshidxloc = mesh$idx$loc, ## RTMB: 1-based indexing !
             X          = as.matrix(X),
             c0 = spde$c0,
             g1 = spde$g1,
             g2 = spde$g2
)


parameters <- list(beta      = c(-5.0,0,0,0,0),
                   log_tau   = -2.0,
                   log_kappa = 2.5,
                   log_omega = -1,
                   x         = rep(0.0, nrow(spde$c0)) )

f <- function(parms) {
  getAll(parms, data)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)
  omega <- exp(log_omega)  ## Parameter of Weibull distribution
  
  Q <- tau^2*(kappa^4 * c0 + 2 * kappa^2 * g1 + g2)        ## GMRF prior
  nll <- -dgmrf(x, 0, Q, log=TRUE)        ## Negative log likelihood
  eta <- X%*%beta  +  x[meshidxloc]
  lambda <- exp(eta)
  for (i in 1:length(notcens)) {
    if (notcens[i])
      nll <- nll - dweibull(time[i],shape=omega,scale=1./lambda[i], log = TRUE)
    else
      nll <- nll - log(1.0-pweibull(time[i],shape=omega,scale=1./lambda[i])) # Survival prob.
  }
  nu <- 1.0            ## nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
  rho <- sqrt(8*nu)/kappa  ## Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
  ADREPORT(rho)
  
  return(nll)
}

obj <- MakeADFun(f, parameters, random="x")
opt <- nlminb(obj$par, obj$fn, obj$gr)


## Calculate standard deviations, and extract correlation distance (rho)
Rep <- sdreport(obj)
rho_est <- summary(Rep,"report")
