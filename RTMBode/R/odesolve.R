setTape <- function(F, parms) {
    ptr <- RTMB:::.pointer(environment(F)$mod)
    X <- RTMB:::ptr_getx(ptr)
    Y <- RTMB:::ptr_gety(ptr)
    nstate <- attr(Y, "size")
    nparms <- attr(X, "size") - (nstate + 1L)
    .Call(set_pointers, ptr, X, Y, nstate, nparms)
    .Call(set_parms, parms)
}

## State augmentation (inital state and parameter derivatives)
## DEFINITION:
## - The tape represents the generator, i.e. the mapping (t,y,parms) -> f(t,y)
## - The augmented tape must also represent the generator!
## - It follows that the augmented tape must have:
##    Domain = 1 + length(augstate) + length(parms)
##    Range  = length(augstate)
## - We want the augmented tape to only represent the minimum number of states necessary:
## LAYOUT:
##   ORDER 0: nstate
##   ORDER 1: nstate + nstate * (nstate + nparms)
##   ORDER 2: nstate + nstate * (nstate + nparms) + nstate * (nstate + nparms)^2
addInfo <- function(F, times) {
    domain <- environment(F)$mod$domain()
    range <- environment(F)$mod$range()
    ## ORDER 0
    info <- list(
        order = 0,
        nparms = domain - range - 1,
        nstate = range,
        naugstate = range,
        state = head(tail(F$par(), -1), range),
        augstate = head(tail(F$par(), -1), range)
    )
    info$par <- F$par()
    info$state <- head(tail(info$par, -1), info$nstate)
    info$parms <- tail(info$par, info$nparms)
    info$times <- times
    attr(F, "info") <- info
    F
}
augment <- function(F, oldpar=F$par()) {
    info <- attr(F, "info")
    if (is.null(info)) stop("'F' must have 'info' attribute")
    newinfo <- info ## Updated info must be attached to output !
    newinfo$order <- with(info, order + 1)
    newinfo$naugstate <- with(info, naugstate + nstate*(nparms+nstate)^(order+1))
    ## When differentiating the system equation wrt a parameter, we
    ## need to locate the derivative in the augmented state vector
    ## (messy!).
    ## In general - see above 'LAYOUT' (i=state index, k=parameter index):
    dotIndex <- function(i, k) {
        gs <- nstate * (nstate + nparms)^(0:(info$order+1)) ## 'Group sizes'
        stride <- rep(gs, gs)
        i + stride[i] * k
    }
    nstate <- info$nstate
    nparms <- info$nparms
    augstate <- c(y=info$state,
                  dy=as.numeric(diag(nstate)),
                  dy=rep(0, newinfo$naugstate - nstate - nstate^2 ))
    newinfo$augstate <- augstate
    newpar <- c(t=0, augstate, info$parms)
    dotTab <- outer(1:info$naugstate, 1:(nstate+nparms), dotIndex)
    g <- function(x) {
        ## Step 1. Evaluate the old tape (to get time derivatives of old states)
        ## Pointers to extract old parameter from new
        iold <- seq_along(x)
        iold <- c(head(iold, 1+info$naugstate) , tail(iold, nparms))
        xold <- x[iold]
        res1 <- F(xold) ## oldstate time derivative
        ## Step 2. Append time derivatives of the new states
        nappend <- newinfo$naugstate - info$naugstate ## Number of states to append
        ntail <- with(info, nstate * (nstate + nparms)^order) ## Old states for which we need derivatives
        Fproj <- MakeTape(function(x) tail(F(x), ntail), F$par())
        J <- Fproj$jacfun(sparse=TRUE)(xold) ## ntail x length(xold)
        ## We need two parts of this jacobian: (Using layout of 'xold')
        tmp <- seq_along(xold)[-1] ## without 't'
        Jparms <- J[ , tail(tmp, nparms)]  ## Derivatives of system equation wrt parms
        Jstate <- J[ , head(tmp, -nparms)] ## Derivatives of system equation wrt states
        xdot <- array(x[dotTab+1], dim(dotTab)) ## xold[-1] x (nstate + nparms)
        res2 <- Jstate %*% xdot ## ntail x (nstate+nparms)
        res2[,-(1:nstate)] <- res2[,-(1:nstate)] + Jparms
        c(res1, as.numeric(res2))
    }
    G <- MakeTape(g, newpar)
    G$simplify()
    attr(G, "info") <- newinfo
    G
}

ODEadjoint <- function(F, ...) {
    info <- attr(F, "info")
    if (is.null(info)) stop("Please 'addInfo()' to 'F'")
    times <- info$times
    ##nstate <- length(state)
    nparms <- info$nparms
    sol <- NULL
    solcols <-
        local( tail( seq_len(naugstate)+1,  nstate * (nstate + nparms)^order), info)
    updateSolution <- function(x) {
        if (identical(attr(sol, "x"), x))
            return (NULL)
        state <- head(x, -nparms)
        augstate <- info$augstate
        augstate[seq_along(state)] <- state
        parms <- tail(x, nparms)
        ## We should always setTape before calling ode:
        setTape(F, parms)
        sol <<- deSolve::ode(augstate,
                             times,
                             func = "desolve_derivs",
                             parms=NULL,
                             dllname = "RTMBode",
                             ...)
        ## cat("Calling sol "); print(ncol(sol))
        attr(sol, "x") <<- x
        NULL
    }
    Df <- NULL ## Not yet initialized
    f <- function(x) {
        updateSolution(x)
        sol[, solcols]
    }
    reverse <- function(x, y, w) {
        if (is.null(Df))
            Df <<- ODEadjoint(augment(F), ...)
        t(w) %*% matrix(Df(x), nrow=length(y))
    }
    RTMB::ADjoint(f, reverse, "ODE")
}

func2tape <- function(func, y, parms) {
    nstate <- length(y)
    if (is.list(parms)) {
        skeleton <- as.relistable(parms)
        parms <- unlist(skeleton)
    } else {
        skeleton <- NULL
    }
    x <- numeric(1 + nstate + length(parms))
    names(x) <- c("t", names(y), names(parms))
    MakeTape(function(typ) {
        t <- typ[1]
        yp <- typ[-1]
        y <- yp[seq_len(nstate)]
        p <- yp[-seq_len(nstate)]
        if (!is.null(skeleton)) {
            p <- relist(p, skeleton)
        }
        func(t, y, p)[[1]]
    }, x)
}

##' ODE solver via 'deSolve' with 'RTMB' autodiff capabilities.
##'
##' This `ode` solver is essentially a wrapper around the corresponding function from the `deSolve` package. It adds the following extra features:
##' - Autodiff adjoint code so that ODE solving can be used as part of general gradient based optimization.
##' - Faster ODE solving using `RTMB` 'tapes' to eliminate R interpreter overhead.
##' @param y the initial (state) values for the ODE system, a vector (see `deSolve` package).
##' @param times time sequence for which output is wanted (see `deSolve` package).
##' @param func an R-function that computes the values of the derivatives in the ODE system (see `deSolve` package).
##' @param parms parameters passed to `func` (see `deSolve` package).
##' @param method the integrator to use (see `deSolve` package).
##' @param ... additional arguments passed to the integrator (see `deSolve` package).
##' @return Solution matrix with time as the first column (see `deSolve` package).
##' @examples
##' require(RTMB)
##' ## Lotka-Volterra example from 'deSolve' manual
##' LVmod <- function(Time, State, Pars) {
##'     with(as.list(c(State, Pars)), {
##'         Ingestion <- rIng * Prey * Predator
##'         GrowthPrey <- rGrow * Prey * (1 - Prey/K)
##'         MortPredator <- rMort * Predator
##'         dPrey <- GrowthPrey - Ingestion
##'         dPredator <- Ingestion * assEff - MortPredator
##'         return(list(c(dPrey, dPredator)))
##'     })
##' }
##' pars <- c(rIng = 0.2, # /day, rate of ingestion
##'           rGrow = 1.0, # /day, growth rate of prey
##'           rMort = 0.2 , # /day, mortality rate of predator
##'           assEff = 0.5, # -, assimilation efficiency
##'           K = 10) # mmol/m3, carrying capacity
##' yini <- c(Prey = 1, Predator = 2)
##' times <- seq(0, 200, by = 1)
##' ## Simulate ODE with measurement noise
##' set.seed(1)
##' obs <- deSolve::ode(func = LVmod, y = yini, parms = pars, times = times)[,-1]
##' obs <- obs + rnorm(length(obs), sd=1)
##' ## Likelihood function
##' likfun <- function(p) {
##'     getAll(p)
##'     obs <- OBS(obs)
##'     sol <- ode(func = LVmod, y = yini, parms = pars, times = times, atol=1e-8, rtol=1e-8)
##'     obs %~% dnorm(mean=sol[,-1], sd=sdobs)
##' }
##' ## Initial guess
##' p <- list(pars=pars*1.5, yini=yini*1.5, sdobs=1.5)
##' ## Parameter estimation
##' obj <- MakeADFun(likfun, p, silent=TRUE)
##' system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
##' (sdr <- sdreport(obj))
##' ## as.list(sdr, "Est")
##' ## as.list(sdr, "Std")
ode <- function (y, times, func, parms, method=NULL, ...) {
    F <- if (inherits(func, "Tape"))
             func
         else
             func2tape(func, y, parms)
    "c" <- ADoverload("c")
    while (is.list(parms)) parms <- do.call("c", parms)
    ad.case <- inherits(y, "advector") || inherits(parms, "advector")
    if (!ad.case) {
        setTape(F, parms)
        deSolve::ode(y = y,
                     times = times,
                     func = "desolve_derivs",
                     parms = NULL,
                     method = method,
                     dllname = "RTMBode",
                     ...)
    } else {
        F <- addInfo(F, times=times)
        F <- ODEadjoint(F, method=method, ...) ## Attach adjoint code
        x <- c(y, parms)
        ans <- F(x)
        ans <- cbind(times, ans)
        colnames(ans) <- c("time", names(y))
        ans
    }
}
