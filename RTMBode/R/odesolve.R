Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X + Y * Z
    dY <- b * (Y - Z)
    dZ <- -X * Y + c * Y - Z
    list(c(dX, dY, dZ))
  })
}

## Example tape
extape <- function() {
    parameters <- c(a = -8/3, b = -10, c =  28)
    state <- c(X = 1, Y = 1, Z = 1)
    MakeTape(function(typ) {
        t <- typ[1]
        yp <- typ[-1]
        y <- yp[1:3]
        p <- yp[-(1:3)]
        Lorenz(t,y,p)[[1]]
    } ,  c(t=0, state, parameters) )
}

setTape <- function(F) {
    ##domain <- environment(F)$mod$domain()
    range <- environment(F)$mod$range()
    x <- F$par()
    ##y <- F(F$par())
    y <- double(range)
    ##ptr <- RTMB:::.pointer(environment(F)$mod)
    ptr <- environment(F)$mod$ptrTMB()$ptr
    .Call(set_pointers, ptr, x, y)
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

makef <- function(F) {
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
        ## Hack to update internal parms
        ## (FIXME: not great when there's an active ad_context !)
        tmp <- c(times[1], augstate, parms)
        F(tmp)
        ## We should always setTape before calling ode:
        setTape(F)
        sol <<- deSolve::ode(augstate,
                             times,
                             func = "desolve_derivs", parms=NULL,
                             dllname = "RTMBode")
        ## cat("Calling sol "); print(ncol(sol))
        attr(sol, "x") <<- x
        NULL
    }
    Df <- NULL ## Not yet initialized
    f <- function(x) {
        if (inherits(x, "advector")) {
            DataEval(makef(F), x)
        } else {
            updateSolution(x)
            sol[, solcols]
        }
    }
    attr(f, "reverse") <- function(x, y, w) {
        if (is.null(Df))
            Df <<- makef(augment(F))
        t(w) %*% matrix(Df(x), nrow=length(y))
    }
    attr(f, "name") <- "ODE"
    f
}
