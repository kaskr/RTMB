################################################################################
## This file contains:
## - Experimental features. Not yet exported!
################################################################################

## Interface to newton::Tag
Tag <- function(x) {
    if (inherits(x, "advector"))
        LowRankTag(x)
    else
        x
}

## Post vectorization heuristic (Move to TMB if useful)
vectorize <- function(obj, verbose=FALSE) {
    if (verbose) info1 <- TMB:::info(obj$env$ADFun)
    TMB:::TransformADFunObject(obj$env$ADFun, method = "accumulation_tree_split") 
    TMB:::TransformADFunObject(obj$env$ADFun, method = "reorder_sub_expressions")
    TMB:::TransformADFunObject(obj$env$ADFun, method = "optimize")
    TMB:::TransformADFunObject(obj$env$ADFun, method = "fuse_and_replay")
    if (verbose) info2 <- TMB:::info(obj$env$ADFun)
    if (verbose) {
        cat("Before:\n")
        print(unlist(info1[-1]))
        cat("After:\n")
        print(unlist(info2[-1]))        
    }
    invisible(obj)
}

## Decompose F(x)=L(T(x)) where L is *linear* and *maximal*
## Or apply an 'inner' scale transformation S(F(x)) using linearity of L.
term_split <- function(F, scale=NULL) {
  ptr <- .pointer(environment(F)$mod)
  nodes <- get_term_nodes(ptr)
  L <- .copy(environment(F)$mod)
  T <- .copy(environment(F)$mod)
  vars <- op2var(ptr, nodes)
  setinvIndex(.pointer(L), vars)
  inactivate(.pointer(L), nodes)
  L <- .expose(L)
  L$simplify("eliminate")
  setdepIndex(.pointer(T), vars)
  T <- .expose(T)
  T$simplify("eliminate")
  ## Linearize:
  ## L(x) = L(x0) + L'(x0) * (x - x0)
  zero <- rep(0, length(L$par()))
  L0 <- L(zero)
  J <- L$jacfun(sparse=TRUE)(zero)
  if (!is.null(scale)) {
    J <- AD(J)
    J@x[] <- scale[J@i+1L] * J@x
    ans <- MakeTape(function(x) scale * L0 + AD(J) %*% T(x), T$par())
    return(ans)
  }
  L <- MakeTape(function(x) L0 + AD(J) %*% x, L$par())
  list(T=T, L=L)
}
