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
