################################################################################
## This file contains:
## - Method to transform sparse Hessian such that all data terms are PD
################################################################################

## 'posfun' that works for matrices or scalars
posfun <- function(x, p=.5) {
    if (is.matrix(x)) {
        p * math_absm(x) + (1-p) * x
    } else {
        p * abs(x) + (1-p) * x
    }
}

## - Use RTMBp:::Term to mark terms to be corrected
## - Automatically handle the sign
## - Limitation: Must have same pattern as original
getPositiveHessian <- function(obj, p = .5) {
    ## Get pattern of existing hessian
    h <- obj$env$spHess(random=TRUE)
    ## Get tape of the existing hessian
    H <- GetTape(obj, "ADHess")
    ## Get corresponding pointer
    Hptr <- .pointer(environment(H)$mod)
    ## Locate independent variables on the tape
    invold <- getinvIndex(Hptr)
    ## Locate all data terms on the tape
    i <- findIndex(Hptr, "TermOp1")
    ## - Set the data terms as new independent variables
    ## - And make the old InvOp place holders 'replay persistent'
    setinvIndex(Hptr, i)
    InvPersistent(Hptr, TRUE)
    ## Get sparse Jacobian wrt all these data terms
    ## NOTE:
    ##   rows = Hessian non zeros
    ##   cols = data terms
    HH <- H$jacfun(sparse=TRUE)
    ## Restore H tape to original state
    setinvIndex(Hptr, invold)
    InvPersistent(Hptr, FALSE)
    ## We can get the term signs (cool?)
    signs <- HH$par()
    if (length(unique(signs)) != 1)
        warning("Multiple term signs: ", paste(unique(signs), collapse=" "))
    ## Find the InvOp placeholders on HH tape
    HHptr <- .pointer(environment(HH)$mod)
    inv_ <- findIndex(HHptr, "InvOp_")
    ## Determine and set new inv_index
    invnew <- integer(length(invold))
    invnew[order(invold)] <- inv_
    setinvIndex(HHptr, invnew)
    ## Remove persistent InvOp from HH tape
    InvPersistent(HHptr, FALSE)
    ## Helper that modifies individual columns of 'HH'
    modify <- function(jac, col) {
        ## indices into jac@i and jac@x defining current column 'col'
        ind <- (jac@p[col]+1):(jac@p[col+1])
        ## row indices of jacobian column = pointers into h@x of hessian
        tmp <- jac@i[ind]+1
        ## Absolute row indices of the block (zero based)
        i <- h@i[tmp]
        ## We can only handle two cases: dense or diag
        dim <- length(unique(i)) ## Sub matrix dimension
        ## Get the non-zeros of the block (multiplied by the sign)
        xnz <- jac@x[ind] * signs[col]
        if (dim == length(i)) {
            ## Diagonal case (easy elementwise posfun)
            posfun(xnz, p)
        }
        else if ( (dim * dim - dim) * .5 + dim == length(i)) {
            ## Dense case (spectral posfun)
            X <- advector(matrix(0, dim, dim))
            lt <- lower.tri(X, TRUE)
            X[lt] <- xnz
            X[!lt] <- t(X)[!lt] ## Shouldn't be needed
            posfun(X, p)[lt]
        }
        else stop("Can only handle 'diag' or 'dense' data terms")
    }
    ## R function that tapes the new Hessian non-zeros (same pattern as original)
    newHess <- function(x) {
        ## ---------------- Eval H *without* data
        Hptr <- .pointer(environment(H)$mod)
        TermsZero(Hptr, TRUE)
        ans1 <- H(x)
        TermsZero(Hptr, FALSE)
        ## ---------------- Eval H with *only* data
        jac <- HH(x)
        mcols <- list()
        for (j in 1:ncol(jac)) {
            mcols[[j]] <- modify(jac, j)
        }
        ans <- unlist(mcols)
        class(ans) <- "advector"
        jac@x[] <- ans
        one <- rep(1, ncol(jac))
        ans2 <- jac %*% one
        ## ---------------- Add the two
        ans1+ans2
    }
    MakeTape(newHess, H$par())
}

## Get a function to swap hessians
## Usage:
##    obj$env$altHess <- altHessFun(obj)
altHessFun <- function(obj, p=.5) {
    posHess <- getPositiveHessian(obj, p)
    ADHess2 <- environment(posHess)$mod$ptrTMB()
    attr(ADHess2, "rcpp") <- environment(posHess)$mod ## 'PROTECTS' ptrTMB !
    ADHess2$DLL <- "RTMBp"
    e <- environment(obj$env$spHess)
    ADHess <- e$ADHess
    attributes(ADHess2$ptr) <- attributes(ADHess$ptr)
    function(altHess=FALSE) {
        e$ADHess <- if (altHess) ADHess2 else ADHess
        NULL
    }
}
