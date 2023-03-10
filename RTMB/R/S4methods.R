setMethod("show", "advector", function(object) print.advector(object) )

setAs("sparseMatrix", "adsparse",
      function(from) {
          x <- from
          x <- as(x, "generalMatrix")
          x <- as(x, "CsparseMatrix")
          new("adsparse", x=advector(x@x), i=x@i, p=x@p, Dim=x@Dim)
      })

##setClassUnion("advector_castable", c("advector", "numeric"))

##' @describeIn ADmatrix AD matrix exponential
setMethod("expm", "advector", function(x) math_expm(x))
##' @describeIn ADmatrix AD matrix exponential
setMethod("expm", "adsparse", function(x) math_expm(x))

## Methods sparseMatrix -> adsparse

setMethod("Ops",
          signature("sparseMatrix", "advector"),
          function(e1, e2) callGeneric( as(e1, "adsparse") , e2) )

setMethod("Ops",
          signature("advector", "sparseMatrix"),
          function(e1, e2) callGeneric( e1, as(e2, "adsparse") ) )

setMethod("Ops",
          signature("sparseMatrix", "adsparse"),
          function(e1, e2) callGeneric( as(e1, "adsparse") , e2) )

setMethod("Ops",
          signature("adsparse", "sparseMatrix"),
          function(e1, e2) callGeneric( e1, as(e2, "adsparse") ) )

## Methods adsparse
##' @describeIn ADmatrix AD sparse matrix dimension
setMethod("dim", "adsparse", function(x) x@Dim)

setMethod("Ops",
          signature("ad", "adsparse"),
          function(e1, e2) SparseArith2(advector(e1), e2, .Generic) )
setMethod("Ops",
          signature("adsparse", "ad"),
          function(e1, e2) SparseArith2(e1, advector(e2), .Generic) )
setMethod("Ops",
          signature("adsparse", "adsparse"),
          function(e1, e2) {
              if (!identical(e1@Dim, e2@Dim))
                  stop("non-conformable arguments")
              SparseArith2(e1, e2, .Generic)
          })

##' @describeIn ADmatrix AD matrix multiply
setMethod("%*%",
          signature("adsparse", "ad"),
          function(x, y) {
              y <- as.matrix(advector(y))
              if ( ncol(x) != nrow(y) )
                  stop("non-conformable arguments")
              SparseArith2(x, y, .Generic)
          })
##' @describeIn ADmatrix AD matrix multiply
setMethod("%*%",
          signature("ad", "adsparse"),
          function(x, y) {
              x <- as.matrix(advector(x))
              if ( ncol(x) != nrow(y) )
                  stop("non-conformable arguments")
              SparseArith2(x, y, .Generic)
          })
##' @describeIn ADmatrix AD matrix multiply
setMethod("%*%",
          signature("adsparse", "adsparse"),
          function(x, y) {
              if ( ncol(x) != nrow(y) )
                  stop("non-conformable arguments")
              SparseArith2(x, y, .Generic)
          })
##' @describeIn ADmatrix AD matrix multiply
setMethod("%*%",
          signature("ad", "ad"),
          function(x, y) {
              x <- as.matrix(x)
              y <- as.matrix(y)
              matmul(advector(x), advector(y))
          })
##' @describeIn ADmatrix AD matrix multiply
setMethod("tcrossprod", signature("advector"),
          function(x, y=NULL) {if (is.null(y)) y <- x; x %*% t(y)} )
##' @describeIn ADmatrix AD matrix multiply
setMethod( "crossprod", signature("advector"),
          function(x, y=NULL) {if (is.null(y)) y <- x; t(x) %*% y} )
##' @describeIn ADmatrix AD matrix inversion and solve
##' @param a matrix
##' @param b matrix, vector or missing
setMethod("solve",
          signature("advector", "ad."),
          function(a, b) {
              a <- as.matrix(a)
              ans <- matinv(advector(a))
              if (!missing(b)) {
                  b <- as.matrix(b)
                  ans <- ans %*% b
              }
              ans
          })
##' @describeIn ADmatrix AD matrix (or array) colsums
setMethod("colSums", signature(x="advector"),
          function(x) { apply(x, seq_len(length(dim(x)))[-1L], sum) } )
##' @describeIn ADmatrix AD matrix (or array) rowsums
setMethod("rowSums", signature("advector"),
          function(x) { apply(x, 1L, sum) } )
##' @describeIn ADmatrix AD matrix column bind
##' @param ... As \link[base]{cbind}
cbind.advector <- function (...) {
    args <- lapply(list(...), advector)
    ans <- do.call("cbind", lapply(args, unclass))
    class(ans) <- "advector"
    asS4(ans)
}
##' @describeIn ADmatrix AD matrix row bind
rbind.advector <- function (...) {
    args <- lapply(list(...), advector)
    ans <- do.call("rbind", lapply(args, unclass))
    class(ans) <- "advector"
    asS4(ans)
}

## Show general idea which is automated in 'distributions.R'
## First we generate the version we want for AD types (dot signifies 'default argument')
##' @describeIn Distributions AD implementation of \link[stats]{dnorm}
setMethod("dnorm", signature("ad", "ad.", "ad.", "logical."),
          function(x, mean, sd, log) {
              r <- (x - mean) / sd
              ans <- - .5 * r * r - log(sqrt(2*pi)) - log(sd)
              if (log) ans else exp(ans)
          })
## This matches 'too much', so we fix by adding a specialization:
##' @describeIn Distributions Default method
setMethod("dnorm", signature("num", "num.", "num.", "logical."),
          function(x, mean, sd, log) {
              stats::dnorm(x, mean, sd, log)
          })
## For S4 generics we add the OSA version like this:
##' @describeIn Distributions OSA implementation
setMethod("dnorm", "osa", function(x, mean, sd, log) {
    dGenericOSA(.Generic, x=x, mean=mean, sd=sd, log=log)
})
## For S4 generics we add the simref version like this:
##' @describeIn Distributions Simulation implementation. Modifies \code{x} and returns zero.
setMethod("dnorm", "simref", function(x, mean, sd, log) {
    ## works when x, mean or sd are simref
    if (inherits(mean, "simref")) {
        x <- x - mean
        mean <- 0
    }
    if (inherits(sd, "simref")) {
        x <- x / sd
        sd <- 1
    }
    dGenericSim(.Generic, x=x, mean=mean, sd=sd, log=log)
})

##' @describeIn Distributions Minimal AD implementation of \link[stats]{plogis}
setMethod("plogis", c("advector", "missing", "missing", "missing", "missing"),
          function(q) 1 / (1 + exp(-q) ) )
##' @describeIn Distributions Minimal AD implementation of \link[stats]{qlogis}
setMethod("qlogis", c("advector", "missing", "missing", "missing", "missing"),
          function(p) log( p / ( 1 - p ) ) )

## 'diag' needs patching.
## - base::diag works fine for AD matrix input (diagonal extraction and replacement)
## - However, matrix construction has issues

##' @describeIn ADconstruct Equivalent of \link[base]{diag}
##' @param x As \link[base]{diag}
##' @param nrow As \link[base]{diag}
##' @param ncol As \link[base]{diag}
setMethod("diag", signature(x="num.", nrow="num.", ncol="num."),
          function(x, nrow, ncol) {
              ans <- callNextMethod()
              if (ad_context()) ans <- advector(ans)
              ans
          })

##' @describeIn ADconstruct Equivalent of \link[base]{diag}
setMethod("diag", signature(x="advector", nrow="ANY", ncol="ANY"),
          function(x, nrow, ncol) {
              ## Diagonal extraction: base::diag works fine
              if (length(dim(x)) >= 2)
                  return(callNextMethod())
              ## Matrix creation
              ans <- advector(base::diag(seq_along(x), nrow=nrow, ncol=ncol))
              diag(ans) <- x
              ans
          })

## Constructors that need 'magic'
##' @describeIn ADconstruct Equivalent of \link[base]{numeric}
##' @param length As \link[base]{numeric}
setMethod("numeric", signature(length="num."),
          function(length) {
              ans <- callNextMethod()
              if (ad_context()) ans <- advector(ans)
              ans
          })
##' @describeIn ADconstruct Equivalent of \link[base]{matrix}
##' @param data As \link[base]{matrix}
##' @param byrow As \link[base]{matrix}
setMethod("matrix", signature(data="num."),
          function(data, nrow, ncol, byrow, dimnames) {
              ans <- callNextMethod()
              if (ad_context()) ans <- advector(ans)
              ans
          })
##' @describeIn ADconstruct Equivalent of \link[base]{matrix}
setMethod("matrix", signature(data="advector"),
          function(data, nrow, ncol, byrow, dimnames) {
              ans <- callNextMethod()
              class(ans) <- "advector"
              ans
          })
##' @describeIn ADconstruct Equivalent of \link[base]{array}
##' @param data As \link[base]{array}
##' @param dim As \link[base]{array}
##' @param dimnames As \link[base]{array}
setMethod("array", signature(data="num."),
          function(data, dim, dimnames) {
              ans <- callNextMethod()
              if (ad_context()) ans <- advector(ans)
              ans
          })
##' @describeIn ADapply As \link[base]{apply}
##' @param X As \link[base]{apply}
##' @param MARGIN As \link[base]{apply}
##' @param FUN As \link[base]{apply}
##' @param ... As \link[base]{apply}
setMethod("apply", signature(X="advector"),
          function (X, MARGIN, FUN, ...)  {
              ans <- callNextMethod()
              if (is.complex(ans))
                  class(ans) <- "advector"
              ans
          })
##' @describeIn ADapply As \link[base]{sapply}
##' @param simplify As \link[base]{sapply}
##' @param USE.NAMES As \link[base]{sapply}
setMethod("sapply", signature(X="advector"),
          function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
              ans <- callNextMethod()
              if (is.complex(ans))
                  class(ans) <- "advector"
              ans
          })

##' @describeIn ADvector Equivalent of \link[base]{ifelse}
##' @param test \code{logical} vector
##' @param yes \code{advector}
##' @param no \code{advector}
setMethod("ifelse", signature(test="num", yes="ad", no="ad"),
          function(test, yes, no) {
              yes <- advector(yes)
              no <- advector(no)
              ans <- callNextMethod()
              class(ans) <- "advector"
              ans
          })
##' @describeIn ADvector Default method
setMethod("ifelse", signature(test="num", yes="num", no="num"),
          function(test, yes, no) {
              base::ifelse(test, yes, no)
          })

##' @describeIn ADvector Equivalent of \link[base]{outer}
##' @param X As \link[base]{outer}
##' @param Y As \link[base]{outer}
setMethod("outer", signature(X="advector", Y="advector", FUN="missing"),
          function (X, Y) outer(X, Y, function(x, y) x * y))
