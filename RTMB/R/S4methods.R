setMethod("show", "advector", function(object) print.advector(object) )

setAs("sparseMatrix", "adsparse",
      function(from) {
          x <- from
          x <- as(x, "generalMatrix")
          x <- as(x, "CsparseMatrix")
          new("adsparse", x=advector(x@x), i=x@i, p=x@p, Dim=x@Dim)
      })

##setClassUnion("advector_castable", c("advector", "numeric"))


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

setMethod("Ops",
          signature("ad", "adsparse"),
          function(e1, e2) SparseArith2(advector(e1), e2, .Generic) )
setMethod("Ops",
          signature("adsparse", "ad"),
          function(e1, e2) SparseArith2(e1, advector(e2), .Generic) )
setMethod("Ops",
          signature("adsparse", "adsparse"),
          function(e1, e2) SparseArith2(e1, e2, .Generic) )

setMethod("%*%",
          signature("ad", "ad"),
          function(x, y) {
              x <- as.matrix(x)
              y <- as.matrix(y)
              matmul(advector(x), advector(y))
          })

setMethod("tcrossprod", signature("advector"),
          function(x, y=NULL) {if (is.null(y)) y <- x; x %*% t(y)} )
setMethod( "crossprod", signature("advector"),
          function(x, y=NULL) {if (is.null(y)) y <- x; t(x) %*% y} )

## Show general idea which is automated in 'distributions.R'
## First we generate the version we want for AD types (dot signifies 'default argument')
setMethod("dnorm", signature("ad", "ad.", "ad.", "logical."),
          function(x, mean, sd, log) {
              r <- (x - mean) / sd
              ans <- - .5 * r * r - log(sqrt(2*pi)) - log(sd)
              if (log) ans else exp(ans)
          })
## This matches 'too much', so we fix by adding a specialization:
setMethod("dnorm", signature("num", "num.", "num.", "logical."),
          function(x, mean, sd, log) {
              stats::dnorm(x, mean, sd, log)
          })
## For S4 generics we add the OSA version like this:
setMethod("dnorm", "osa", function(x, mean, sd, log) {
    dGenericOSA(.Generic, x=x, mean=mean, sd=sd, log=log)
})
## For S4 generics we add the simref version like this:
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

## 'diag' needs patching.
## - base::diag works fine for AD matrix input (diagonal extraction and replacement)
## - However, matrix construction has issues

setMethod("diag", signature(x="num.", nrow="num.", ncol="num."),
          function(x, nrow, ncol) {
              ans <- callNextMethod()
              if (ad_context()) ans <- advector(ans)
              ans
          })

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

setMethod("apply", signature(X="advector"),
          function (X, MARGIN, FUN, ...)  {
              ans <- callNextMethod()
              if (is.complex(ans))
                  class(ans) <- "advector"
              ans
})
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
setMethod("ifelse", signature(test="logical", yes="advector", no="advector"),
          function(test, yes, no) {
              ans <- callNextMethod()
              class(ans) <- "advector"
              ans
          })
