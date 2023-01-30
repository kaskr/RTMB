setClass("advector") ## Virtual class

setClass("adsparse",
         slots=c(x="advector", i="integer", p="integer", Dim="integer"))

setAs("sparseMatrix", "adsparse",
      function(from) {
          x <- from
          x <- as(x, "generalMatrix")
          x <- as(x, "CsparseMatrix")
          new("adsparse", x=advector(x@x), i=x@i, p=x@p, Dim=x@Dim)
      })

##setClassUnion("advector_castable", c("advector", "numeric"))

## Helpers to setMethod
## Match numeric like objects that are 'AD castable' via advector()
setClassUnion("num", c("array", "numeric", "logical"))
setClassUnion("num.", c("num", "missing"))
setClassUnion("ad",  c("advector", "num"))
setClassUnion("ad.", c("advector", "num."))
setClassUnion("logical.", c("logical", "missing"))


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
              matmul(advector(x), advector(y), method="atomic")
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
