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

##setClassUnion("advector_convertable", c("advector", "numeric"))
setClassUnion("advector_convertable", c("advector", "array", "numeric", "logical"))

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
          signature("advector_convertable", "adsparse"),
          function(e1, e2) SparseArith2(advector(e1), e2, .Generic) )
setMethod("Ops",
          signature("adsparse", "advector_convertable"),
          function(e1, e2) SparseArith2(e1, advector(e2), .Generic) )
setMethod("Ops",
          signature("adsparse", "adsparse"),
          function(e1, e2) SparseArith2(e1, e2, .Generic) )

setMethod("%*%",
          signature("advector_convertable", "advector_convertable"),
          function(x, y) {
              x <- as.matrix(x)
              y <- as.matrix(y)
              matmul(advector(x), advector(y), method="atomic")
          })

setMethod("tcrossprod", signature("advector"),
          function(x, y=NULL) {if (is.null(y)) y <- x; x %*% t(y)} )
setMethod( "crossprod", signature("advector"),
          function(x, y=NULL) {if (is.null(y)) y <- x; t(x) %*% y} )

Sig3 <- signature("advector_convertable", "advector_convertable", "advector_convertable", "logical")
setMethod("dnorm", Sig3,
          function(x, mean = 0, sd = 1, log = FALSE) {
              r <- (x - mean) / sd
              ans <- - .5 * r * r - log(sqrt(2*pi)) - log(sd)
              if (log) ans else exp(ans)
          })
