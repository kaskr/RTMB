setClass("advector") ## Virtual class

setMethod("Ops",
          signature("sparseMatrix", "advector"),
          function(e1, e2) callGeneric( magic(e1, TRUE) , e2) )

setMethod("Ops",
          signature("advector", "sparseMatrix"),
          function(e1, e2) callGeneric( e1, magic(e2, TRUE) ) )

setClass("adsparse",
         slots=c(x="advector", i="integer", p="integer", Dim="integer"))

setAs("sparseMatrix", "adsparse",
      function(from) {
          x <- from
          x <- as(x, "generalMatrix")
          x <- as(x, "CsparseMatrix")
          new("adsparse", x=advector(x@x), i=x@i, p=x@p, Dim=x@Dim)
      })
