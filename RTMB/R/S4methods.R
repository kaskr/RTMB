setClass("advector") ## Virtual class

setMethod("Ops",
          signature("sparseMatrix", "advector"),
          function(e1, e2) callGeneric( magic(e1, TRUE) , e2) )

setMethod("Ops",
          signature("advector", "sparseMatrix"),
          function(e1, e2) callGeneric( e1, magic(e2, TRUE) ) )
