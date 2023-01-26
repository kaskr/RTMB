setClass("advector") ## Virtual class

setMethod("*",
          signature("sparseMatrix", "advector"),
          function(e1, e2) magic(e1, TRUE) * e2)

setMethod("*",
          signature("advector", "sparseMatrix"),
          function(e1, e2) e2 * e1 )

setMethod("/",
          signature("sparseMatrix", "advector"),
          function(e1, e2) magic(e1, TRUE) / e2)
