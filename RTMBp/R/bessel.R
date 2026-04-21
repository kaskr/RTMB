##' @describeIn Distributions AD implementation of \link[base]{besselK}
##' @param expon.scaled See \link[base]{besselK}
setMethod("besselK",
          signature(x = "ad", nu = "ad", expon.scaled = "ANY"),
          function(x, nu, expon.scaled) {
            x <- advector (x)
            nu <- advector (nu)
            distr_besselK (x, nu, expon.scaled)
          })
##' @describeIn Distributions Default method
setMethod("besselK",
          signature(x = "num", nu = "num", expon.scaled = "ANY"),
          function(x, nu, expon.scaled) {
            base::besselK (x=x, nu=nu, expon.scaled)
          })
##' @describeIn Distributions AD implementation of \link[base]{besselI}
setMethod("besselI",
          signature(x = "ad", nu = "ad", expon.scaled = "ANY"),
          function(x, nu, expon.scaled) {
            x <- advector (x)
            nu <- advector (nu)
            distr_besselI (x, nu, expon.scaled)
          })
##' @describeIn Distributions Default method
setMethod("besselI",
          signature(x = "num", nu = "num", expon.scaled = "ANY"),
          function(x, nu, expon.scaled) {
            base::besselI (x=x, nu=nu, expon.scaled)
          })
##' @describeIn Distributions AD implementation of \link[base]{besselJ}
setMethod("besselJ",
          signature(x = "ad", nu = "ad"),
          function(x, nu) {
            x <- advector (x)
            nu <- advector (nu)
            distr_besselJ (x, nu)
          })
##' @describeIn Distributions Default method
setMethod("besselJ",
          signature(x = "num", nu = "num"),
          function(x, nu) {
            base::besselJ (x=x, nu=nu)
          })
##' @describeIn Distributions AD implementation of \link[base]{besselY}
setMethod("besselY",
          signature(x = "ad", nu = "ad"),
          function(x, nu) {
            x <- advector (x)
            nu <- advector (nu)
            distr_besselY (x, nu)
          })
##' @describeIn Distributions Default method
setMethod("besselY",
          signature(x = "num", nu = "num"),
          function(x, nu) {
            base::besselY (x=x, nu=nu)
          })
