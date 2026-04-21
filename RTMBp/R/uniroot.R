##' AD one-dimensional root finding.
##'
##' Univariate root finding extending R's native \link[stats]{uniroot} function to work in both *standard* and *AD* evaluation modes.
##' @rdname ADuniroot
##' @name ADuniroot
##' @aliases uniroot,ANY-method
##' @param f,interval,...,lower,upper,f.lower,f.upper,extendInt,check.conv,tol,maxiter,trace See \link[stats]{uniroot}.
##' @return List with component `"root"`.
setMethod("uniroot",
          c(f="ANY"),
          function(f, interval, ..., lower, upper, f.lower, f.upper, extendInt, check.conv, tol, maxiter, trace) {
            if (!ad_context()) {
              stats::uniroot(f=f, interval=interval, ..., lower=lower, upper=upper, f.lower=f.lower, f.upper=f.upper, extendInt=extendInt, check.conv=check.conv, tol=tol, maxiter=maxiter, trace=trace)
            } else {
              if (!missing(extendInt)) stop("Argument 'extendInt' not available in AD mode")
              ## Pass ...
              f <- dotfun(f, ...) ## See integrate.R
              ## Tape function
              F <- MakeTape(f, numeric(1))
              F$simplify()
              decompose_refs(.pointer(environment(F)$mod))
              F$reorder() ## Move x parameters up front
              parms <- resolve_refs(.pointer(environment(F)$mod))
              ptr <- .pointer(environment(F)$mod)
              ## Config
              cfg <- list(tol=as.double(tol), maxiter=as.integer(maxiter))
              ans <- list(root = ad_uniroot(ptr, lower, upper, parms, cfg))
              ans
            }
          })
