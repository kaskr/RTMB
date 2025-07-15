##' @describeIn Distributions Dirichlet. Calculate density.
ddirichlet <- function(x, alpha, log = FALSE) {
  K <- length(alpha)
  if (length(x) != K)
    stop("x[] and alpha[] must be equal length vectors.")
  r <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
  if (log)
    r
  else exp(r)
  if (log) r else exp(r)
}

## First we generate the version we want for AD types (dot signifies 'default argument')
##' @describeIn Distributions AD implementation of \link[brms]{ddirichlet}
setMethod("ddirichlet", signature("ad", "ad", "logical."),
          function(x, alpha, log) {
            # if (!is.matrix(x)) {
            #   x <- matrix(x, nrow = 1)
            # }
            # if (!is.matrix(alpha)) {
            #   alpha <- matrix(alpha, nrow(x), length(alpha), byrow = TRUE)
            # }
            # if (nrow(x) == 1L && nrow(alpha) > 1L) {
            #   x <- repl(x, nrow(alpha))
            #   x <- do_call(rbind, x)
            # } else if (nrow(x) > 1L && nrow(alpha) == 1L) {
            #   alpha <- repl(alpha, nrow(x))
            #   alpha <- do_call(rbind, alpha)
            # }
            # if (isTRUE(any(x < 0))) {
            #   stop("x must be non-negative.")
            # }
            # if (!isTRUE(all.equal(rowSums(x), rep(1, nrow(x))))) {
            #   stop("x must sum to 1 per row.")
            # }
            # if (isTRUE(any(alpha <= 0))) {
            #   stop("alpha must be positive.")
            # }
            # r <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) + rowSums((alpha - 1) * log(x))
            K <- length(alpha)
            if (length(x) != K)
              stop("x[] and alpha[] must be equal length vectors.")
            r <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
            if (log)
              r
            else exp(r)
          })

## This matches 'too much', so we fix by adding a specialization:
##' @describeIn Distributions Default method
setMethod("ddirichlet", signature("num", "num", "logical."),
          function(x, alpha, log) {
            K <- length(alpha)
            if (length(x) != K)
              stop("x[] and alpha[] must be equal length vectors.")
            r <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
            if (log)
              r
            else exp(r)
          })

## For S4 generics we add the OSA version like this:
##' @describeIn Distributions OSA implementation
setMethod("ddirichlet", "osa", function(x, alpha, log) {
  # prob <- prob / sum(prob)
  # if (is.null(size)) {
  #   size <- sum(x@x)
  # }
  ## Factorize in successive binomials
  perm <- order(attr(x@keep, "ord")) ## FIXME: Make extractor in osa.R ?
  x <- x[perm]
  alpha <- alpha[perm]
  ## Binomial parameters
  "c" <- ADoverload("c")
  cumsum0 <- function(x) c(0, cumsum(x[-length(x)]))
  size <- size - cumsum0(x@x)
  size <- getValues(advector(size)) ## Not variable
  prob <- prob / (1 - cumsum0(prob))
  sum(dbinom(x, size, prob, log = TRUE))
})

## For S4 generics we add the simref version like this:
##' @describeIn Distributions Simulation implementation. Modifies \code{x} and returns zero.
setMethod("ddirichlet", "simref", function(x, alpha, log) {
  nrep <- 1 ## ddirichlet is not vectorized
  K <- length(alpha)
  # x <- matrix(rgamma(K * nrep, alpha), ncol = K, byrow = TRUE)
  # sm <- tcrossprod(x, rep(1, K))
  # x/as.vector(sm)
  r <- rgamma(K * nrep, alpha)
  x[] <- r / sum(r)
  rep(0, nrep)
})

## To prevent unintendend usage, we change the default method.
InvalidMethod <- function(which=-2) {
  fr <- sys.frame(which)
  cl <- match.call(sys.function(which), sys.call(which), FALSE)
  cls <- sapply(names(cl[-1]), function(nm)class(get(nm, envir=fr)))
  cl[-1] <- cls
  c("Unexpected combination of classes used in AD context:\n", deparse(cl))
}

##' @describeIn Distributions Default implementation that checks for invalid usage.
setMethod("ddirichlet", signature("ANY", "ANY", "ANY"),
          function (x, alpha, log) {
            if (ad_context()) {
              stop(InvalidMethod())
            }
            stats::dmultinom(x=x, size=size, prob=prob, log=log)
          })

##' @describeIn Distributions Dirichlet-multinomial. Calculate density.
ddirmult <- function(x, alpha, log = FALSE) {
  K <- length(alpha)
  if (length(x) != K)
    stop("x[] and alpha[] must be equal length vectors.")
  sum_alpha <- sum(alpha)
  sum_x <- sum(x)
  r <- (lgamma(sum_x + 1) + sum(lgamma(x + alpha)) + lgamma(sum_alpha)) - 
    (sum(lgamma(x + 1)) + sum(lgamma(alpha)) + lgamma(sum_alpha + sum_x))
  if (log)
    r
  else exp(r)
  if (log) r else exp(r)
}

## First we generate the version we want for AD types (dot signifies 'default argument')
##' @describeIn Distributions AD implementation of ddirmult
setMethod("ddirmult", signature("ad", "ad", "logical."),
          function(x, alpha, log) {
            K <- length(alpha)
            if (length(x) != K)
              stop("x[] and alpha[] must be equal length vectors.")
            # if (!is.matrix(x)) {
            #   x <- matrix(x, nrow = 1)
            # }
            # if (!is.matrix(alpha)) {
            #   alpha <- matrix(alpha, nrow(x), length(alpha), byrow = TRUE)
            # }
            # sum_alpha <- rowSums(alpha)
            # sum_x <- rowSums(x)
            # r <- (lgamma(sum_x + 1) + rowSums(lgamma(x + alpha)) + lgamma(sum_alpha)) - 
            #   (rowSums(lgamma(x + 1)) + rowSums(lgamma(alpha)) + lgamma(sum_alpha + sum_x))
            sum_alpha <- sum(alpha)
            sum_x <- sum(x)
            r <- (lgamma(sum_x + 1) + sum(lgamma(x + alpha)) + lgamma(sum_alpha)) - 
              (sum(lgamma(x + 1)) + sum(lgamma(alpha)) + lgamma(sum_alpha + sum_x))
            if (log)
              r
            else exp(r)
          })

## This matches 'too much', so we fix by adding a specialization:
##' @describeIn Distributions Default method
setMethod("ddirmult", signature("num", "num", "logical."),
          function(x, alpha, log) {
            K <- length(alpha)
            if (length(x) != K)
              stop("x[] and alpha[] must be equal length vectors.")
            sum_alpha <- sum(alpha)
            sum_x <- sum(x)
            r <- (lgamma(sum_x + 1) + sum(lgamma(x + alpha)) + lgamma(sum_alpha)) - 
              (sum(lgamma(x + 1)) + sum(lgamma(alpha)) + lgamma(sum_alpha + sum_x))
            if (log)
              r
            else exp(r)
          })

## For S4 generics we add the OSA version like this:
##' @describeIn Distributions OSA implementation
setMethod("ddirmult", "osa", function(x, alpha, log) {
  # prob <- prob / sum(prob)
  # if (is.null(size)) {
  #   size <- sum(x@x)
  # }
  ## Factorize in successive binomials
  perm <- order(attr(x@keep, "ord")) ## FIXME: Make extractor in osa.R ?
  x <- x[perm]
  alpha <- alpha[perm]
  ## Binomial parameters
  "c" <- ADoverload("c")
  cumsum0 <- function(x) c(0, cumsum(x[-length(x)]))
  size <- size - cumsum0(x@x)
  size <- getValues(advector(size)) ## Not variable
  prob <- prob / (1 - cumsum0(prob))
  sum(dbinom(x, size, prob, log = TRUE))
})

## For S4 generics we add the simref version like this:
##' @describeIn Distributions Simulation implementation. Modifies \code{x} and returns zero.
setMethod("ddirmult", "simref", function(x, alpha, log) {
  nrep <- 1 ## ddirmult is not vectorized
  K <- length(alpha)
  r <- rgamma(K * nrep, alpha)
  prob <- r / sum(r)
  x[] <- stats::rmultinom(n = nrep, size = sum(alpha), prob = prob)
  rep(0, nrep)
})

## To prevent unintendend usage, we change the default method.
InvalidMethod <- function(which=-2) {
  fr <- sys.frame(which)
  cl <- match.call(sys.function(which), sys.call(which), FALSE)
  cls <- sapply(names(cl[-1]), function(nm)class(get(nm, envir=fr)))
  cl[-1] <- cls
  c("Unexpected combination of classes used in AD context:\n", deparse(cl))
}
##' @describeIn Distributions Default implementation that checks for invalid usage.
setMethod("ddirmult", signature("ANY", "ANY", "ANY"),
          function (x, alpha, log) {
            if (ad_context()) {
              stop(InvalidMethod())
            }
            stats::dmultinom(x=x, size=size, prob=prob, log=log)
          })
