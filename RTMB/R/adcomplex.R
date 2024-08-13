################################################################################
## This file contains:
## - Complex number arithmetic
## - Complex FFT
################################################################################

##' AD complex numbers
##'
##' A limited set of complex number operations can be used when constructing AD tapes. The available methods are listed in this help page.
##'
##' @rdname ADcomplex
##' @name ADcomplex
##' @examples
##' ## Tape using complex operations
##' F <- MakeTape(function(x) {
##'   x <- as.complex(x)
##'   y <- exp( x * ( 1 + 2i ) )
##'   c(Re(y), Im(y))
##' }, numeric(1))
##' F
##' F(1)
##' ## Complex FFT on the tape
##' G <- MakeTape(function(x) sum(Re(fft(x))), numeric(3))
##' G$simplify()
##' G$print()
NULL

setClass("adcomplex",
         slots=c(real="advector", imag="advector"))
##' @describeIn ADcomplex Construct \code{adcomplex} vector
##' @param real Real part
##' @param imag Imaginary part
##' @return Object of class \code{"adcomplex"}.
adcomplex <- function(real, imag=rep(advector(0), length(real))) {
    real <- advector(real)
    dim(imag) <- dim(real)
    new("adcomplex", real=real, imag=imag)
}
##' @describeIn ADcomplex As \link[base]{complex}
##' @param x An object of class \code{'adcomplex'}
##' @param y An object of class \code{'adcomplex'}
##' @param z An object of class \code{'adcomplex'}
##' @param object An object of class \code{'adcomplex'}
Re.adcomplex <- function(z) z@real
##' @describeIn ADcomplex As \link[base]{complex}
Im.adcomplex <- function(z) z@imag
##' @describeIn ADcomplex Print method
setMethod("show", "adcomplex", function (object)  {
    cat("class='adcomplex'\n")
    y <- complex(real=getValues(Re(object)), imaginary=getValues(Im(object)))
    dim(y) <- dim(object)
    print(y)
})
##' @describeIn ADcomplex As \link[base]{dim}
dim.adcomplex <- function(x) dim(Re(x))
##' @describeIn ADcomplex As \link[base]{dim}
##' @param value Replacement value
"dim<-.adcomplex" <- function(x, value) { dim(x@real) <- dim(x@imag) <- value; x }
##' @describeIn ADcomplex As \link[base]{[}
##' @param ... As \link[base]{[}
"[.adcomplex" <- function(x, ...) {
    adcomplex( Re(x)[...], Im(x)[...] )
}
##' @describeIn ADcomplex As \link[base]{[<-}
"[<-.adcomplex" <- function(x, ..., value) {
    x@real[...] <- Re(value)
    x@imag[...] <- Im(value)
    x
}
##' @describeIn ADcomplex As \link[base]{t}
t.adcomplex <- function(x) adcomplex(t(Re(x)), t(Im(x)))
##' @describeIn ADcomplex As \link[base]{length}
length.adcomplex <- function(x) length(Re(x))
##' @describeIn ADcomplex As \link[base]{complex}
Conj.adcomplex <- function(z) adcomplex(Re(z), -Im(z))
##' @describeIn ADcomplex As \link[base]{complex}
Mod.adcomplex <- function(z) sqrt(Re(z)*Re(z)+Im(z)*Im(z))
##' @describeIn ADcomplex As \link[base]{complex}
"+.adcomplex" <- function(x, y) {
    adcomplex(Re(x)+Re(y), Im(x)+Im(y))
}
##' @describeIn ADcomplex As \link[base]{complex}
"-.adcomplex" <- function(x, y) {
    if (missing(y)) {y <- x; x <- 0}
    adcomplex(Re(x)-Re(y), Im(x)-Im(y))
}
##' @describeIn ADcomplex As \link[base]{complex}
"*.adcomplex" <- function(x, y) {
    adcomplex(Re(x)*Re(y) - Im(x)*Im(y), Re(x)*Im(y) + Im(x)*Re(y))
}
recip <- function(x) adcomplex(1/Re(x * Conj(x))) * Conj(x)
##' @describeIn ADcomplex As \link[base]{complex}
"/.adcomplex" <- function(x, y) x * recip(y)
##' @describeIn ADcomplex As \link[base]{complex}
exp.adcomplex <- function(x) {
    s <- exp(Re(x))
    adcomplex(s*cos(Im(x)), s*sin(Im(x)))
}
##' @describeIn ADcomplex As \link[base]{complex}
sqrt.adcomplex <- function(x) {
    M <- Mod(x)
    s <- sign(Im(x))
    adcomplex(sqrt(.5*(Re(x)+M)), s * sqrt(.5*(-Re(x)+M)))
}
unsplit <- function(z) {
    dim(z) <- NULL
    as.vector(t(cbind(Re(z), Im(z))))
}
resplit <- function(x) {
    dim(x) <- c(2, length(x)/2)
    if (inherits(x, "advector")) {
        adcomplex(real=x[1,], imag=x[2,])
    } else {
        complex(real=x[1,], imaginary=x[2,])
    }
}
##' @describeIn ADcomplex Fast Fourier Transform equivalent to \link[stats]{fft}. Notably this is the **multivariate** transform when `x` is an array.
##' @param inverse As \link[stats]{fft}
setMethod("fft", "adcomplex",
          function(z, inverse) {
              if (!ad_context()) {
                  stop("No active AD context")
              }
              d <- dim(z)
              if (is.null(d)) d0 <- length(z) else d0 <- d
              z <- unsplit(z)
              ans <- fft_complex(z, d0, inverse)
              ans <- resplit(ans)
              dim(ans) <- d
              ans
          })
##' @describeIn ADcomplex If real input is supplied it is first converted to complex.
setMethod("fft", "advector",
          function(z, inverse) {
              fft(adcomplex(z), inverse)
          })

#####################################

##' @describeIn ADcomplex As \link[base]{rep}
rep.adcomplex <- function(x,...)
    adcomplex(rep(Re(x),...),
              rep(Im(x),...))

##' @describeIn ADcomplex Apply for each of real/imag
as.vector.adcomplex <- function(x, mode="any")
    adcomplex(as.vector(Re(x), mode),
              as.vector(Im(x), mode))

##' @describeIn ADcomplex Apply for real
is.matrix.adcomplex <- function(x) is.matrix(Re(x))

##' @describeIn ADcomplex Apply for each of real/imag
as.matrix.adcomplex <- function(x, ...) adcomplex(as.matrix(Re(x)),
                                                  as.matrix(Im(x)))

##' @describeIn ADcomplex Complex matrix multiply
setMethod("%*%", "adcomplex", function(x, y) {
    adcomplex(Re(x)%*%Re(y) - Im(x)%*%Im(y),
              Re(x)%*%Im(y) + Im(x)%*%Re(y))
})

##' @describeIn ADcomplex Complex matrix inversion and solve
setMethod("solve", "adcomplex", function(a, b) {
    A <- Re(a); B <- Im(a)
    Ainv <- solve(A)
    Y1 <- solve(A + B %*% Ainv %*% B)
    Y2 <- -Ainv %*% B %*% Y1
    ans <- adcomplex(Y1, Y2)
    if (!missing(b)) {
        ## b <- as.matrix(b)
        ans <- ans %*% b
    }
    ans
})

##' @describeIn ADcomplex Apply for each of real/imag
setMethod("colSums", "adcomplex",
          function(x) adcomplex(colSums(Re(x)),
                                colSums(Im(x))))

##' @describeIn ADcomplex Apply for each of real/imag
setMethod("rowSums", "adcomplex",
          function(x) adcomplex(rowSums(Re(x)),
                                rowSums(Im(x))))

##' @describeIn ADcomplex Apply for each of real/imag
setMethod("diag", "adcomplex",
          function(x) adcomplex(diag(Re(x)),
                                diag(Im(x))))
