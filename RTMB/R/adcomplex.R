setClass("adcomplex",
         slots=c(real="advector", imag="advector"))
adcomplex <- function(real, imag=rep(advector(0), length(real))) {
    real <- advector(real)
    new("adcomplex", real=real, imag=imag)
}
Re.adcomplex <- function(x) x@real
Im.adcomplex <- function(x) x@imag
Conj.adcomplex <- function(x) adcomplex(Re(x), -Im(x))
Mod.adcomplex <- function(x) sqrt(Re(x)*Re(x)+Im(x)*Im(x))
"+.adcomplex" <- function(x, y) {
    adcomplex(Re(x)+Re(y), Im(x)+Im(y))
}
"-.adcomplex" <- function(x, y) {
    if (missing(y)) {y <- x; x <- 0}
    adcomplex(Re(x)-Re(y), Im(x)-Im(y))
}
"*.adcomplex" <- function(x, y) {
    adcomplex(Re(x)*Re(y) - Im(x)*Im(y), Re(x)*Im(y) + Im(x)*Re(y))
}
recip <- function(x) adcomplex(1/Re(x * Conj(x))) * Conj(x)
"/.adcomplex" <- function(x, y) x * recip(y)
exp.adcomplex <- function(x) {
    s <- exp(Re(x))
    adcomplex(s*cos(Im(x)), s*sin(Im(x)))
}
sqrt.adcomplex <- function(x) {
    M <- Mod(x)
    s <- x/abs(x) ## FIXME: AD sign
    adcomplex(sqrt(.5*(Re(x)+M)), sqrt(.5*(-Re(x)+M)))
}
unsplit <- function(z) as.vector(t(cbind(Re(z), Im(z))))
resplit <- function(x) {
    dim(x) <- c(2, length(x)/2)
    adcomplex(x[1,], x[2,])
}
setMethod("fft", "advector",
          function(z, inverse) {
              ## if (!ad_context()) { ## Workaround
              ##     F <- .MakeTape(function(...)
              ##         fft(z,inverse),numeric(0))
              ##     return advector(F$eval(numeric(0)))
              ## }
              d <- dim(z)
              z <- unsplit(z)
              ans <- fft_complex(z, d)
              ans <- resplit(ans)
              dim(ans) <- d
              ans
          })
