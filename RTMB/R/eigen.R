## Eigen function - general case
eigen_rescaled <- function(X) {
    X[] <- as.complex(X)
    e <- eigen(X, symmetric=FALSE)
    D <- e[["values"]]
    V <- e[["vectors"]]
    ## Select locally unique version by scaling with a complex sign
    s <- colSums(V)
    s[s == 0i] <- 1
    s <- Mod(s) / s
    V <- V %*% diag(s)
    ## Rescale to get t(V)%*%V=I
    s <- sqrt(colSums(V*V))
    s <- 1 / s ##s <- Mod(s) / s
    V <- V %*% diag(s)
    ## output
    cbind(D, V)
}
eigen_rescaled_adj <- function(X, Y, dY) {
    D <- Y[,1]
    V <- Y[,-1, drop=FALSE]
    dD <- dY[,1]
    dV <- dY[,-1, drop=FALSE]
    Vinv <- solve(V)
    Diff <- t(outer(D, D, "-"))
    Delta <- Conj(Diff) / (Diff * Conj(Diff) + 1e-300)
    diag(Delta) <- 0
    a <- colSums(dV * V)
    dX_ <- diag(dD) + Delta * ( (t(V) %*% dV) - t( a * (t(V) %*% V) ) )
    t(Vinv) %*% dX_ %*% t(V)
}

## Eigen function - symmetric case
eigen_symmetric <- function(x) {
    e <- base::eigen(x, symmetric=TRUE)
    cbind(e$val, e$vec)
}
eigen_symmetric_adj <- function(x, y, dy) {
    D <- y[,1]
    U <- y[,-1, drop=FALSE]
    dD <- dy[,1]
    dU <- dy[,-1, drop=FALSE]
    Diff <- outer(D, D, "-")
    Delta <- sign(Diff) / (abs(Diff) + 1e-300)
    diag(Delta) <- 0
    S <- (t(dU) %*% U) * Delta + diag(dD)
    S <- S + t(S)
    S <- U %*% S %*% t(U)
    diag(S) <- diag(S) * 0.5
    S[upper.tri(S)] <- 0
    S
}

## - General complex
## - General real
## - Symmetric complex (Hermitian)
eigen_rescaled_atomic <- ADjoint(eigen_rescaled,
                                 eigen_rescaled_adj,
                                 complex=TRUE)
## - Symmetric real
eigen_symmetric_atomic <- ADjoint(eigen_symmetric,
                                  eigen_symmetric_adj,
                                  complex=FALSE)

##' @describeIn ADmatrix General AD eigen decomposition for complex matrices. Note that argument \code{symmetric} is **not** auto-detected so **must** be specified.
##' @return List (vectors/values) with \code{adcomplex} components.
##' @param symmetric Logical; Is input matrix symmetric (Hermitian) ?
##' @param only.values Ignored
##' @param EISPACK Ignored
setMethod("eigen", "adcomplex",
          function (x, symmetric, only.values, EISPACK) {
              if (missing(symmetric))
                  stop("RTMB does not auto detect argument 'symmetric'. Please specify")
              x <- as.matrix(x)
              if (symmetric) {
                  U <- upper.tri(x)
                  x[U] <- Conj(t(x)[U])
              }
              y <- eigen_rescaled_atomic(x)
              D <- y[,1]
              V <- y[,-1,drop=FALSE]
              lngt <- sqrt(colSums(Re(V * Conj(V))))
              scale <- 1 / lngt
              V <- V * rep(scale, each=nrow(V))
              structure(list(values=D, vectors=V), class="eigen")
          })

##' @describeIn ADmatrix AD eigen decomposition for real matrices. The non-symmetric case is redirected to the \code{adcomplex} method. Note that argument \code{symmetric} is **not** auto-detected so **must** be specified.
##' @return List (vectors/values) with \code{advector} components in symmetric case and \code{adcomplex} components otherwise.
setMethod("eigen", "advector",
          function (x, symmetric, only.values, EISPACK) {
              if (missing(symmetric))
                  stop("RTMB does not auto detect argument 'symmetric'. Please specify")
              if (symmetric) {
                  y <- eigen_symmetric_atomic(x)
                  structure(list(values=y[,1],
                                 vectors=y[,-1]), class="eigen")
              } else {
                  x <- adcomplex(x)
                  eigen(x, symmetric, only.values, EISPACK)
              }
          })

##' @describeIn ADmatrix AD svd decomposition for real matrices.
##' @param nu Ignored
##' @param nv Ignored
##' @param LINPACK Ignored
setMethod("svd", "advector",
          function (x, nu, nv, LINPACK = FALSE) {
              if (!missing(nu)) stop("'nu' argument not implemented")
              if (!missing(nv)) stop("'nv' argument not implemented")
              trans <- (nrow(x) > ncol(x))
              if (trans) x <- t(x)
              e <- eigen(x%*%t(x), symmetric=TRUE)
              u <- e$vec
              d <- sqrt(e$val)
              v <- t((1/d)*(t(u)%*%x))
              ans <- list(d=d, u=u, v=v)
              if (trans) ans[2:3] <- ans[3:2]
              ans
          })
