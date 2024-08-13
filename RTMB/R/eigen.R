## Eigen function - general case
eigen_rescaled <- function(X) {
    X[] <- as.complex(X)
    n <- sqrt(length(X))
    dim(X) <- c(n,n)
    e <- eigen(X, symmetric=FALSE)
    D <- e[["values"]]
    V <- e[["vectors"]]
    ## Select locally unique version by scaling with a complex sign
    s <- colSums(V)
    s <- Mod(s) / s
    V <- V %*% diag(s)
    ## Rescale to get t(V)%*%V=I
    s <- sqrt(colSums(V*V))
    s <- 1 / s ##s <- Mod(s) / s
    V <- V %*% diag(s)
    ## output
    c(D, V)
}
eigen_rescaled_adj <- function(X, Y, dY) {
    n <- sqrt(length(X))
    i <- 1:n ## altrep
    j <- (n+1):(n+n*n) ## altrep
    D <- Y[i]; dD <- dY[i]
    V <- Y[j]; dV <- dY[j]
    dim(V) <- dim(dV) <- c(n,n)
    Vinv <- solve(V)
    Delta <- 1 / t( outer(D, D, "-") )
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
    n <- sqrt(length(x))
    D <- y[,1]
    U <- y[,-1, drop=FALSE]
    dD <- dy[,1]
    dU <- dy[,-1, drop=FALSE]
    Delta <- 1 / outer(D, D, "-")
    diag(Delta) <- 0
    ## Eigen value reverse update
    M1 <- U %*% diag(dD) %*% t(U)
    ## Eigen vector reverse update
    M2 <- U %*% ( (t(dU) %*% U ) * Delta ) %*% t(U)
    ## Total update
    M <- M1 + M2
    ## correction for lower triangle repr. of symmetric matrix
    M <- M+t(M); diag(M) <- diag(M)*.5; M[upper.tri(M)] <- 0
    M
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
              x <- as.matrix(x)
              if (symmetric) {
                  U <- upper.tri(x)
                  x[U] <- Conj(t(x)[U])
              }
              n <- nrow(x)              
              y <- eigen_rescaled_atomic(x)
              D <- y[1:n]
              V <- y[(n+1):(n+n*n)]
              dim(V) <- c(n,n)
              lngt <- sqrt(colSums(Re(V * Conj(V))))
              V <- V * adcomplex(rep(lngt, each=n))
              structure(list(values=D, vectors=V), class="eigen")
          })

##' @describeIn ADmatrix AD eigen decomposition for real matrices. The non-symmetric case is redirected to the \code{adcomplex} method. Note that argument \code{symmetric} is **not** auto-detected so **must** be specified.
##' @return List (vectors/values) with \code{advector} components in symmetric case and \code{adcomplex} components otherwise.
setMethod("eigen", "advector",
          function (x, symmetric, only.values, EISPACK) {
              if (symmetric) {
                  y <- eigen_symmetric_atomic(x)
                  structure(list(values=y[,1],
                                 vectors=y[,-1]), class="eigen")
              } else {
                  x <- adcomplex(x)
                  eigen(x, symmetric, only.values, EISPACK)
              }
          })
