## Eigen function
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

eigen_rescaled_atomic <- ADjoint(eigen_rescaled,
                                 eigen_rescaled_adj,
                                 complex=TRUE)

setMethod("eigen", "adcomplex",
          function (x, symmetric, only.values, EISPACK) {
              ## x <- as.matrix(x) ## FIXME: Not working for adcomplex
              if (symmetric) {
                  U <- upper.tri(Re(x)) ## FIXME: Not working for adcomplex
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

setMethod("eigen", "advector",
          function (x, symmetric, only.values, EISPACK) {
              x <- adcomplex(x)
              ans <- eigen(x, symmetric, only.values, EISPACK)
              if (symmetric) {
                  ans$values <- Re(ans$values)
                  ans$vectors <- Re(ans$vectors)
              }
              ans
          })
