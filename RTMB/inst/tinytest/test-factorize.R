## Test matrix factorizations

## From example(eigen)
base::eigen(cbind(c(0, 1i), c(-1i, 0)), symmetric=TRUE)

m <- structure(c(-0.626, 0.184, -0.836, 1.595, 0.33, -0.82, 0.487,
                 0.738, 0.576, -0.305, 1.512, 0.39, -0.621, -2.215,
                 1.125, -0.045, -0.016, 0.944, 0.821, 0.594, 0.919,
                 0.782, 0.075, -1.989, 0.62 ),
               dim = c(5L, 5L))

x <- m+(1-diag(5))*1i+diag(5)*2
base::eigen(x, symmetric=TRUE) ## Not actually hermitian => should access lower triangle only
