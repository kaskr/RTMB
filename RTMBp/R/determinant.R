## c(log(abs(det(x))), sign(det(x)))
logdet_atomic <- ADjoint(
    function(x) {
        unlist(determinant(x, logarithm=TRUE))
    },
    function(x, y, dy) {
        t(solve(x)) * dy[1]
    },
    name = "logdet")

##' @describeIn ADmatrix AD log determinant
##' @param logarithm Not used
##' @examples
##' F <- MakeTape(det, diag(2)) ## Indirectly available via 'determinant'
##' F$jacobian(matrix(1:4,2))
determinant.advector <- function(x, logarithm = TRUE, ...) {
    if (!logarithm) stop("'logarithm' must be 'TRUE'")
    y <- logdet_atomic(x)
    return(structure(list(
        modulus = structure(y[1], logarithm = TRUE),
        sign = y[2]),
        class = "det"))
}
