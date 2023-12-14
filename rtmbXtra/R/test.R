test <- function() {
    RTMB::MakeTape(function(x)rtmbXtra:::SparseSquare(x*Matrix::.symDiagonal(10))@x,2)
}
