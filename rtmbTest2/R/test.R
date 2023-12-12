test <- function() {
    RTMB::MakeTape(function(x)rtmbTest2:::SparseSquare(x*Matrix::.symDiagonal(10))@x,2)
}
