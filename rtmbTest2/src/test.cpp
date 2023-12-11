// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// [[Rcpp::export]]
SEXP SparseSquare(SEXP x) {
  Eigen::SparseMatrix<ad> X = SparseInput(x);
  Eigen::SparseMatrix<ad> Y = X * X;
  return SparseOutput(Y);
}
