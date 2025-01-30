// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// [[Rcpp::export]]
ADrep subset(ADrep x,
             ADrep i) {
  return atomic::subset<ad>(x, i);
}

// [[Rcpp::export]]
ADrep findInterval(ADrep x,
                   ADrep i) {
  return atomic::findInterval<ad>(x, i);
}

// [[Rcpp::export]]
ADrep order(ADrep x) {
  return atomic::order<ad>(x);
}

// [[Rcpp::export]]
ADrep sort(ADrep x) {
  return atomic::sort<ad>(x);
}
