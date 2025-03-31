// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// [[Rcpp::export]]
ADrep subset_ad(ADrep x,
                ADrep i) {
  return atomic::subset<ad>(x, i);
}

// [[Rcpp::export]]
ADrep findInterval_ad(ADrep x,
                      ADrep i) {
  return atomic::findInterval<ad>(x, i);
}

// [[Rcpp::export]]
ADrep order_ad(ADrep x) {
  return atomic::order<ad>(x);
}

// [[Rcpp::export]]
ADrep sort_ad(ADrep x) {
  return atomic::sort<ad>(x);
}
