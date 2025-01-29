// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// [[Rcpp::export]]
ADrep subset(ADrep x,
             ADrep i) {
  vector<ad> ans = atomic::subset<ad>(x, i);
  return ADrep(ans.data(), ans.data() + ans.size());
}

// [[Rcpp::export]]
ADrep findInterval(ADrep x,
                   ADrep i) {
  vector<ad> ans = atomic::findInterval<ad>(x, i);
  return ADrep(ans.data(), ans.data() + ans.size());
}
