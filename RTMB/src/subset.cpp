// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// [[Rcpp::export]]
ADrep subset(ADrep x,
             ADrep i) {
  Eigen::Map<Eigen::Array<ad, -1, 1> > X(adptr(x), x.size());
  Eigen::Map<Eigen::Array<ad, -1, 1> > Y(adptr(i), i.size());
  vector<ad> ans = atomic::subset(vector<ad>(X), vector<ad>(Y));
  return ADrep(ans.data(), ans.data() + ans.size());
}
