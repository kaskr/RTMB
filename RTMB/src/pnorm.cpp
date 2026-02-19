// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

/********************************************************************
 * Adding 'log.p' argument missing in TMB
 ********************************************************************/

TMB_ATOMIC_VECTOR_FUNCTION(
                           // ATOMIC_NAME
                           log_pnorm_both
                           ,
                           // OUTPUT_DIM
                           2
                           ,
                           // ATOMIC_DOUBLE
                           Rf_pnorm_both(tx[0], &(ty[0]), &(ty[1]), 2, 1);
                           ,
                           // ATOMIC_REVERSE
                           Type f = dnorm(tx[0], Type(0), Type(1), true);
                           px[0] = py[0]*exp(f-ty[0]) - py[1]*exp(f-ty[1]);
                           )

// [[Rcpp::export]]
ADrep distr_log_pnorm(ADrep q, bool lower_tail = true) {
  int n = q.size();
  ADrep ans(n);
  const ad* X = adptr(q);
  ad* Y = adptr(ans);
  CppAD::vector<ad> tmp(1);
  int tail = !lower_tail;
  for (int i=0; i<n; i++) {
    tmp[0] = X[i];
    Y[i] = log_pnorm_both(tmp)[tail];
  }
  return ans;
}

// [[Rcpp::export]]
ADrep distr_pnorm (ADrep q, ADrep mean, ADrep sd) {
  int n1=q.size();
  int n2=mean.size();
  int n3=sd.size();
  int nmax = std::max({n1, n2, n3});
  int nmin = std::min({n1, n2, n3});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(q); const ad* X2 = adptr(mean); const ad* X3 = adptr(sd);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = pnorm(X1[i % n1], X2[i % n2], X3[i % n3]);
  if (n == n1) SHALLOW_DUPLICATE_ATTRIB(ans, q);
  return ans;
}
