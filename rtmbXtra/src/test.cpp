// [[Rcpp::depends(TMB)]]
#include <RTMB.h>

namespace adaptive {
template<class T>
T logspace_gamma(const T &x) {
  /* Tradeoff: The smaller x the better approximation *but* the higher
     risk of psigamma() overflow.
     > identical(lgamma(exp(-150)), 150)
     [1] TRUE
  */
  if (x < -150)
    return -x;
  else
    return lgamma(exp(x));
}
}
TMB_BIND_ATOMIC(logspace_gamma, 1, adaptive::logspace_gamma(x[0]))
template<class Type>
Type logspace_gamma(Type x) {
  CppAD::vector<Type> args(2); // Last index reserved for derivative order
  args[0] = x;
  args[1] = 0;
  return logspace_gamma(args)[0];
}

/* Based on examples in RTMB/src/distributions.cpp */
#define VECTORIZE_UNARY(FUN)                    \
size_t n = x.size();                            \
const ad* X = adptr(x);                         \
ADrep ans(n);                                   \
ad* Y = adptr(ans);                             \
for (size_t i=0; i < n; i++) {                  \
  Y[i] = FUN(X[i]);                             \
}                                               \
return ans;


// [[Rcpp::export]]
ADrep logspace_gamma(ADrep x) {
  VECTORIZE_UNARY(logspace_gamma);
}

extern "C" {
  /* See 'R-API: entry points to C-code' (Writing R-extensions) */
  double Rf_logspace_sub (double logx, double logy);
  void   Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
}
/* y(x) = logit_invcloglog(x) := log( exp(exp(x)) - 1 ) = logspace_sub( exp(x), 0 )

   y'(x) = exp(x) + exp(x-y) = exp( logspace_add(x, x-y) )

*/
TMB_ATOMIC_VECTOR_FUNCTION(
                           // ATOMIC_NAME
                           logit_invcloglog
                           ,
                           // OUTPUT_DIM
                           1,
                           // ATOMIC_DOUBLE
                           ty[0] = Rf_logspace_sub(exp(tx[0]), 0.);
                           ,
                           // ATOMIC_REVERSE
                           px[0] = exp( logspace_add(tx[0], tx[0]-ty[0]) ) * py[0];
                           )
template<class Type>
Type logit_invcloglog(Type x) {
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return logit_invcloglog(tx)[0];
}

// [[Rcpp::export]]
ADrep logit_invcloglog(ADrep x) {
  VECTORIZE_UNARY(logit_invcloglog);
}

/* y(x) = logit_pnorm(x) := logit( pnorm(x) ) =
   pnorm(x, lower.tail=TRUE,  log.p=TRUE) -
   pnorm(x, lower.tail=FALSE, log.p=TRUE)

   y'(x) = dnorm(x) * ( (1+exp(y)) + (1+exp(-y)) )

*/
double logit_pnorm(double x) {
  double log_p_lower, log_p_upper;
  Rf_pnorm_both(x, &log_p_lower, &log_p_upper, 2 /* both tails */, 1 /* log_p */);
  return log_p_lower - log_p_upper;
}
TMB_ATOMIC_VECTOR_FUNCTION(
                           // ATOMIC_NAME
                           logit_pnorm
                           ,
                           // OUTPUT_DIM
                           1,
                           // ATOMIC_DOUBLE
                           ty[0] = logit_pnorm(tx[0])
                           ,
                           // ATOMIC_REVERSE
                           Type zero = 0;
                           Type tmp1 = logspace_add(zero, ty[0]);
                           Type tmp2 = logspace_add(zero, -ty[0]);
                           Type tmp3 = logspace_add(tmp1, tmp2);
                           Type tmp4 = dnorm(tx[0], Type(0), Type(1), true) + tmp3;
                           px[0] = exp( tmp4 ) * py[0];
                           )
template<class Type>
Type logit_pnorm(Type x) {
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return logit_pnorm(tx)[0];
}

// [[Rcpp::export]]
ADrep logit_pnorm(ADrep x) {
  VECTORIZE_UNARY(logit_pnorm);
}
