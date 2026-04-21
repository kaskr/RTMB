// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

/********************************************************************
 * Adding exponentially scaled versions missing in TMB
 ********************************************************************/
namespace atomic {
TMB_BIND_ATOMIC(bessel_k_expo,
		11,
		bessel_utils::bessel_k(x[0], x[1], 2.) )
TMB_BIND_ATOMIC(bessel_i_expo,
		11,
		bessel_utils::bessel_i(x[0], x[1], 2.) )
}
template<class Type>
Type besselK_expo(Type x, Type nu) {
  CppAD::vector<Type> args(3);
  args[0] = x;
  args[1] = nu;
  args[2] = 0; // derivative order
  return atomic::bessel_k_expo(args)[0];
}
template<class Type>
Type besselI_expo(Type x, Type nu) {
  CppAD::vector<Type> args(3);
  args[0] = x;
  args[1] = nu;
  args[2] = 0; // derivative order
  return atomic::bessel_i_expo(args)[0];
}

// [[Rcpp::export]]
ADrep distr_besselK (ADrep x, ADrep nu, bool expo=false) {
  int n1=x.size();
  int n2=nu.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = (expo ?
                                  besselK_expo(X1[i % n1], X2[i % n2]) :
                                  besselK     (X1[i % n1], X2[i % n2]));
  return ans;
}

// [[Rcpp::export]]
ADrep distr_besselI (ADrep x, ADrep nu, bool expo=false)
{
  int n1=x.size();
  int n2=nu.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = (expo ?
                                  besselI_expo(X1[i % n1], X2[i % n2]) :
                                  besselI     (X1[i % n1], X2[i % n2]));
  return ans;
}

// [[Rcpp::export]]
ADrep distr_besselJ (ADrep x, ADrep nu)
{
  int n1=x.size();
  int n2=nu.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = besselJ(X1[i % n1], X2[i % n2]);
  return ans;
}

// [[Rcpp::export]]
ADrep distr_besselY (ADrep x, ADrep nu)
{
  int n1=x.size();
  int n2=nu.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(x); const ad* X2 = adptr(nu);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = besselY(X1[i % n1], X2[i % n2]);
  return ans;
}
