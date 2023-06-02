// [[Rcpp::depends(TMB)]]

/*
  ====================================================================

  This file implements _complex AD_.

  WARNING: Not to confuse with the RComplex representation of AD types

  ====================================================================
*/

#include "RTMB.h"

// Some extra interface helper functions
typedef std::complex<ad> adcplx;

// [[Rcpp::export]]
Rcpp::ComplexVector Arith2_complex(const Rcpp::ComplexVector &x,
                                   const Rcpp::ComplexVector &y,
                                   std::string op) {
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  size_t nx = x.size(), ny = y.size();
  size_t n = (std::min(nx, ny) > 0 ? std::max(nx, ny) : 0);
  Rcpp::ComplexVector z(n);
  // AD complex
  adcplx* X = (adcplx*) adptr(x);
  adcplx* Y = (adcplx*) adptr(y);
  adcplx* Z = (adcplx*) adptr(z);
  size_t ncplx = n / 2;
#define CALL(OP)                                                \
  for (size_t i=0; i<ncplx; i++) Z[i] = X[i % nx] OP Y[i % ny];
  if (!op.compare("+")) CALL(+)
  else if (!op.compare("-")) CALL(-)
  else if (!op.compare("*")) CALL(*)
  else if (!op.compare("/")) CALL(/)
  else if (!op.compare("^")) {
    for (size_t i=0; i<n; i++) Z[i] = pow(X[i % nx] , Y[i % ny]);
  }
  else Rcpp::stop("'%s' not implemented", op.c_str());
#undef CALL
  return as_advector(z);
}

// [[Rcpp::export]]
Rcpp::ComplexVector Math1_complex(const Rcpp::ComplexVector &x, std::string op) {
  CHECK_INPUT(x);
  size_t n = x.size();
  Rcpp::ComplexVector y(n);
  // AD complex
  adcplx* X = (adcplx*) adptr(x);
  adcplx* Y = (adcplx*) adptr(y);
  size_t ncplx = n / 2;
#define CALL(OP) for (size_t i=0; i<ncplx; i++) Y[i] = OP ( X[i] )
#define CUMC(OP) for (size_t i=1; i<ncplx; i++) Y[i] = Y[i-1] OP X[i];
  if (!op.compare("abs")) CALL(std::abs);
  else if (!op.compare("sqrt")) CALL(sqrt);
  else if (!op.compare("exp")) CALL(exp);
  else if (!op.compare("log")) CALL(log);
  else if (!op.compare("cos")) CALL(cos);
  else if (!op.compare("sin")) CALL(sin);
  else if (!op.compare("tan")) CALL(tan);
  else if (!op.compare("acos")) CALL(acos);
  else if (!op.compare("asin")) CALL(asin);
  else if (!op.compare("atan")) CALL(atan);
  else if (!op.compare("cosh")) CALL(cosh);
  else if (!op.compare("sinh")) CALL(sinh);
  else if (!op.compare("tanh")) CALL(tanh);
  // FIXME:
  // else if (!op.compare("acosh")) CALL(acosh);
  // else if (!op.compare("asinh")) CALL(asinh);
  // else if (!op.compare("atanh")) CALL(atanh);
  // else if (!op.compare("lgamma")) CALL(lgamma);
  // else if (!op.compare("gamma")) CALL(rtmb_gamma);
  else if (!op.compare("cumsum")) {
    if (n > 0) { Y[0] = X[0]; CUMC(+); }
  }
  else if (!op.compare("cumprod")) {
    if (n > 0) { Y[0] = X[0]; CUMC(*); }
  }
  else Rf_error("'%s' not implemented", op.c_str());
#undef CALL
#undef CUMC
  return as_advector(y);
}
