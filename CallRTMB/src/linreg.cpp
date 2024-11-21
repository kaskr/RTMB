// Simple linear regression.
#include <Rcpp.h>
#define TMBAD_INDEX_TYPE uint64_t
#define TMB_LIB_INIT R_init_CallRTMB
#include <TMB.hpp>
#include <CallRTMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Get R function to be called via RTMB
  CallRTMB<Type> foo(getListElement(data, "foo"));
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));
  Type nll = -sum(dnorm(Y, foo(a+b*x), exp(logSigma), true));
  return nll;
}
