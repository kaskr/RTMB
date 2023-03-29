// RTMB configuration to include by TMB.h

#include <Rcpp.h>
// Any failed internal assertion sends Rcpp::exception
#define TMB_ABORT Rcpp::stop("TMB unexpected")
// Catch *all* std exceptions (not just bad_alloc)
#define TMB_CATCH catch(std::exception& excpt)
