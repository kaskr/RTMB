// RTMB configuration to include by TMB.h

#include <Rcpp.h>
// Any failed internal assertion sends Rcpp::exception
#define TMB_ABORT Rcpp::stop("TMB unexpected")
// Catch *all* std exceptions (not just bad_alloc)
#define TMB_CATCH catch(std::exception& excpt)
// Do not include TMB's thread-safe workarounds
#ifdef _OPENMP
#define TMB_HAVE_THREAD_SAFE_R
#endif
// TMB FIXME: Some occurances of ASSERT and ASSERT2
#undef  ASSERT
#define ASSERT(x) TMBAD_ASSERT(x)
#undef  ASSERT2
#define ASSERT2(x, msg) TMBAD_ASSERT2(x, msg)
