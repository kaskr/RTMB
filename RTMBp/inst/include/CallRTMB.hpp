/*
  =================== Calling RTMB from a normal TMB model ==================

  DESCRIPTION: This file can be included in a normal TMB C++ template
  to allow calling R functions from the template. For this to work,
  the header of the TMB template should look something like:

  #include <Rcpp.h>
  #define TMBAD_INDEX_TYPE uint64_t
  #include <TMB.hpp>
  #include <CallRTMB.hpp>

  Inside the objective function you can receive functions from R like this:

  // Get function 'func' from .GlobalEnv
  CallRTMB<Type> F(Rcpp::Function("func"));
  // Get function 'func' from TMB data
  CallRTMB<Type> F(getListElement(data, "func"));
*/

// FIXME: Allow user to replace this include by "RTMB.h" and move definitions to other compile unit
#include "RTMB_stubs.cpp"

// Permanently change the global pointer of this DLL to that of RTMB.
void RTMB_pointerSet() {
  typedef SEXP(*FUN_PTR)(SEXP);
  FUN_PTR RTMB_getSetGlobalPtr =
    (FUN_PTR) R_FindSymbol("getSetGlobalPtr", "RTMB", NULL);
  if (RTMB_getSetGlobalPtr == NULL)
    Rcpp::stop("RTMB namespace must be loaded");
  SEXP rtmb_ptr = PROTECT(RTMB_getSetGlobalPtr(R_NilValue));
  getSetGlobalPtr(rtmb_ptr);
  UNPROTECT(1);
}
// In case of an R exception
void RTMB_exception_cleanup() {
  TMBad::global* glob = TMBad::get_glob();
  if (glob != NULL) {
    glob->ad_stop();
    delete glob;
  }
}
// Throw exception on illegal parallelization
void RTMB_check_valid_parallelization() {
    if (config.nthreads > 1 && config.autopar == 0)
      Rcpp::stop("Calling R from TMB with parallization requires 'autopar=TRUE'");
}

template <class Type>
struct CallRTMB_ { };
template<>
struct CallRTMB_<double> {
  // Construct function argument
  Rcpp::NumericVector rtmb_arg(const vector<double> &x) {
    return
      Rcpp::NumericVector(x.data(), x.data() + x.size());
  }
  // Construct function result
  vector<double> rtmb_res(Rcpp::RObject x_) {
    Rcpp::NumericVector x(x_);
    return
      Eigen::Map<Eigen::ArrayXd> ( (double*) x.begin(), x.size() );
  }
};
template<>
struct CallRTMB_<ad> {
  // Construct function argument
  ADrep rtmb_arg(const vector<ad> &x) {
    return ADrep(x.data(), x.data() + x.size());
  }
  // Construct function result
  vector<ad> rtmb_res(Rcpp::RObject x_) {
    ADrep x = ADrep(x_);
    return Eigen::Map<Eigen::Array<ad, -1, 1> > (x.adptr(), x.size());
  }
};
template <class Type>
struct CallRTMB : CallRTMB_<Type> {
  typedef CallRTMB_<Type> Base;
  typedef vector<Type> Vec;
  Rcpp::Function F;
  CallRTMB(Rcpp::Function F) : F(F) { }
#define EVALUATOR(...)                                          \
  Vec ans;                                                      \
  BEGIN_RCPP                                                    \
  try {                                                         \
    RTMB_check_valid_parallelization();                         \
    RTMB_pointerSet();                                          \
    Rcpp::RObject tmp = F(__VA_ARGS__);                         \
    ans = Base::rtmb_res(tmp);                                  \
  }                                                             \
  catch (...) {                                                 \
    RTMB_exception_cleanup();                                   \
    throw;                                                      \
  }                                                             \
  VOID_END_RCPP                                                 \
  return ans;
  // 1 argument
  Vec operator()(const Vec &x1) {
    EVALUATOR(Base::rtmb_arg(x1));
  }
  // 2 argument
  Vec operator()(const Vec &x1,
                 const Vec &x2) {
    EVALUATOR(Base::rtmb_arg(x1),
              Base::rtmb_arg(x2));
  }
  // 3 argument
  Vec operator()(const Vec &x1,
                 const Vec &x2,
                 const Vec &x3) {
    EVALUATOR(Base::rtmb_arg(x1),
              Base::rtmb_arg(x2),
              Base::rtmb_arg(x3));
  }
  // 4 argument
  Vec operator()(const Vec &x1,
                 const Vec &x2,
                 const Vec &x3,
                 const Vec &x4) {
    EVALUATOR(Base::rtmb_arg(x1),
              Base::rtmb_arg(x2),
              Base::rtmb_arg(x3),
              Base::rtmb_arg(x4));
  }
  #undef EVALUATOR
};
