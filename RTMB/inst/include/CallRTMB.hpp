/*
  =================== Calling RTMB from a normal TMB model ==================

  DESCRIPTION: This file can be included in a normal TMB C++ template
  to allow calling R functions from the template. For this to work,
  the header of the TMB template should look something like:

  #include <Rcpp.h>
  #define TMBAD_INDEX_TYPE uint64_t
  #include <TMB.hpp>
  #include <CallRTMB.hpp>

  Inside the objective function you can retrieve functions from R like this:

  // Get function 'func' from .GlobalEnv
  CallRTMB<Type> F(Rcpp::Function("func"));
  // Get function 'func' from TMB data
  CallRTMB<Type> F(getListElement(data, "func"));
*/

// From RTMB API
typedef TMBad::ad_aug ad;
// From RTMB API
Rcomplex ad2cplx(const ad &x) {
  static_assert(std::is_same<TMBad::Index, uint64_t>::value,
                "Please add: #define TMBAD_INDEX_TYPE uint64_t");
  static_assert(sizeof(ad) == sizeof(Rcomplex),
                "ad size must match Rcomplex");
  Rcomplex* px = (Rcomplex*)(&x);
  return *px;
}
// From RTMB API
ad cplx2ad(const Rcomplex &x) {
  static_assert(std::is_same<TMBad::Index, uint64_t>::value,
                "Please add: #define TMBAD_INDEX_TYPE uint64_t");
  static_assert(sizeof(ad) == sizeof(Rcomplex),
                "ad size must match Rcomplex");
  ad* px = (ad*)(&x);
  return *px;
}
// From RTMB API
Rcpp::ComplexVector as_advector(Rcpp::ComplexVector x) {
  x = Rf_asS4(x, TRUE, FALSE);
  x.attr("class") = "advector";
  return x;
}
// Enable/disable RTMB replay to this DLL
void RTMB_pointerSwap() {
  typedef SEXP(*FUN_PTR)(SEXP);
  static FUN_PTR RTMB_getSetGlobalPtr =
    (FUN_PTR) R_FindSymbol("getSetGlobalPtr", "RTMB", NULL);
  SEXP      ptr = PROTECT(     getSetGlobalPtr(R_NilValue));
  SEXP rtmb_ptr = PROTECT(RTMB_getSetGlobalPtr(R_NilValue));
  getSetGlobalPtr(rtmb_ptr);
  RTMB_getSetGlobalPtr(ptr);
  UNPROTECT(2);
}
// In case of an R exception
void RTMB_exception_cleanup() {
  RTMB_pointerSwap();
  TMBad::global* glob = TMBad::get_glob();
  if (glob != NULL) {
    glob->ad_stop();
    delete glob;
  }
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
  vector<double> rtmb_res(const Rcpp::NumericVector &x) {
    return
      Eigen::Map<Eigen::ArrayXd> ( (double*) x.begin(), x.size() );
  }
};
template<>
struct CallRTMB_<ad> {
  // Construct function argument
  Rcpp::ComplexVector rtmb_arg(const vector<ad> &x) {
    Rcpp::ComplexVector ans(x.size());
    for (size_t i=0; i < (size_t) x.size(); i++) {
      ans[i] = ad2cplx(x[i]);
    }
    return as_advector(ans);
  }
  // Construct function result
  vector<ad> rtmb_res(const Rcpp::ComplexVector &x) {
    vector<ad> ans(x.size());
    for (size_t i=0; i < (size_t) x.size(); i++) {
      ans[i] = cplx2ad(x[i]);
    }
    return ans;
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
  RTMB_pointerSwap();                                           \
  try {                                                         \
    if (config.nthreads > 1)                                    \
      Rcpp::stop("Calling R from TMB requires nthreads=1");     \
    ans = Base::rtmb_res(F(__VA_ARGS__));                       \
  }                                                             \
  catch (...) {                                                 \
    RTMB_exception_cleanup();                                   \
    throw;                                                      \
  }                                                             \
  RTMB_pointerSwap();                                           \
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
