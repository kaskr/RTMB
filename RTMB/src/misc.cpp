// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// Workarounds needed for parallel case:
// R_CheckStack() is junk when called from non-master thread.
// We can effectively disable the check by setting 'R_CStackDir = 0'.
#ifdef _OPENMP
extern int	R_CStackDir;
#endif
struct CStackWorkaround {
#ifdef _OPENMP
  int old;
  void begin() {
    old = R_CStackDir;
    if (omp_get_thread_num() != 0) R_CStackDir = 0;
  }
  void end() {
    R_CStackDir = old;
  }
#else
  void begin() { }
  void end() { }
#endif
};

namespace TMBad {
// Operator that evaluates an R function taking single numeric scalar as input
struct EvalOp : global::DynamicOperator< 1 , -1 > {
  static const bool is_linear = true;
  static const bool have_input_size_output_size = true; // FIXME: Should give compile time error if 'false'
  static const bool add_forward_replay_copy = true;
  Rcpp::Function F;
  size_t n;
  Index input_size()  const { return 1; }
  Index output_size() const { return n; }
  EvalOp (Rcpp::Function F, size_t n) : F(F), n(n) { }
  void forward(ForwardArgs<double> &args) {
#ifdef _OPENMP
#pragma omp critical
    {
#endif
      CStackWorkaround R;
      R.begin();
      Rcpp::NumericVector i = Rcpp::NumericVector::create(args.x(0));
      SEXP y = F(i);
      if ((size_t) LENGTH(y) != n) {
        R.end();
        Rcpp::stop("Wrong output length");
      }
      if (Rf_isReal(y)) {
        double* py = REAL(y);
        for (size_t i=0; i<n; i++) { args.y(i) = py[i]; }
      } else if (Rf_isInteger(y)) {
        int* py = INTEGER(y);
        for (size_t i=0; i<n; i++) { args.y(i) = py[i]; }
      } else {
        R.end();
        Rcpp::stop("EvalOp: Function must return 'real' or 'integer'");
      }
      R.end();
#ifdef _OPENMP
    }
#endif
  }
  template <class Type> void forward(ForwardArgs<Type> &args) {
    TMBAD_ASSERT(false);
  }
  template <class Type> void reverse(ReverseArgs<Type> &args) {
    // Void derivs
  }
  const char* op_name() {return "EvalOp";}
  void print(TMBad::global::print_config cfg) {
    Rcout << cfg.prefix;
    Rcout << "F=" << F << " ";
    Rcout << "n=" << n << "\n";
  }
};
}

// [[Rcpp::export]]
Rcpp::ComplexVector TapedEval(Rcpp::Function F, Rcpp::ComplexVector i) {
  if (!ad_context()) Rcpp::stop("TapedSubset requires an active ad context");
  if (!is_adscalar(i)) Rcpp::stop("TapedSubset requires ad scalar 'i'");
  ad i_ = ScalarInput(i);
  // Test eval to get n
  Rcpp::NumericVector i_test = Rcpp::NumericVector::create(i_.Value());
  Rcpp::NumericVector y_test = F(i_test);
  size_t n = LENGTH(y_test);
  // Add to tape
  std::vector<ad> x(1, i_);
  std::vector<ad> y = TMBad::global::Complete<TMBad::EvalOp>(F, n) (x);
  // Pass to R
  Rcpp::ComplexVector ans(n);
  for (size_t j=0; j < n; j++) {
    ans[j] = ad2cplx(y[j]);
  }
  DUPLICATE_ATTRIB(ans, y_test);
  return as_advector(ans);
}
