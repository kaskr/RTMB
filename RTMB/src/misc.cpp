// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// Workarounds needed for parallel case:
// R_CheckStack() is junk when called from non-master thread.
// We can effectively disable the check by setting 'R_CStackDir = 0'.
#ifdef _OPENMP
extern int R_CStackDir;
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
/*
  EvalOp: Operator that evaluates an R function taking single numeric
  scalar as input.

  Memory management:

  - 'Rcpp::Function' automatically PROTECTs/UNPROTECTs. This is in
    general useful for temporary liftime objects, but we have to be
    careful when placing a 'Rcpp::Function' as a member of an object
    with 'unlimited lifetime' (e.g. an AD operator). If the object is
    copied (from R), it'll invoke PROTECTs that are not followed by
    the corresponding UNPROTECTs. A simple solution is to place
    'Rcpp::Function' in a shared pointer. That effectively disables
    the Rcpp memory management, while still taking advantage of the
    Rcpp cleanup code triggered when the last reference dies. I have
    verified that cleanup code is indeed triggered when the tapes
    (obj) *and* R function (F) are removed from the R
    workspace. Removing just F is not enough (nice).

  - Evaluating F from C++ is known to not work in parallel. We use
    'omp critical' to make sure F is never evaluated by two threads at
    the same time. However, that's not enough. R occasionally runs
    'R_StackCheck()' which works in a highly thead unsafe manner. It
    will fail for any other thread than the master. Inspection of
    'R_StackCheck()' reveals that it can be disabled by temporarily
    setting 'R_CStackDir = 0'. FIXME: We could probably work out the
    rest of these 'R_CStack...' variables in order for
    'R_StackCheck()' to work reliably for all threads. Don't know if
    this is worth while though.
*/
struct EvalOp : global::DynamicOperator< -1 , -1 > {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  std::shared_ptr<Rcpp::Function> Fptr;
  size_t m, n;
  Index input_size()  const { return m; }
  Index output_size() const { return n; }
  EvalOp (Rcpp::Function F, size_t m, size_t n) : Fptr(std::make_shared<Rcpp::Function>(F)), m(m), n(n) { }
  void forward(ForwardArgs<double> &args) {
#ifdef _OPENMP
#pragma omp critical
    {
#endif
      CStackWorkaround R;
      R.begin();
      Rcpp::NumericVector i(m);
      for (size_t l=0; l<m; l++) i[l] = args.x(l);
      SEXP y = (*Fptr)(i);
      // FIXME: Any Rcpp way of handling arbitrary output? For now doing PROTECT/UNPROTECT manually...
      PROTECT(y);
      if ((size_t) LENGTH(y) != n) {
        R.end();
	UNPROTECT(1); // y
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
	UNPROTECT(1); // y
        Rcpp::stop("EvalOp: Function must return 'real' or 'integer'");
      }
      R.end();
      UNPROTECT(1); // y
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
    Rcout << "F=" << *Fptr << " ";
    Rcout << "n=" << n << "\n";
  }
};
}

// [[Rcpp::export]]
Rcpp::ComplexVector TapedEval(Rcpp::Function F, Rcpp::ComplexVector i) {
  if (!ad_context()) Rcpp::stop("TapedSubset requires an active ad context");
  CHECK_INPUT(i);
  size_t m = i.size();
  ad* pi = adptr(i);
  // Test eval to get n
  Rcpp::NumericVector i_test(m);
  for (size_t l=0; l<m; l++) i_test[l] = pi[l].Value();
  Rcpp::NumericVector y_test = F(i_test);
  size_t n = LENGTH(y_test);
  // Add to tape
  std::vector<ad> x(pi, pi + m);
  std::vector<ad> y = TMBad::global::Complete<TMBad::EvalOp>(F, m, n) (x);
  // Pass to R
  Rcpp::ComplexVector ans(n);
  for (size_t j=0; j < n; j++) {
    ans[j] = ad2cplx(y[j]);
  }
  DUPLICATE_ATTRIB(ans, y_test);
  return as_advector(ans);
}
