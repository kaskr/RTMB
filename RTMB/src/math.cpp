// [[Rcpp::depends(TMB)]]
// #include <Rcpp.h>
// #include "TMB.h"
#include "RTMB.h"

// [[Rcpp::export]]
Rcpp::ComplexVector Arith2(const Rcpp::ComplexVector &x,
                           const Rcpp::ComplexVector &y,
                           std::string op) {
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  size_t nx = x.size(), ny = y.size();
  size_t n = (std::min(nx, ny) > 0 ? std::max(nx, ny) : 0);
  bool do_vectorize =
    tape_config.ops_vectorize() && (nx==1 || nx==n) && (ny==1 || ny==n);
  Rcpp::ComplexVector z(n);
  ad* X = adptr(x);
  ad* Y = adptr(y);
  ad* Z = adptr(z);
#define CALL(OP)                                                \
  {                                                             \
    if (!do_vectorize) {                                        \
      for (size_t i=0; i<n; i++) Z[i] = X[i % nx] OP Y[i % ny]; \
    } else {                                                    \
      TMBad::ad_segment S =                                     \
        TMBad::ad_segment(X, nx) OP TMBad::ad_segment(Y, ny);   \
      for (size_t i=0; i<n; i++) Z[i] = S[i];                   \
    }                                                           \
  }
#define COMPARISON(CONDEXP)                                             \
  {                                                                     \
    if (tape_config.compare_forbid()) {                                 \
      Rcpp::stop("Comparison is generally unsafe for AD types");        \
    } else if (tape_config.compare_taped()) {                           \
      for (size_t i=0; i<n; i++)                                        \
        Z[i] = CONDEXP(X[i % nx], Y[i % ny], ad(1.), ad(0.));           \
    } else if (tape_config.compare_allow()) {                           \
      Rcpp::stop("tape_config.compare_allow is not handled here");      \
    } else                                                              \
      Rcpp::stop("Nothing selected by tape_config.compare_* !");        \
  }
  if (!op.compare("+")) CALL(+)
  else if (!op.compare("-")) CALL(-)
  else if (!op.compare("*")) CALL(*)
  else if (!op.compare("/")) CALL(/)
  else if (!op.compare("^")) {
    for (size_t i=0; i<n; i++) Z[i] = pow(X[i % nx] , Y[i % ny]);
  }
  else if (!op.compare("==")) COMPARISON(CondExpEq)
  else if (!op.compare("!=")) COMPARISON(CondExpNe)
  else if (!op.compare(">=")) COMPARISON(CondExpGe)
  else if (!op.compare("<=")) COMPARISON(CondExpLe)
  else if (!op.compare(">"))  COMPARISON(CondExpGt)
  else if (!op.compare("<"))  COMPARISON(CondExpLt)
  else Rf_error("'%s' not implemented", op.c_str());
#undef CALL
  return as_advector(z);
}

/* Math:

   • ‘abs’, ‘sign’, ‘sqrt’,
   ‘floor’, ‘ceiling’, ‘trunc’,
   ‘round’, ‘signif’
   
   • ‘exp’, ‘log’, ‘expm1’, ‘log1p’,
   ‘cos’, ‘sin’, ‘tan’,
   ‘cospi’, ‘sinpi’, ‘tanpi’,
   ‘acos’, ‘asin’, ‘atan’

   ‘cosh’, ‘sinh’, ‘tanh’,
   ‘acosh’, ‘asinh’, ‘atanh’

   • ‘lgamma’, ‘gamma’, ‘digamma’, ‘trigamma’
   
   • ‘cumsum’, ‘cumprod’, ‘cummax’, ‘cummin’
*/
ad rtmb_gamma(ad x) { return exp(lgamma(x)); }

// [[Rcpp::export]]
Rcpp::ComplexVector Math1(const Rcpp::ComplexVector &x, std::string op) {
  CHECK_INPUT(x);
  size_t n = x.size();
  bool do_vectorize =
    tape_config.math_vectorize() && n>1;
  Rcpp::ComplexVector y(n);
  ad* X = adptr(x); // FIXME: TMBad::ad_segment(const *)
  ad* Y = adptr(y);
#define CALL(OP) for (size_t i=0; i<n; i++) Y[i] = OP ( X[i] )
#define VCALL(OP)                                               \
  {                                                             \
    if (!do_vectorize) { CALL(OP); }                            \
    else {                                                      \
      TMBad::ad_segment S = OP( TMBad::ad_segment(X, n) );      \
      for (size_t i=0; i<n; i++) Y[i] = S[i];                   \
    }                                                           \
  }
#define CUMC(OP) for (size_t i=1; i<n; i++) Y[i] = Y[i-1] OP X[i];
  if (!op.compare("abs")) CALL(fabs);
  else if (!op.compare("sqrt")) VCALL(sqrt)
  else if (!op.compare("exp")) VCALL(exp)
  else if (!op.compare("log")) VCALL(log)
  else if (!op.compare("expm1")) VCALL(expm1)
  else if (!op.compare("log1p")) VCALL(log1p)
  else if (!op.compare("cos")) VCALL(cos)
  else if (!op.compare("sin")) VCALL(sin)
  else if (!op.compare("tan")) VCALL(tan)
  else if (!op.compare("acos")) VCALL(acos)
  else if (!op.compare("asin")) VCALL(asin)
  else if (!op.compare("atan")) VCALL(atan)
  else if (!op.compare("cosh")) VCALL(cosh)
  else if (!op.compare("sinh")) VCALL(sinh)
  else if (!op.compare("tanh")) VCALL(tanh)
  // FIXME:
  // else if (!op.compare("acosh")) VCALL(acosh)
  // else if (!op.compare("asinh")) VCALL(asinh)
  // else if (!op.compare("atanh")) VCALL(atanh)
  else if (!op.compare("lgamma")) CALL(lgamma);
  else if (!op.compare("gamma")) CALL(rtmb_gamma);
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

// [[Rcpp::export]]
Rcpp::ComplexVector Reduce1(const Rcpp::ComplexVector &x, std::string op) {
  CHECK_INPUT(x);
  size_t n = x.size();
  Rcpp::ComplexVector y(1);
  ad ans = 0;
#define REDUCE(OP) for (size_t i=0; i<n; i++) ans = ans OP cplx2ad(x[i]);
  if (!op.compare("+")) {
    if ( !tape_config.sum_vectorize() ) {
      ans = 0.; REDUCE(+);
    } else {
      ad* X = adptr(x);
      ans = TMBad::sum(TMBad::ad_segment(X, n));
    }
  } else if (!op.compare("*")) {
    ans = 1.; REDUCE(*);
  }
  else Rf_error("'%s' not implemented", op.c_str());
#undef REDUCE
  y[0] = ad2cplx(ans);
  return as_advector(y);
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix matmul (const Rcpp::ComplexMatrix &x,
                            const Rcpp::ComplexMatrix &y) {
  if (x.ncol() != y.nrow())
    Rcpp::stop("non-conformable arguments");
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  ConstMapMatrix X = MatrixInput(x);
  ConstMapMatrix Y = MatrixInput(y);
  Rcpp::ComplexMatrix Z;
  if ( tape_config.matmul_plain() )
    Z = MatrixOutput(X * Y);
  else if ( tape_config.matmul_atomic() )
    Z = MatrixOutput(atomic::matmul(matrix<ad>(X), matrix<ad>(Y)));
  else if ( tape_config.matmul_TMBad() ) {
    if (!ad_context())
      Rcpp::stop("tape_config.matmul_TMBad() requires an active AD context");
    Z = MatrixOutput(TMBad::matmul(matrix<ad>(X), matrix<ad>(Y)));
  }
  else
    Rcpp::stop("Nothing selected by tape_config.matmul_* !");
  return Z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix matinv (const Rcpp::ComplexMatrix &x) {
  if (x.ncol() != x.nrow())
    Rcpp::stop("Expected a square matrix");
  CHECK_INPUT(x);
  ConstMapMatrix X = MatrixInput(x);
  return MatrixOutput(atomic::matinv(matrix<ad>(X)));
}

template<class nlDensity>
Rcpp::ComplexVector colApply (const Rcpp::ComplexMatrix &x,
                              nlDensity &F,
                              bool give_log) {
  ConstMapMatrix X((ad*) x.begin(), x.nrow(), x.ncol());
  Rcpp::ComplexVector z(x.ncol());
  for (int j=0; j < X.cols(); j++) {
    ad ans = -F(vector<ad>(X.col(j)));
    if (!give_log) ans = exp(ans);
    z[j] = ad2cplx(ans);
  }
  return as_advector(z);
}
template<class nlDensity>
Rcpp::ComplexVector colApply2(const Rcpp::ComplexMatrix &x,
                              const Rcpp::ComplexVector &keep,
                              nlDensity &F,
                              bool give_log) {
  ConstMapMatrix X((ad*) x.begin(), x.nrow(), x.ncol());
  ConstMapMatrix K((ad*) keep.begin(), x.nrow(), x.ncol());
  Rcpp::ComplexVector z(x.ncol());
  for (int j=0; j < X.cols(); j++) {
    ad ans = -F(vector<ad>(X.col(j)), vector<ad>(K.col(j)));
    if (!give_log) ans = exp(ans);
    z[j] = ad2cplx(ans);
  }
  return as_advector(z);
}

// [[Rcpp::export]]
Rcpp::ComplexVector dmvnorm0 (const Rcpp::ComplexMatrix &x,
                              const Rcpp::ComplexMatrix &s,
                              bool give_log,
                              SEXP keep = R_NilValue) {
  if (s.ncol() != s.nrow())
    Rcpp::stop("cov matrix must be square");
  if (x.nrow() != s.nrow())
    Rcpp::stop("non-conformable arguments");
  CHECK_INPUT(x);
  CHECK_INPUT(s);
  ConstMapMatrix S((ad*) s.begin(), s.nrow(), s.ncol());
  auto nldens = density::MVNORM(matrix<ad>(S), tape_config.mvnorm_atomic() );
  if (Rf_isNull(keep)) {
    return colApply(x, nldens, give_log);
  } else {
    Rcpp::ComplexVector k (keep);
    if (x.size() != k.size())
      Rcpp::stop("x / keep non-conformable arguments");
    CHECK_INPUT(k);
    return colApply2(x, k, nldens, give_log);
  }
}

// [[Rcpp::export]]
Rcpp::ComplexVector dgmrf0 (const Rcpp::ComplexMatrix &x,
                            Rcpp::S4 q,
                            bool give_log) {
  // TMB FIXME:
  //   newton::log_determinant<TMBad::global::ad_aug> (H=...) at /TMB/include/tmbutils/newton.hpp:1191
  // adds to tape!
  if (!ad_context())
    Rcpp::stop("'dgmrf0' currently requires an active tape");
  if (!is_adsparse(q))
    Rcpp::stop("Precision matrix must be sparse");
  Rcpp::IntegerVector Dim = q.slot("Dim");
  if (Dim[0] != Dim[1])
    Rcpp::stop("Precision matrix must be square");
  if (x.nrow() != Dim[0])
    Rcpp::stop("non-conformable arguments");
  CHECK_INPUT(x);
  CHECK_INPUT(q.slot("x"));
  Eigen::SparseMatrix<ad> Q = SparseInput(q);
  auto nldens = density::GMRF(Q);
  return colApply(x, nldens, give_log);
}

// [[Rcpp::export]]
SEXP SparseArith2(SEXP x,
                  SEXP y,
                  std::string op) {
  SEXP z;
  // Sparse OP Sparse
  if (is_adsparse(x) && is_adsparse(y)) {
    Eigen::SparseMatrix<ad> X = SparseInput(x);
    Eigen::SparseMatrix<ad> Y = SparseInput(y);
    if (!op.compare("+"))      z = SparseOutput(X + Y);
    else if (!op.compare("-")) z = SparseOutput(X - Y);
    else if (!op.compare("%*%")) z = SparseOutput(X * Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  // scalar OP Sparse
  else if (is_adscalar(x) && is_adsparse(y)) {
    ad X = ScalarInput(x);
    Eigen::SparseMatrix<ad> Y = SparseInput(y);
    if (!op.compare("*"))      z = SparseOutput(X * Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  // Sparse OP scalar
  else if (is_adsparse(x) && is_adscalar(y)) {
    Eigen::SparseMatrix<ad> X = SparseInput(x);
    ad Y = ScalarInput(y);
    if (!op.compare("*"))      z = SparseOutput(X * Y);
    else if (!op.compare("/"))      z = SparseOutput(X / Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  // Sparse OP Dense
  else if (is_adsparse(x) && is_admatrix(y)) {
    Eigen::SparseMatrix<ad> X = SparseInput(x);
    ConstMapMatrix          Y = MatrixInput(y);
    if (!op.compare("%*%")) z = MatrixOutput(X * Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  // Dense OP Sparse
  else if (is_admatrix(x) && is_adsparse(y)) {
    ConstMapMatrix          X = MatrixInput(x);
    Eigen::SparseMatrix<ad> Y = SparseInput(y);
    if (!op.compare("%*%")) z = MatrixOutput(X * Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  else Rf_error("Wrong use of 'SparseArith2'");
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix math_expm (SEXP x) {
  matrix<ad> X;
  if (is_adsparse(x)) {
    X = SparseInput(x);
  } else if (is_admatrix(x)) {
    X = MatrixInput(x);
  } else {
    Rcpp::stop("expm: Expected matrix-like input");
  }
  if (X.rows() != X.cols())
    Rcpp::stop("expm: Expected square matrix");
  return MatrixOutput(expm(X));
}
