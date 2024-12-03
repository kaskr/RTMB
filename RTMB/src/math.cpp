// [[Rcpp::depends(TMB)]]
// #include <Rcpp.h>
// #include "TMB.h"
#include "RTMB.h"
extern tape_config_t tape_config;

// [[Rcpp::export]]
ADrep Arith2(ADrep x,
             ADrep y,
             std::string op) {
  size_t nx = x.size(), ny = y.size();
  size_t n = (std::min(nx, ny) > 0 ? std::max(nx, ny) : 0);
  bool do_vectorize =
    tape_config.ops_vectorize() && (nx==1 || nx==n) && (ny==1 || ny==n);
  ADrep z(n);
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
  // Object determining attrib of result.
  // FIXME: Not quite accurate - check what R src does
  SEXP attrib_from = ( ny > nx || ny == 0 ? y : x );
  SHALLOW_DUPLICATE_ATTRIB(z, attrib_from);
  return z;
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
ADrep Math1(ADrep x, std::string op) {
  size_t n = x.size();
  bool do_vectorize =
    tape_config.math_vectorize() && n>1;
  ADrep y(n);
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
  else if (!op.compare("sign")) CALL(sign);
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
  else if (!op.compare("acosh")) CALL(acosh);
  else if (!op.compare("asinh")) CALL(asinh);
  else if (!op.compare("atanh")) CALL(atanh);
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
  SHALLOW_DUPLICATE_ATTRIB(y, x);
  return y;
}

// atan2 is not in any group !
// [[Rcpp::export]]
ADrep math_atan2 (ADrep y, ADrep x) {
  int n1=y.size();
  int n2=x.size();
  int nmax = std::max({n1, n2});
  int nmin = std::min({n1, n2});
  int n = (nmin == 0 ? 0 : nmax);
  ADrep ans(n);
  const ad* X1 = adptr(y); const ad* X2 = adptr(x);
  ad* Y = adptr(ans);
  for (int i=0; i<n; i++) Y[i] = TMBad::atan2(X1[i % n1], X2[i % n2]);
  return ans;
}

// [[Rcpp::export]]
ADrep Reduce1(ADrep x, std::string op) {
  size_t n = x.size();
  ADrep y(1);
  ad& ans = y.adptr()[0];
  ad* X = adptr(x);
#define REDUCE(OP) for (size_t i=0; i<n; i++) ans = ans OP X[i];
#define REDUCE2(OP) for (size_t i=1; i<n; i++) ans = OP(ans, X[i]);
  if (!op.compare("+")) {
    if ( !tape_config.sum_vectorize() ) {
      ans = 0.; REDUCE(+);
    } else {
      ans = TMBad::sum(TMBad::ad_segment(X, n));
    }
  } else if (!op.compare("*")) {
    ans = 1.; REDUCE(*);
  } else if (!op.compare("min")) {
    if (n == 0) Rcpp::stop("Length must be positive");
    ans = X[0]; REDUCE2(TMBad::min);
  } else if (!op.compare("max")) {
    if (n == 0) Rcpp::stop("Length must be positive");
    ans = X[0]; REDUCE2(TMBad::max);
  }
  else Rf_error("'%s' not implemented", op.c_str());
#undef REDUCE
#undef REDUCE2
  return y;
}

// [[Rcpp::export]]
ADrep matmul (ADrep x,
              ADrep y) {
  if (x.ncol() != y.nrow())
    Rcpp::stop("non-conformable arguments");
  ConstMapMatrix X = MatrixInput(x);
  ConstMapMatrix Y = MatrixInput(y);
  ADrep Z;
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
ADrep matinv (ADrep x) {
  if (x.ncol() != x.nrow())
    Rcpp::stop("Expected a square matrix");
  ConstMapMatrix X = MatrixInput(x);
  return MatrixOutput(atomic::matinv(matrix<ad>(X)));
}

template<class nlDensity>
ADrep colApply (ADrep x,
                nlDensity &F,
                bool give_log) {
  ConstMapMatrix X = MatrixInput(x);
  ADrep z(x.ncol()); ad* Z = adptr(z);
  for (int j=0; j < X.cols(); j++) {
    ad ans = -F(vector<ad>(X.col(j)));
    if (!give_log) ans = exp(ans);
    Z[j] = ans;
  }
  return z;
}
template<class nlDensity>
ADrep colApply2(ADrep x,
                ADrep keep,
                nlDensity &F,
                bool give_log) {
  ConstMapMatrix X = MatrixInput(x);
  ConstMapMatrix K = MatrixInput(keep);
  ADrep z(x.ncol()); ad* Z = adptr(z);
  for (int j=0; j < X.cols(); j++) {
    ad ans = -F(vector<ad>(X.col(j)), vector<ad>(K.col(j)));
    if (!give_log) ans = exp(ans);
    Z[j] = ans;
  }
  return z;
}

// [[Rcpp::export]]
ADrep dmvnorm0 (ADrep x,
                ADrep s,
                bool give_log,
                SEXP keep = R_NilValue) {
  if (s.ncol() != s.nrow())
    Rcpp::stop("cov matrix must be square");
  if (x.nrow() != s.nrow())
    Rcpp::stop("non-conformable arguments");
  ConstMapMatrix S = MatrixInput(s);
  auto nldens = density::MVNORM(matrix<ad>(S), tape_config.mvnorm_atomic() );
  if (Rf_isNull(keep)) {
    return colApply(x, nldens, give_log);
  } else {
    ADrep k (keep);
    if (x.size() != k.size())
      Rcpp::stop("x / keep non-conformable arguments");
    return colApply2(x, k, nldens, give_log);
  }
}

// [[Rcpp::export]]
ADrep dgmrf0 (ADrep x,
              Rcpp::RObject q,
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
  if (x.nrow() != (size_t) Dim[0])
    Rcpp::stop("non-conformable arguments");
  Eigen::SparseMatrix<ad> Q = SparseInput(q);
  auto nldens = density::GMRF(Q);
  return colApply(x, nldens, give_log);
}

// [[Rcpp::export]]
Rcpp::RObject SparseArith2(Rcpp::RObject x,
                           Rcpp::RObject y,
                           std::string op) {
  Rcpp::RObject z;
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
    else if (!op.compare("+")) z = MatrixOutput(X + Y);
    else if (!op.compare("-")) z = MatrixOutput(X - Y);
    else if (!op.compare("*")) z = SparseOutput(X.cwiseProduct(Y));
    else Rf_error("'%s' not implemented", op.c_str());
  }
  // Dense OP Sparse
  else if (is_admatrix(x) && is_adsparse(y)) {
    ConstMapMatrix          X = MatrixInput(x);
    Eigen::SparseMatrix<ad> Y = SparseInput(y);
    if (!op.compare("%*%")) z = MatrixOutput(X * Y);
    else if (!op.compare("+")) z = MatrixOutput(X + Y);
    else if (!op.compare("-")) z = MatrixOutput(X - Y);
    else if (!op.compare("*")) z = SparseOutput(Y.cwiseProduct(X));
    else Rf_error("'%s' not implemented", op.c_str());
  }
  else Rf_error("Wrong use of 'SparseArith2'");
  return z;
}

// [[Rcpp::export]]
Rcpp::RObject Dense2Sparse(ADrep x) {
  matrix<ad> X = MatrixInput(x);
  Eigen::SparseMatrix<ad> Y = asSparseMatrix(X);
  return SparseOutput(Y);
}

#define MATH_MATRIX_FUNCTION(MFUN)                      \
  ADrep math_ ## MFUN (Rcpp::RObject x) {               \
  matrix<ad> X;                                         \
  if (is_adsparse(x)) {                                 \
    X = SparseInput(x);                                 \
  } else if (is_admatrix(x)) {                          \
    X = MatrixInput(x);                                 \
  } else {                                              \
    Rcpp::stop( #MFUN ": Expected matrix-like input");  \
  }                                                     \
  if (X.rows() != X.cols())                             \
    Rcpp::stop( #MFUN ": Expected square matrix");      \
  return MatrixOutput(atomic::MFUN(X));                 \
}

// [[Rcpp::export]]
ADrep math_expm(Rcpp::RObject x);
MATH_MATRIX_FUNCTION(expm)

// [[Rcpp::export]]
ADrep math_sqrtm(Rcpp::RObject x);
MATH_MATRIX_FUNCTION(sqrtm)

// [[Rcpp::export]]
ADrep math_absm(Rcpp::RObject x);
MATH_MATRIX_FUNCTION(absm)

#undef MATH_MATRIX_FUNCTION

// [[Rcpp::export]]
ADrep expATv (Rcpp::RObject AT,
              ADrep v,
              ADrep N,
              Rcpp::List cfg,
              Rcpp::RObject orig) {
  if (!is_adsparse(AT)) Rcpp::stop("Expecting adsparse 'AT'");
  if (!is_adscalar(N)) Rcpp::stop("Expecting adscalar 'N'");
  // Inputs
  matrix<ad> v_ = MatrixInput(v);
  ad N_ = ScalarInput(N);
  // Configuration parameters
  sparse_matrix_exponential::config<ad> cfg_;
#define SET_CONFIG(XXX) if (!Rf_isNull(cfg[#XXX]))      \
  cfg_.XXX = Rcpp::IntegerVector((SEXP) cfg[#XXX])[0]
  SET_CONFIG(Nmax);
  SET_CONFIG(trace);
  SET_CONFIG(warn);
#undef SET_CONFIG
  // Cache functor
  typedef sparse_matrix_exponential::expm_series<ad> expm_t;
  expm_t* F;
  if (!orig.hasAttribute("SparseMatrixExponential")) {
    Eigen::SparseMatrix<ad> AT_ = SparseInput(AT);
    F = new expm_t (AT_, N_, cfg_);
    SEXP ptr = Rcpp::XPtr<expm_t>(F);
    orig.attr("SparseMatrixExponential") = ptr;
  } else {
    SEXP ptr = orig.attr("SparseMatrixExponential");
    F = Rcpp::XPtr<expm_t> (ptr);
  }
  // Evaluate
  matrix<ad> ans(v_.rows(), v_.cols());
  for (int j=0; j<ans.cols(); j++) {
    vector<ad> vec = v_.col(j).array();
    vector<ad> out = (*F)(vec);
    ans.col(j).array() = out;
  }
  return MatrixOutput(ans);
}

// [[Rcpp::export]]
ADrep SparseSolve(Rcpp::RObject s, ADrep x) {
  typedef Eigen::SparseLU<Eigen::SparseMatrix<ad>, Eigen::COLAMDOrdering<int> >
    SparseLU_t;
  SparseLU_t* solver;
  // Cache factorization as Matrix package does
  if (!Rcpp::RObject(s).hasAttribute("SparseLU")) {
    Eigen::SparseMatrix<ad> S = SparseInput(s);
    solver = new SparseLU_t;
    (*solver).compute(S);
    SEXP ptr = Rcpp::XPtr<SparseLU_t>(solver);
    Rcpp::RObject(s).attr("SparseLU") = ptr;
  } else {
    SEXP ptr = Rcpp::RObject(s).attr("SparseLU");
    solver = Rcpp::XPtr<SparseLU_t> (ptr);
  }
  ConstMapMatrix X = MatrixInput(x);
  matrix<ad> ans = (*solver).solve(X);
  return MatrixOutput(ans);
}
