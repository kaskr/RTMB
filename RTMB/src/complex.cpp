// [[Rcpp::depends(TMB)]]
#include <Rcpp.h>
#include "TMB.h"

// Dummy
template<class Type>
Type objective_function<Type>::operator() () {
  // pretend all parameters have been read
  this->index = this->theta.size();
  return 0;
}

typedef TMBad::ad_aug ad;
typedef std::vector<ad> ad_vec;

/* ========================================================================== */
/* ADFun object */
/* ========================================================================== */

// Control the tape (Start / Stop / Print)
void ad_start(TMBad::ADFun<>* adf) {
  adf->glob.ad_start();
}
void ad_stop(TMBad::ADFun<>* adf) {
  adf->glob.ad_stop();
}
void ad_print(TMBad::ADFun<>* adf) {
  adf->print();
}
// Some ADFun object evaluators
std::vector<double> Eval(TMBad::ADFun<>* tp, const std::vector<double> &x) {
  std::vector<double> y = (*tp)(x);
  return y;
}
Rcpp::ComplexVector EvalAD(TMBad::ADFun<>* tp, const Rcpp::ComplexVector &x);
Rcpp::NumericMatrix Jacobian(TMBad::ADFun<>* tp, const std::vector<double> &x) {
  std::vector<double> y = tp->Jacobian(x);
  Rcpp::NumericMatrix Jt(x.size(), y.size() / x.size(), y.begin());
  return Rcpp::transpose(Jt);
}
// Some ADFun object transformations
void JacFun(TMBad::ADFun<>* adf) {
  (*adf) = (*adf).JacFun();
}
void parallelize(TMBad::ADFun<>* adf, int nthreads) {
  (*adf) =  (*adf).parallelize(nthreads);
}
void fuse(TMBad::ADFun<>* adf) {
  (*adf).glob.set_fuse(true);
  (*adf).replay();
  (*adf).glob.set_fuse(false);
}
void optimize(TMBad::ADFun<>* adf) {
  (*adf).optimize();
}
SEXP ptrTMB(TMBad::ADFun<>* pf) {
  SEXP res;
  PROTECT(res=R_MakeExternalPtr((void*) pf,Rf_install("ADFun"),R_NilValue));
  //Rf_setAttrib(res,Rf_install("range.names"),info);
  SEXP ans;
  //Rf_setAttrib(res,Rf_install("par"),par);
  PROTECT(ans=ptrList(res));
  UNPROTECT(2);
  return ans;
}
void Copy(TMBad::ADFun<>* adf, Rcpp::XPtr<TMBad::ADFun<> > other) {
  *adf = *other;
}
// Collect free functions in a module
RCPP_MODULE(mod_adfun) {
  using namespace Rcpp;
  class_<TMBad::ADFun<> >( "adfun" )
  .constructor()
  .method("copy", &Copy)
  .method("start", &ad_start)
  .method("stop",  &ad_stop)
  .method("print", &ad_print)
  .method("eval",  &Eval)
  .method("evalAD",  &EvalAD)
  .method("jacobian", &Jacobian)
  .method("jacfun", &JacFun)
  .method("parallelize", &parallelize)
  .method("fuse", &fuse)
  .method("optimize", &optimize)
  .method("ptrTMB", &ptrTMB)
  ;
}

// Global Tape Configuration
static struct tape_config_t {
  int comparison; // Safe=0 / Taped=1 / Unsafe=2
  int atomic;     // No atomic=0 / Use atomic=1
  int vectorize;  // No vectorize =0 / Use vectorize=1
  tape_config_t() : comparison(0), atomic(1), vectorize(0) {}
  bool matmul_plain    () { return (atomic == 0); }
  bool matmul_atomic   () { return (atomic == 1) && vectorize==0; }
  bool matmul_TMBad    () { return (atomic == 1) && vectorize==1; }
  bool ops_vectorize   () { return vectorize == 1; }
  bool math_vectorize  () { return vectorize == 1; }
  bool sum_vectorize   () { return vectorize == 1; }
  bool compare_safe    () { return comparison == 0; }
  bool compare_taped   () { return comparison == 1; }
  bool compare_unsafe  () { return comparison == 2; }
  bool mvnorm_atomic   () { return (atomic == 1); }
} tape_config;
// [[Rcpp::export]]
Rcpp::List set_tape_config(int comparison=0, int atomic=1, int vectorize=0) {
  tape_config.comparison = comparison;
  tape_config.atomic = atomic;
  tape_config.vectorize = vectorize;
#define GET(name) Rcpp::Named(#name) = tape_config.name()
  return Rcpp::List::create(
                            GET(matmul_plain),
                            GET(matmul_atomic),
                            GET(matmul_TMBad),
                            GET(ops_vectorize),
                            GET(math_vectorize),
                            GET(sum_vectorize),
                            GET(compare_safe),
                            GET(compare_taped),
                            GET(compare_unsafe),
                            GET(mvnorm_atomic));
}

/* ========================================================================== */
/* AD vector object */
/* ========================================================================== */

Rcomplex ad2cplx(const ad &x) {
  static_assert(sizeof(ad) == sizeof(Rcomplex),
                "ad size must match Rcomplex");
  Rcomplex* px = (Rcomplex*)(&x);
  return *px;
}
ad cplx2ad(const Rcomplex &x) {
  static_assert(sizeof(ad) == sizeof(Rcomplex),
                "ad size must match Rcomplex");
  ad* px = (ad*)(&x);
  return *px;
}
ad* adptr(const Rcpp::ComplexVector &x) {
  static_assert(sizeof(ad) == sizeof(Rcomplex),
                "ad size must match Rcomplex");
  ad* px = (x.size() > 0 ? (ad*)(&(x[0])) : NULL);
  return px;
}
bool is_advector (SEXP x) {
  return Rf_inherits(x, "advector");
}
bool is_sparse (SEXP x) {
  return Rf_inherits(x, "adsparse");
}
bool is_scalar (SEXP x) {
  return is_advector(x) && (Rcpp::ComplexVector(x).size() == 1);
}
bool valid(const ad &x) {
  return
    !x.ontape() || x.in_context_stack(x.data.glob);
}
// [[Rcpp::export]]
bool valid(Rcpp::ComplexVector x) {
  for (int i=0; i<x.size(); i++)
    if (!valid(cplx2ad(x[i]))) return false;
  return true;
}
// [[Rcpp::export]]
bool ad_context() {
  return TMBad::get_glob() != NULL;
}

#define CHECK_INPUT(x)                                                  \
if (!is_advector(x))                                                    \
  Rcpp::stop("'" #x "' must be 'advector' (lost class attribute?)" );   \
 if (!valid(Rcpp::ComplexVector(x)))                                    \
  Rcpp::stop("'" #x "' is not a valid 'advector' (constructed using illegal operation?)" );

Rcpp::ComplexVector& as_advector(Rcpp::ComplexVector &x) {
  x.attr("class") = "advector";
  SET_S4_OBJECT(x);
  return x;
}

Rcpp::ComplexVector EvalAD(TMBad::ADFun<>* tp, const Rcpp::ComplexVector &x) {
  CHECK_INPUT(x);
  std::vector<ad> x_( (ad*) x.begin(), (ad*) x.end());
  std::vector<ad> y_ = (*tp)(x_);
  Rcpp::ComplexVector y( (Rcomplex*) (y_.data()), (Rcomplex*) (y_.data() + y_.size()) );
  return as_advector(y);
}

// [[Rcpp::export]]
Rcpp::ComplexVector advec(const Rcpp::NumericVector &x) {
  Rcpp::ComplexVector ans(x.size());
  for (int i=0; i<x.size(); i++) ans[i] = ad2cplx(ad(x[i]));
  return as_advector(ans);
}

// [[Rcpp::export]]
Rcpp::ComplexVector dependent(const Rcpp::ComplexVector &x) {
  CHECK_INPUT(x);
  if (TMBad::get_glob() == NULL)
    Rcpp::stop("No active AD context");
  Rcpp::ComplexVector ans(x.size());
  for (int i=0; i<x.size(); i++) {
    ad xad = cplx2ad(x[i]);
    xad.Dependent();
    ans[i] = ad2cplx(xad);
  }
  return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector independent(const Rcpp::ComplexVector &x) {
  CHECK_INPUT(x);
  if (TMBad::get_glob() == NULL)
    Rcpp::stop("No active AD context");
  Rcpp::ComplexVector ans(x.size());
  for (int i=0; i<x.size(); i++) {
    ad xad = cplx2ad(x[i]);
    xad.Independent();
    ans[i] = ad2cplx(xad);
  }
  return as_advector(ans);
}

/* Arith:

   • ‘"+"’, ‘"-"’, ‘"*"’, ‘"/"’, ‘"^"’, ‘"%%"’, ‘"%/%"’

   • ‘"&"’, ‘"|"’, ‘"!"’

   • ‘"=="’, ‘"!="’, ‘"<"’, ‘"<="’, ‘">="’, ‘">"’
*/

// Sparse case redirect
Rcpp::ComplexVector SparseArith2(const Rcpp::ComplexVector &x,
                                 const Rcpp::ComplexVector &y,
                                 std::string op);
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
  if (!op.compare("+")) CALL(+)
  else if (!op.compare("-")) CALL(-)
  else if (!op.compare("*")) CALL(*)
  else if (!op.compare("/")) CALL(/)
  else if (!op.compare("^")) {
    for (size_t i=0; i<n; i++) Z[i] = pow(X[i % nx] , Y[i % ny]);
  }
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
  else if (!op.compare("cos")) VCALL(cos)
  else if (!op.compare("sin")) VCALL(sin)
  else if (!op.compare("tan")) VCALL(tan)
  else if (!op.compare("acos")) VCALL(acos)
  else if (!op.compare("asin")) VCALL(asin)
  else if (!op.compare("atan")) VCALL(atan)
  else if (!op.compare("lgamma")) CALL(lgamma);
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
Rcpp::NumericVector getValues(const Rcpp::ComplexVector &x) {
  CHECK_INPUT(x);
  Rcpp::NumericVector ans(x.size());
  for (int i=0; i<x.size(); i++) {
    ans[i] = cplx2ad((x)[i]).Value() ;
  }
  return ans;
}

// [[Rcpp::export]]
Rcpp::LogicalVector getVariables(const Rcpp::ComplexVector &x) {
  CHECK_INPUT(x);
  Rcpp::LogicalVector ans(x.size());
  for (int i=0; i<x.size(); i++) {
    ans[i] = !cplx2ad((x)[i]).constant() ;
  }
  return ans;
}

// [[Rcpp::export]]
void dbgprint(const Rcpp::ComplexVector &x) {
  if (!is_advector(x)) Rcpp::stop("'x' must be advector");
  for (int i=0; i<x.size(); i++) {
    ad xi = cplx2ad((x)[i]) ;
    Rcout << "index="
          << xi.index()
          << " union={glob=" << xi.data.glob
          << ", value=" << xi.data.value << "}"
          << " valid=" << valid(xi) << "\n";
  }
}

// [[Rcpp::export]]
Rcpp::ComplexVector matmul (const Rcpp::ComplexMatrix &x,
                            const Rcpp::ComplexMatrix &y,
                            std::string method) {
  typedef Eigen::Map<Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > MapMatrix;
  typedef Eigen::Map<const Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > ConstMapMatrix;
  if (x.ncol() != y.nrow())
    Rcpp::stop("non-conformable arguments");
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  Rcpp::ComplexMatrix z(x.nrow(), y.ncol());
  ConstMapMatrix X((ad*) x.begin(), x.nrow(), x.ncol());
  ConstMapMatrix Y((ad*) y.begin(), y.nrow(), y.ncol());
  MapMatrix Z((ad*) z.begin(), z.nrow(), z.ncol());
  if ( tape_config.matmul_plain() )
    Z = X * Y;
  else if ( tape_config.matmul_atomic() )
    Z = atomic::matmul(matrix<ad>(X), matrix<ad>(Y));
  else if ( tape_config.matmul_TMBad() )
    Z = TMBad::matmul(matrix<ad>(X), matrix<ad>(Y));
  else
    Rf_error("Method '%s' not implemented", method.c_str());
  return as_advector(z);
}

template<class nlDensity>
Rcpp::ComplexVector colApply (const Rcpp::ComplexMatrix &x,
                              nlDensity &F,
                              bool give_log) {
  typedef Eigen::Map<const Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > ConstMapMatrix;
  ConstMapMatrix X((ad*) x.begin(), x.nrow(), x.ncol());
  Rcpp::ComplexVector z(x.ncol());
  for (int j=0; j < X.cols(); j++) {
    ad ans = -F(vector<ad>(X.col(j)));
    if (!give_log) ans = exp(ans);
    z[j] = ad2cplx(ans);
  }
  return as_advector(z);
}

// [[Rcpp::export]]
Rcpp::ComplexVector dmvnorm0 (const Rcpp::ComplexMatrix &x,
                              const Rcpp::ComplexMatrix &s,
                              bool give_log) {
  typedef Eigen::Map<const Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > ConstMapMatrix;
  if (s.ncol() != s.nrow())
    Rcpp::stop("cov matrix must be square");
  if (x.nrow() != s.nrow())
    Rcpp::stop("non-conformable arguments");
  CHECK_INPUT(x);
  CHECK_INPUT(s);
  ConstMapMatrix S((ad*) s.begin(), s.nrow(), s.ncol());
  auto nldens = density::MVNORM(matrix<ad>(S), tape_config.mvnorm_atomic() );
  return colApply(x, nldens, give_log);
}

Eigen::SparseMatrix<ad> SparseInput(Rcpp::S4 x);
// [[Rcpp::export]]
Rcpp::ComplexVector dgmrf0 (const Rcpp::ComplexMatrix &x,
                            Rcpp::S4 q,
                            bool give_log) {
  // TMB FIXME:
  //   newton::log_determinant<TMBad::global::ad_aug> (H=...) at /TMB/include/tmbutils/newton.hpp:1191
  // adds to tape!
  if (!ad_context())
    Rcpp::stop("'dgmrf0' currently requires an active tape");
  if (!is_sparse(q))
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

// ============================== Sparse matrices
Eigen::SparseMatrix<ad> SparseInput(Rcpp::S4 S) {
  Rcpp::ComplexVector x(S.slot("x"));
  CHECK_INPUT(x);
  Rcpp::IntegerVector i = S.slot("i");
  Rcpp::IntegerVector p = S.slot("p");
  Rcpp::IntegerVector Dim = S.slot("Dim");
  return
    Eigen::Map<const Eigen::SparseMatrix<ad> > (Dim[0], // rows()
                                                Dim[1], // cols()
                                                i.size(), // nonZeros()
                                                p.begin(), // outerIndexPtr()
                                                i.begin(), // innerIndexPtr()
                                                (ad*) x.begin(), // data()
                                                NULL); // innerNonZeroPtr();
}

Rcpp::S4 SparseOutput (const Eigen::SparseMatrix<ad> &S) {
  size_t nnz  = S.nonZeros();
  Rcpp::IntegerVector Dim(2);
  Dim[0] = S.rows();
  Dim[1] = S.cols();
  Rcpp::IntegerVector i(S.innerIndexPtr(), S.innerIndexPtr() + nnz);
  Rcpp::IntegerVector p(S.outerIndexPtr(), S.outerIndexPtr() + Dim[1] + 1);
  Rcpp::ComplexVector x( (Rcomplex*) (S.valuePtr()),
                         (Rcomplex*) (S.valuePtr() + nnz));
  Rcpp::S4 ans("adsparse");
  ans.slot("x") = as_advector(x);
  ans.slot("i") = i;
  ans.attr("p") = p;
  ans.attr("Dim") = Dim;
  return ans;
}

ad ScalarInput(SEXP x_) {
  Rcpp::ComplexVector x(x_);
  CHECK_INPUT(x);
  return cplx2ad(x[0]);
}

// [[Rcpp::export]]
Rcpp::S4 SparseArith2(SEXP x,
                      SEXP y,
                      std::string op) {
  Rcpp::S4 z;
  // Sparse OP Sparse
  if (is_sparse(x) && is_sparse(y)) {
    Eigen::SparseMatrix<ad> X = SparseInput(x);
    Eigen::SparseMatrix<ad> Y = SparseInput(y);
    if (!op.compare("+"))      z = SparseOutput(X + Y);
    else if (!op.compare("-")) z = SparseOutput(X - Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  // scalar OP Sparse
  else if (is_scalar(x) && is_sparse(y)) {
    ad X = ScalarInput(x);
    Eigen::SparseMatrix<ad> Y = SparseInput(y);
    if (!op.compare("*"))      z = SparseOutput(X * Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  // Sparse OP scalar
  else if (is_sparse(x) && is_scalar(y)) {
    Eigen::SparseMatrix<ad> X = SparseInput(x);
    ad Y = ScalarInput(y);
    if (!op.compare("*"))      z = SparseOutput(X * Y);
    else if (!op.compare("/"))      z = SparseOutput(X / Y);
    else Rf_error("'%s' not implemented", op.c_str());
  }
  else Rf_error("Wrong use of 'SparseArith2'");
  return z;
}
