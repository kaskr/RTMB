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
// Collect free functions in a module
RCPP_MODULE(mod_adfun) {
  using namespace Rcpp;
  class_<TMBad::ADFun<> >( "adfun" )
  .constructor()
  .method("start", &ad_start)
  .method("stop",  &ad_stop)
  .method("print", &ad_print)
  .method("eval",  &Eval)
  .method("jacobian", &Jacobian)
  .method("jacfun", &JacFun)
  .method("parallelize", &parallelize)
  .method("fuse", &fuse)
  .method("optimize", &optimize)
  .method("ptrTMB", &ptrTMB)
  ;
}

/* ========================================================================== */
/* AD vector object */
/* ========================================================================== */

typedef TMBad::ad_aug ad;
typedef std::vector<ad> ad_vec;

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
bool is_advector (const Rcpp::ComplexVector &x) {
  return
    x.hasAttribute("class") &&
    std::strcmp(x.attr("class"), "advector") == 0;
}
bool valid(const ad &x) {
  return
    !x.ontape() || x.in_context_stack(x.data.glob);
}
// [[Rcpp::export]]
bool valid(const Rcpp::ComplexVector &x) {
  for (int i=0; i<x.size(); i++)
    if (!valid(cplx2ad(x[i]))) return false;
  return true;
}

#define CHECK_INPUT(x)                                                  \
if (!is_advector(x))                                                    \
  Rcpp::stop("'" #x "' must be 'advector' (lost class attribute?)" );   \
if (!valid(x))                                                          \
  Rcpp::stop("'" #x "' is not a valid 'advector' (constructed using illegal operation?)" );

Rcpp::ComplexVector& as_advector(Rcpp::ComplexVector &x) {
  x.attr("class") = "advector";
  return x;
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

// [[Rcpp::export]]
Rcpp::ComplexVector Arith2(const Rcpp::ComplexVector &x,
                           const Rcpp::ComplexVector &y,
                           std::string op) {
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  size_t nx = x.size(), ny = y.size();
  size_t n = (std::min(nx, ny) > 0 ? std::max(nx, ny) : 0);
  Rcpp::ComplexVector z(n);
  const ad* X = adptr(x);
  const ad* Y = adptr(y);
  ad* Z = adptr(z);
#define CALL(OP) for (size_t i=0; i<n; i++) Z[i] = X[i % nx] OP Y[i % ny]
  if (!op.compare("+")) CALL(+);
  else if (!op.compare("-")) CALL(-);
  else if (!op.compare("*")) CALL(*);
  else if (!op.compare("/")) CALL(/);
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
  Rcpp::ComplexVector y(n);
  const ad* X = adptr(x);
  ad* Y = adptr(y);
#define CALL(OP) for (size_t i=0; i<n; i++) Y[i] = OP ( X[i] )
#define CUMC(OP) for (size_t i=1; i<n; i++) Y[i] = Y[i-1] OP X[i];
  if (!op.compare("abs")) CALL(fabs);
  else if (!op.compare("sqrt")) CALL(sqrt);
  else if (!op.compare("exp")) CALL(exp);
  else if (!op.compare("log")) CALL(log);
  else if (!op.compare("cos")) CALL(cos);
  else if (!op.compare("sin")) CALL(sin);
  else if (!op.compare("tan")) CALL(tan);
  else if (!op.compare("acos")) CALL(acos);
  else if (!op.compare("asin")) CALL(asin);
  else if (!op.compare("atan")) CALL(atan);
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
    ans = 0.; REDUCE(+);
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
  if (!method.compare("plain"))
    Z = X * Y;
  else if (!method.compare("atomic"))
    Z = atomic::matmul(matrix<ad>(X), matrix<ad>(Y));
  else if (!method.compare("TMBad"))
    Z = TMBad::matmul(matrix<ad>(X), matrix<ad>(Y));
  else
    Rf_error("Method '%s' not implemented", method.c_str());
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
  Rcpp::ComplexVector z(x.ncol());
  ConstMapMatrix X((ad*) x.begin(), x.nrow(), x.ncol());
  ConstMapMatrix S((ad*) s.begin(), s.nrow(), s.ncol());
  auto nldens = density::MVNORM(matrix<ad>(S));
  for (int j=0; j < X.cols(); j++) {
    ad ans = -nldens(vector<ad>(X.col(j)));
    if (!give_log) ans = exp(ans);
    z[j] = ad2cplx(ans);
  }
  return as_advector(z);
}


// ============================== Sparse matrices
Eigen::SparseMatrix<ad> SparseInput(const Rcpp::ComplexVector &x) {
  CHECK_INPUT(x);
  Rcpp::IntegerVector i = x.attr("i");
  Rcpp::IntegerVector p = x.attr("p");
  Rcpp::IntegerVector Dim = x.attr("Dim");
  return
    Eigen::Map<const Eigen::SparseMatrix<ad> > (Dim[0], // rows()
                                                Dim[1], // cols()
                                                i.size(), // nonZeros()
                                                p.begin(), // outerIndexPtr()
                                                i.begin(), // innerIndexPtr()
                                                (ad*) x.begin(), // data()
                                                NULL); // innerNonZeroPtr();
}

Rcpp::ComplexVector SparseOutput (const Eigen::SparseMatrix<ad> &x) {
  size_t nnz  = x.nonZeros();
  Rcpp::IntegerVector Dim(2);
  Dim[0] = x.rows();
  Dim[1] = x.cols();
  Rcpp::IntegerVector i(x.innerIndexPtr(), x.innerIndexPtr() + nnz);
  Rcpp::IntegerVector p(x.outerIndexPtr(), x.outerIndexPtr() + Dim[1] + 1);
  Rcpp::ComplexVector ans( (Rcomplex*) x.valuePtr(), (Rcomplex*) (x.valuePtr() + nnz));
  ans.attr("i") = i;
  ans.attr("p") = p;
  ans.attr("Dim") = Dim;
  return as_advector(ans);
}

// [[Rcpp::export]]
Rcpp::ComplexVector testSparse(const Rcpp::ComplexVector &x) {
  using TMBad::operator<<;
  Eigen::SparseMatrix<ad> X = SparseInput(x);
  Rcout << X << "\n";
  Eigen::SparseMatrix<ad> Y = X+X;
  Rcout << Y << "\n";
  return SparseOutput(Y);
}

// [[Rcpp::export]]
Rcpp::ComplexVector SparseArith2(const Rcpp::ComplexVector &x, const Rcpp::ComplexVector &y, std::string op) {
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  Eigen::SparseMatrix<ad> X = SparseInput(x);
  Eigen::SparseMatrix<ad> Y = SparseInput(y);
  Rcpp::ComplexVector z;
  if (!op.compare("+"))      z = SparseOutput(X + Y);
  else if (!op.compare("-")) z = SparseOutput(X - Y);
  else if (!op.compare("*")) z = SparseOutput(X * Y);
  else Rf_error("'%s' not implemented", op.c_str());
  return z;
}
