// [[Rcpp::depends(TMB)]]
#include <Rcpp.h>
using Rcpp::Rcout;
using Rcpp::Rcerr;
#define TMBAD_ASSERT2(x,msg)                                            \
if (!(x)) {                                                             \
  Rcerr << "TMBad assertion failed.\n";                                 \
  Rcerr << "The following condition was not met: " << #x << "\n";       \
  Rcerr << "Possible reason: " msg << "\n";                             \
  Rcerr << "For more info run your program through a debugger.\n";      \
  abort();                                                              \
}
#define TMBAD_ASSERT(x) TMBAD_ASSERT2(x,"Unknown")
#include <TMBad/TMBad.hpp>
#include <TMBad/TMBad.cpp>

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
  ;
}

/* ========================================================================== */
/* AD vector object */
/* ========================================================================== */

typedef TMBad::ad_aug ad;
typedef std::vector<ad> ad_vec;

// [[Rcpp::export]]
Rcpp::XPtr< ad_vec > advec(std::vector<double> x) {
  ad_vec* pv = new ad_vec(x.begin(), x.end());
  Rcpp::XPtr< ad_vec > ans(pv);
  return ans;
}
// [[Rcpp::export]]
void Dependent(Rcpp::XPtr< ad_vec > x) {
  Dependent(*x);
}
// [[Rcpp::export]]
void Independent(Rcpp::XPtr< ad_vec > x) {
  Independent(*x);
}

/* Arith:

   • ‘"+"’, ‘"-"’, ‘"*"’, ‘"/"’, ‘"^"’, ‘"%%"’, ‘"%/%"’

   • ‘"&"’, ‘"|"’, ‘"!"’

   • ‘"=="’, ‘"!="’, ‘"<"’, ‘"<="’, ‘">="’, ‘">"’
*/

// [[Rcpp::export]]
Rcpp::XPtr<ad_vec> Arith2(Rcpp::XPtr<ad_vec> x, Rcpp::XPtr<ad_vec> y, std::string op) {
  size_t nx = x->size(), ny = y->size();
  size_t n = std::max(nx, ny);
  ad_vec z(n);
#define CALL(OP) for (size_t i=0; i<n; i++) z[i] = (*x)[i % nx] OP (*y)[i % ny]
  if (!op.compare("+")) CALL(+);
  else if (!op.compare("-")) CALL(-);
  else if (!op.compare("*")) CALL(*);
  else if (!op.compare("/")) CALL(/);
  else if (!op.compare("^")) {
    for (size_t i=0; i<n; i++) z[i] = pow( (*x)[i % nx] , (*y)[i % ny]);
  }
  else Rf_error("Not implemented");
#undef CALL
  ad_vec* pv = new ad_vec (z);
  Rcpp::XPtr< ad_vec > ans(pv);
  return ans;
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
Rcpp::XPtr<ad_vec> Math1(Rcpp::XPtr<ad_vec> x, std::string op) {
  size_t n = x->size();
  ad_vec y(n);
#define CALL(OP) for (size_t i=0; i<n; i++) y[i] = OP( (*x)[i] )
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
  else if (!op.compare("cumsum")) {
    for (size_t i=0; i<n; i++) y[i] = ( i>0 ? y[i-1] : 0) + (*x)[i];
  }
  else Rf_error("Not implemented");
#undef CALL
  ad_vec* pv = new ad_vec (y);
  Rcpp::XPtr< ad_vec > ans(pv);
  return ans;
}

// [[Rcpp::export]]
Rcpp::XPtr<ad_vec> Reduce1(Rcpp::XPtr<ad_vec> x, std::string op) {
  size_t n = x->size();
  ad_vec y(1);
#define REDUCE(OP) for (size_t i=0; i<n; i++) y[0] = y[0] OP (*x)[i];
  if (!op.compare("sum")) {
    y[0] = 0.; REDUCE(+);
  } else if (!op.compare("prod")) {
    y[0] = 1.; REDUCE(*);
  } else
    Rf_error("Not implemented");
#undef REDUCE
  ad_vec* pv = new ad_vec (y);
  Rcpp::XPtr< ad_vec > ans(pv);
  return ans;
}

// [[Rcpp::export]]
Rcpp::XPtr<ad_vec> Subset(Rcpp::XPtr<ad_vec> x, Rcpp::IntegerVector i) {
  ad_vec* pv = new ad_vec (i.size());
  for (int j=0; j<i.size(); j++) (*pv)[j] = (*x)[i[j]];
  Rcpp::XPtr< ad_vec > ans(pv);
  return ans;  
}

// [[Rcpp::export]]
void AppendInplace(Rcpp::XPtr<ad_vec> x, Rcpp::XPtr<ad_vec> y) {
  (*x).insert((*x).end(), (*y).begin(), (*y).end());
}

// [[Rcpp::export]]
Rcpp::IntegerVector Length(Rcpp::XPtr<ad_vec> x) {
  Rcpp::IntegerVector ans(1);
  ans[0] = (*x).size();
  return ans;
}

// [[Rcpp::export]]
void Display(Rcpp::XPtr<ad_vec> x) {
  std::cout << "{";
  for (size_t i=0; i<(*x).size(); i++) {
    if (i>0) Rcout << ", ";
    Rcout << (*x)[i].Value() ;
  }
  std::cout << "}\n";
}
