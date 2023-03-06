#include "RTMB.h"

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
void atomic_transform(TMBad::ADFun<>* adf) {
  *adf = (*adf).atomic();
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
// [[Rcpp::export]]
Rcpp::S4 get_graph(Rcpp::XPtr<TMBad::ADFun<> > adf) {
  // reverse row-major == forward col-major
  TMBad::graph G = (*adf).glob.reverse_graph();
  size_t n = (*adf).glob.opstack.size();
  Rcpp::StringVector names(n);
  for (size_t i=0; i<n; i++) {
    names[i] = (*adf).glob.opstack[i]->op_name();
    std::sort(G.j.begin() + G.p[i], G.j.begin() + G.p[i+1]);
  }
  Rcpp::S4 ans("ngCMatrix");
  ans.slot("i") = Rcpp::IntegerVector(G.j.begin(), G.j.end());
  ans.slot("p") = Rcpp::IntegerVector(G.p.begin(), G.p.end());
  ans.slot("Dim") = Rcpp::IntegerVector::create(n, n);
  ans.slot("Dimnames") = Rcpp::List::create(names, names);
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
  .method("atomic", &atomic_transform)
  .method("ptrTMB", &ptrTMB)
  ;
}

tape_config_t::tape_config_t() : comparison(0), atomic(1), vectorize(0) {}
bool tape_config_t::matmul_plain() { return (atomic == 0); }
bool tape_config_t::matmul_atomic   () { return (atomic == 1) && vectorize==0; }
bool tape_config_t::matmul_TMBad    () { return (atomic == 1) && vectorize==1; }
bool tape_config_t::ops_vectorize   () { return vectorize == 1; }
bool tape_config_t::math_vectorize  () { return vectorize == 1; }
bool tape_config_t::sum_vectorize   () { return vectorize == 1; }
bool tape_config_t::compare_forbid  () { return comparison == 0; }
bool tape_config_t::compare_taped   () { return comparison == 1; }
bool tape_config_t::compare_allow   () { return comparison == 2; }
bool tape_config_t::mvnorm_atomic   () { return (atomic == 1); }
tape_config_t tape_config;

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
                            GET(compare_forbid),
                            GET(compare_taped),
                            GET(compare_allow),
                            GET(mvnorm_atomic));
}
// [[Rcpp::export]]
bool compare_allow() { return tape_config.compare_allow(); }


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
bool is_adsparse (SEXP x) {
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
