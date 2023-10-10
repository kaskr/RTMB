#include "RTMB.h"

// Dummy
template<class Type>
Type objective_function<Type>::operator() () {
  // pretend all parameters have been read
  this->index = this->theta.size();
  return 0;
}
// Force instantiation (needed on some platforms)
template ad objective_function<ad>::operator() ();
template double objective_function<double>::operator() ();

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
int GetDomain(TMBad::ADFun<>* adf) {
  return (*adf).Domain();
}
int GetRange(TMBad::ADFun<>* adf) {
  return (*adf).Range();
}
Rcpp::NumericVector GetDomainVec(TMBad::ADFun<>* adf) {
  std::vector<double> ans = (*adf).DomainVec();
  return Rcpp::NumericVector(ans.begin(), ans.end());
}
Rcpp::NumericVector GetRangeVec(TMBad::ADFun<>* adf) {
  std::vector<double> ans = (*adf).RangeVec();
  return Rcpp::NumericVector(ans.begin(), ans.end());
}
// Some ADFun object transformations
void JacFun(TMBad::ADFun<>* adf) {
  // Get dimensions
  TMBad::Index n = adf->Domain();
  TMBad::Index m = adf->Range();
  // *Transposed* jacobian
  (*adf) = (*adf).JacFun();
  // Safety check - just in case
  if ((*adf).glob.dep_index.size() != m * n)
    Rcpp::stop("Invalid jacobian tape");
  // Transpose result
  Eigen::Map< Eigen::Array<TMBad::Index, -1, -1> >
    Dep ((*adf).glob.dep_index.data(), n, m);
  Eigen::Array<TMBad::Index, -1, -1> DepT = Dep.transpose();
  DepT.resize(n, m);
  Dep = DepT;
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
void eliminate(TMBad::ADFun<>* adf) {
  (*adf).eliminate();
}
void atomic_transform(TMBad::ADFun<>* adf) {
  *adf = (*adf).atomic();
}
// Defined in misc.cpp
void laplace_transform(TMBad::ADFun<>* adf, std::vector<TMBad::Index> random, SEXP config);
void newton_transform(TMBad::ADFun<>* adf, std::vector<TMBad::Index> random, SEXP config);
SEXP ptrTMB(TMBad::ADFun<>* pf) {
  SEXP ans;
#ifdef _OPENMP
  TMBad::ADFun<>* pf_cpy = new TMBad::ADFun<>();
  std::swap (*pf, *pf_cpy); // Take ownership of pf
  vector<TMBad::ADFun<>* > ppf(1);
  ppf[0] = pf_cpy;
  parallelADFun<double>* paf = new parallelADFun<double> (ppf); // 'paf' now owns memory
  SEXP ptr = Rcpp::XPtr< parallelADFun<double> >(paf, true, Rf_install("parallelADFun"));
  ans = Rcpp::List::create(Rcpp::Named("ptr") = ptr);
#else
  SEXP ptr = Rcpp::XPtr< TMBad::ADFun<> >(pf, false, Rf_install("ADFun"));
  ans = Rcpp::List::create(Rcpp::Named("ptr") = ptr);
#endif
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
// [[Rcpp::export]]
Rcpp::DataFrame get_df(Rcpp::XPtr<TMBad::ADFun<> > adf) {
  Rcpp::NumericVector values((*adf).glob.values.begin(),
                             (*adf).glob.values.end());
  Rcpp::NumericVector derivs((*adf).glob.derivs.begin(),
                             (*adf).glob.derivs.end());
  if (derivs.size() == 0) {
    derivs = Rcpp::NumericVector(values.size(), NA_REAL);
  }
  std::vector<TMBad::Index> v2o = (*adf).glob.var2op();
  Rcpp::IntegerVector node(v2o.begin(), v2o.end());
  size_t n = (*adf).glob.opstack.size();
  Rcpp::StringVector names(n);
  for (size_t i=0; i<n; i++) {
    names[i] = (*adf).glob.opstack[i]->op_name();
  }
  return
    Rcpp::DataFrame::create( Rcpp::Named("OpName") = names[node],
                             Rcpp::Named("Node") = node,
                             Rcpp::Named("Value") = values,
                             Rcpp::Named("Deriv") = derivs );
}
// [[Rcpp::export]]
void get_node(Rcpp::XPtr<TMBad::ADFun<> > adf, int node) {
  if ( (node < 0) || ( (*adf).glob.opstack.size() <= (size_t) node) )
       Rcpp::stop("'node' out of bounds");
  (*adf).glob.subgraph_cache_ptr();
  TMBad::Index n1 = (*adf).glob.opstack[node]->input_size();
  TMBad::Index n2 = (*adf).glob.opstack[node]->output_size();
  // Get node inputs
  TMBad::Args<> args((*adf).glob.inputs);
  args.ptr = (*adf).glob.subgraph_ptr[node];
  TMBad::Dependencies node_inputs;
  (*adf).glob.opstack[node]->dependencies(args, node_inputs);
  // Check dependencies
  if (node_inputs.I.size() != 0)
    Rcpp::stop("'get_node' currently cannot handle interval inputs");
  if (node_inputs.size() != n1)
    Rcpp::stop("Node input size mismatch");
  // Which of these node inputs are constant ?
  (*adf).glob.dep_index = node_inputs; // Pretend node inputs are tape outputs
  std::vector<bool> active_inputs = (*adf).activeRange();
  // New opstack
  TMBad::global::operation_stack opstack;
  opstack.push_back((*adf).glob.getOperator<TMBad::global::NullOp2>(0, n1));
  opstack.push_back((*adf).glob.opstack[node]->copy());
  // New inv index
  std::vector<TMBad::Index>
    inv_index = TMBad::which<TMBad::Index>(active_inputs);
  // New dep index
  std::vector<TMBad::Index> dep_index(n2);
  for (size_t i=0; i<n2; i++) dep_index[i] = n1+i;
  // New inputs
  std::vector<TMBad::Index> inputs(n1);
  for (size_t i=0; i<n1; i++) inputs[i] = i;
  // New values
  std::vector<TMBad::Scalar> values(n1 + n2);
  for (size_t i=0; i<n1; i++)
    values[i] = (*adf).glob.values[node_inputs[i]];
  // swap
  std::swap((*adf).glob.opstack, opstack);
  std::swap((*adf).glob.inv_index, inv_index);
  std::swap((*adf).glob.dep_index, dep_index);
  std::swap((*adf).glob.inputs, inputs);
  std::swap((*adf).glob.values, values);
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
  .method("domain", &GetDomain)
  .method("range", &GetRange)
  .method("domainvec", &GetDomainVec)
  .method("rangevec", &GetRangeVec)
  .method("jacfun", &JacFun)
  .method("parallelize", &parallelize)
  .method("fuse", &fuse)
  .method("optimize", &optimize)
  .method("eliminate", &eliminate)
  .method("atomic", &atomic_transform)
  .method("laplace", &laplace_transform)
  .method("newton", &newton_transform)
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
bool is_admatrix (SEXP x) {
  return is_advector(x) && Rcpp::ComplexVector(x).hasAttribute("dim");
}
bool is_adsparse (SEXP x) {
  return Rf_inherits(x, "adsparse");
}
bool is_adscalar (SEXP x) {
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
    if (!xad.constant())
      Rcpp::stop("Dependent 'advector' cannot be set as independent");
    xad.Independent();
    ans[i] = ad2cplx(xad);
  }
  return as_advector(ans);
}

// ============================== Dense matrices
ConstMapMatrix MatrixInput(const Rcpp::ComplexMatrix &x) {
  return ConstMapMatrix ((ad*) x.begin(), x.nrow(), x.ncol());
}
Rcpp::ComplexMatrix MatrixOutput(const matrix<ad> &X) {
  Rcpp::ComplexMatrix z(X.rows(), X.cols());
  MapMatrix Z((ad*) z.begin(), z.nrow(), z.ncol());
  Z = X;
  // FIXME: z = as_advector(z);
  z.attr("class") = "advector";
  SET_S4_OBJECT(z);
  return z;
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
