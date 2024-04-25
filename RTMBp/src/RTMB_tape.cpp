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
void ad_print(TMBad::ADFun<>* adf, int depth=0) {
  TMBad::print_config cfg;
  cfg.depth = depth;
  adf->print(cfg);
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
Rcpp::S4 SpJacFun(Rcpp::XPtr<TMBad::ADFun<> > adf) {
  TMBad::Sparse<TMBad::ADFun<> > sadf = adf->SpJacFun();
  Rcpp::S4 ans("ngTMatrix");
  ans.slot("i") = Rcpp::IntegerVector(sadf.i.begin(), sadf.i.end());
  ans.slot("j") = Rcpp::IntegerVector(sadf.j.begin(), sadf.j.end());
  ans.slot("Dim") = Rcpp::IntegerVector::create(sadf.m, sadf.n);
  TMBad::ADFun<>* pf = new TMBad::ADFun<>(sadf);
  ans.attr("tape") = Rcpp::XPtr<TMBad::ADFun<> > (pf);
  return ans;
}
// [[Rcpp::export]]
void RangeProj(Rcpp::XPtr<TMBad::ADFun<> > adf, Rcpp::IntegerVector i) {
  std::vector<TMBad::Index>& di = (*adf).glob.dep_index;
  Rcpp::IntegerVector di_ = Rcpp::IntegerVector(di.begin(), di.end());
  di_ = di_[i]; // Rcpp handles out-of-bound error
  di = std::vector<TMBad::Index>(di_.begin(), di_.end());
}
// [[Rcpp::export]]
Rcpp::IntegerVector find_op_by_name(Rcpp::XPtr<TMBad::ADFun<> > adf, Rcpp::String name) {
  std::vector<TMBad::Index> nodes = find_op_by_name(adf->glob, name.get_cstring());
  return Rcpp::IntegerVector(nodes.begin(), nodes.end());
}
// [[Rcpp::export]]
Rcpp::IntegerVector op2var(Rcpp::XPtr<TMBad::ADFun<> > adf, Rcpp::IntegerVector nodes) {
  std::vector<TMBad::Index> var =
    adf->glob.op2var( std::vector<TMBad::Index>(nodes.begin(), nodes.end()) ) ;
  return Rcpp::IntegerVector(var.begin(), var.end());
}
// [[Rcpp::export]]
Rcpp::IntegerVector findIndex(Rcpp::XPtr<TMBad::ADFun<> > adf, Rcpp::String name) {
  std::vector<TMBad::Index> nodes = find_op_by_name(adf->glob, name.get_cstring());
  std::vector<TMBad::Index> var = adf->glob.op2var( nodes ) ;
  return Rcpp::IntegerVector(var.begin(), var.end());
}
// [[Rcpp::export]]
void setinvIndex(Rcpp::XPtr<TMBad::ADFun<> > adf, Rcpp::IntegerVector index) {
  adf->glob.inv_index = std::vector<TMBad::Index>(index.begin(), index.end());
}
// [[Rcpp::export]]
Rcpp::IntegerVector getinvIndex(Rcpp::XPtr<TMBad::ADFun<> > adf) {
  return Rcpp::IntegerVector(adf->glob.inv_index.begin(), adf->glob.inv_index.end());
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
#define SET(name) if (name != -1) tape_config.name = name;
  SET(comparison);
  SET(atomic);
  SET(vectorize);
#undef SET
  // Current settings
#define GET(name) Rcpp::Named(#name) = tape_config.name
  return Rcpp::List::create(
                            GET(comparison),
                            GET(atomic),
                            GET(vectorize));
#undef GET
}
// [[Rcpp::export]]
Rcpp::List get_tape_config() {
  // Current derived flags
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
#undef GET
}
// [[Rcpp::export]]
bool compare_allow() { return tape_config.compare_allow(); }

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

// CCallable
void ptr_forward(TMBad::ADFun<>* adf) {
  adf->forward();
}

// [[Rcpp::export]]
Rcpp::XPtr<double> ptr_getx (Rcpp::XPtr<TMBad::ADFun<> > adf) {
  std::vector<TMBad::Index> idx = (*adf).glob.inv_index;
  if (idx.size() == 0) Rcpp::stop("Tape has no inputs");
  for (size_t i=1; i<idx.size(); i++) {
    if (idx[i] - idx[i-1] != 1)
      Rcpp::stop("Tape has Non-consecutive inputs");
  }
  double* ptr = (*adf).glob.values.data() + (*adf).glob.inv_index[0];
  Rcpp::XPtr<double> ans = Rcpp::XPtr<double> (ptr, false); // No finalizer!
  ans.attr("size") = Rcpp::IntegerVector::create(idx.size());
  return ans;
}

// [[Rcpp::export]]
Rcpp::XPtr<double> ptr_gety (Rcpp::XPtr<TMBad::ADFun<> > adf) {
  std::vector<TMBad::Index> idx = (*adf).glob.dep_index;
  if (idx.size() == 0) Rcpp::stop("Tape has no outputs");
  for (size_t i=1; i<idx.size(); i++) {
    if (idx[i] - idx[i-1] != 1)
      Rcpp::stop("Tape has Non-consecutive outputs");
  }
  double* ptr = (*adf).glob.values.data() + (*adf).glob.dep_index[0];
  Rcpp::XPtr<double> ans = Rcpp::XPtr<double> (ptr, false); // No finalizer!
  ans.attr("size") = Rcpp::IntegerVector::create(idx.size());
  return ans;
}
