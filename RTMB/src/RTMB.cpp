#include "RTMB.h"

inline Rcpp::ComplexVector unwrap(ADrep x) {
  return Rcpp::ComplexVector(x);
}
inline Rcpp::ComplexMatrix unwrap_matrix(ADrep x) {
  return Rcpp::ComplexMatrix(x);
}
void ADrep::setclass () {
  Base::operator=( Rf_asS4(*this, TRUE, FALSE) );
  (*this).attr("class") = "advector";
}
// Default CTOR
ADrep::ADrep () {}
// Input
ADrep::ADrep (Rcpp::RObject x) : Rcpp::RObject(x) {
  CHECK_INPUT(*this);
}
ADrep::ADrep (SEXP x) : ADrep(Rcpp::RObject(x)) { }
// Output
ADrep::ADrep (size_t n) {
  Base::operator=( Rcpp::ComplexVector(n) );
  this -> setclass();
}
ADrep::ADrep (const ad* begin, const ad* end) {
  Base::operator=( Rcpp::ComplexVector((Rcomplex*) begin, (Rcomplex*) end) );
  this -> setclass();
}
ADrep::ADrep (size_t n, size_t m) {
  Base::operator=( Rcpp::ComplexMatrix(n, m) );
  this -> setclass();
}
ad* ADrep::adptr() {
  static_assert(sizeof(ad) == sizeof(Rcomplex),
                "ad size must match Rcomplex");
  Rcpp::ComplexVector x = unwrap(*this);
  ad* px = (x.size() > 0 ? (ad*)(x.begin()) : NULL);
  return px;
}
size_t ADrep::size() { return unwrap(*this).size(); }
size_t ADrep::nrow() { return unwrap_matrix(*this).nrow(); }
size_t ADrep::ncol() { return unwrap_matrix(*this).ncol(); }
ADrep::operator vector<ad>() {
  return Eigen::Map<Eigen::Array<ad, -1, 1> > ((*this).adptr(), (*this).size());
}
ADrep::ADrep (const vector<ad> &x) : ADrep(x.data(), x.data() + x.size()) { }


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
ad* adptr(ADrep x) {
  return x.adptr();
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
  return is_advector(x) &&
    (Rcpp::ComplexVector(x).size() == 1) &&
    !Rcpp::ComplexVector(x).hasAttribute("dim");
}
bool valid(const ad &x) {
  return
    !x.ontape() || x.in_context_stack(x.data.glob);
}
// [[Rcpp::export]]
bool valid(ADrep x) {
  ad* X = x.adptr();
  size_t n = x.size();
  for (size_t i=0; i<n; i++)
    if (!valid(X[i])) return false;
  return true;
}
// [[Rcpp::export]]
bool ad_context() {
  return TMBad::get_glob() != NULL;
}

// [[Rcpp::export]]
Rcpp::ComplexVector& as_advector(Rcpp::ComplexVector &x) {
  x = Rf_asS4(x, TRUE, FALSE); // Was: SET_S4_OBJECT(x);
  x.attr("class") = "advector";
  return x;
}

// ============================== Dense matrices
ConstMapMatrix MatrixInput(ADrep x) {
  return ConstMapMatrix (x.adptr(), x.nrow(), x.ncol());
}
ADrep MatrixOutput(const matrix<ad> &X) {
  ADrep z(X.rows(), X.cols());
  MapMatrix Z(z.adptr(), z.nrow(), z.ncol());
  Z = X;
  return z;
}

// ============================== Sparse matrices
Eigen::SparseMatrix<ad> SparseInput(Rcpp::RObject S_) {
  Rcpp::S4 S(S_);
  ADrep x(Rcpp::RObject(S.slot("x")));
  Rcpp::IntegerVector i = S.slot("i");
  Rcpp::IntegerVector p = S.slot("p");
  Rcpp::IntegerVector Dim = S.slot("Dim");
  return
    Eigen::Map<const Eigen::SparseMatrix<ad> > (Dim[0], // rows()
                                                Dim[1], // cols()
                                                i.size(), // nonZeros()
                                                p.begin(), // outerIndexPtr()
                                                i.begin(), // innerIndexPtr()
                                                x.adptr(), // data()
                                                NULL); // innerNonZeroPtr();
}

Rcpp::RObject SparseOutput (const Eigen::SparseMatrix<ad> &S) {
  size_t nnz  = S.nonZeros();
  Rcpp::IntegerVector Dim(2);
  Dim[0] = S.rows();
  Dim[1] = S.cols();
  Rcpp::IntegerVector i(S.innerIndexPtr(), S.innerIndexPtr() + nnz);
  Rcpp::IntegerVector p(S.outerIndexPtr(), S.outerIndexPtr() + Dim[1] + 1);
  ADrep x( S.valuePtr(), S.valuePtr() + nnz);
  Rcpp::S4 ans("adsparse");
  ans.slot("x") = x;
  ans.slot("i") = i;
  ans.attr("p") = p;
  ans.attr("Dim") = Dim;
  return Rcpp::RObject(ans);
}

ad ScalarInput(SEXP x_) {
  ADrep x(x_);
  return adptr(x)[0];
}

TMBad::ad_segment ad_segment(ADrep x) {
  if (!ad_context())
    Rcpp::stop("'ad_segment' requires an active ad context");
  TMBad::global* cur_glob = TMBad::get_glob();
  size_t before = cur_glob->opstack.size();
  TMBad::ad_segment ans(x.adptr(), x.size());
  size_t after = cur_glob->opstack.size();
  if (before != after) {
    // Check safe mutation rules
    // - Variables can safely mutate to a different index within same context.
    // - Constants can safely mutate to the root context.
    bool root = cur_glob->parent_glob == NULL; // Active context = root ?
    bool safe = root || ans.all_on_active_tape(x.adptr(), x.size());
    if (safe) { // mutate
      size_t n = x.size();
      ad* xp = x.adptr();
      for (size_t i=0; i<n; i++) xp[i] = ans[i];
    }
  }
  return ans;
}
