#include "RTMB.h"

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
  return is_advector(x) &&
    (Rcpp::ComplexVector(x).size() == 1) &&
    !Rcpp::ComplexVector(x).hasAttribute("dim");
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

// [[Rcpp::export]]
Rcpp::ComplexVector& as_advector(Rcpp::ComplexVector &x) {
  x = Rf_asS4(x, TRUE, FALSE); // Was: SET_S4_OBJECT(x);
  x.attr("class") = "advector";
  return x;
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
  z = Rf_asS4(z, TRUE, FALSE); // Was: SET_S4_OBJECT(z);
  z.attr("class") = "advector";
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

