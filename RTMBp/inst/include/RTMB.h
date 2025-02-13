/* ========================================================================== */
/* RTMB interface (R <-> C++) */
/* ========================================================================== */

#include <Rcpp.h>
#include "TMB.h"

// AD type shorthands
typedef TMBad::ad_aug ad;
typedef Eigen::Map<Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > MapMatrix;
typedef Eigen::Map<const Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > ConstMapMatrix;

// AD representation
struct ADrep : Rcpp::RObject {
  typedef Rcpp::RObject Base;
  // Grant valid object by setting the class attribute and S4 bit
  void setclass();
  // Default CTOR
  ADrep ();
  // Input and check validity
  ADrep (Rcpp::RObject x);
  // Output vector
  ADrep (size_t n);
  ADrep (const ad* begin, const ad* end);
  // Output matrix
  ADrep (size_t n, size_t m);
  // Pointer to first element
  ad* adptr();
  // Number of elements
  size_t size();
  // Dimension *if* a matrix
  size_t nrow();
  size_t ncol();
};

// AD scalar conversions ( Rcomplex <-> ad )
Rcomplex ad2cplx(const ad &x);
ad       cplx2ad(const Rcomplex &x);
ad       ScalarInput(SEXP x_);
// AD vector view of R representation
ad* adptr(ADrep x);
Rcpp::ComplexVector& as_advector(Rcpp::ComplexVector &x);
// Dense matrix input/output
ConstMapMatrix MatrixInput  (ADrep x);
ADrep          MatrixOutput (const matrix<ad> &X);
// Sparse matrix input/output
Eigen::SparseMatrix<ad> SparseInput  (Rcpp::RObject x);
Rcpp::RObject           SparseOutput (const Eigen::SparseMatrix<ad> &S);

// Tests
bool is_advector (SEXP x);
bool is_admatrix (SEXP x);
bool is_adsparse (SEXP x);
bool is_adscalar (SEXP x);
bool valid(const ad &x);
bool valid(ADrep x);
bool ad_context();

// AD input validity check
#define CHECK_INPUT(x)                                                  \
if (!is_advector(x))                                                    \
  Rcpp::stop("'" #x "' must be 'advector' (lost class attribute?)" );   \
 if (!valid(x))                                                         \
  Rcpp::stop("'" #x "' is not a valid 'advector' (constructed using illegal operation?)" );

// Global Tape Configuration
struct tape_config_t {
  int comparison; // Safe=0 / Taped=1 / Unsafe=2
  int atomic;     // No atomic=0 / Use atomic=1
  int vectorize;  // No vectorize =0 / Use vectorize=1
  tape_config_t();
  bool matmul_plain    ();
  bool matmul_atomic   ();
  bool matmul_TMBad    ();
  bool ops_vectorize   ();
  bool math_vectorize  ();
  bool sum_vectorize   ();
  bool compare_forbid  ();
  bool compare_taped   ();
  bool compare_allow   ();
  bool mvnorm_atomic   ();
};

// Tweak RcppExports
void ptr_forward(TMBad::ADFun<>* adf);
#define RTMB_CCALLABLES                                                 \
  TMB_CCALLABLES("RTMB")                                                \
  R_RegisterCCallable("RTMB", "ptr_forward", (DL_FUNC) &ptr_forward);
