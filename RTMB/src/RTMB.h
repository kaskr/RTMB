#include <Rcpp.h>
#include "TMB.h"


typedef TMBad::ad_aug ad;
typedef Eigen::Map<Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > MapMatrix;
typedef Eigen::Map<const Eigen::Matrix<ad, Eigen::Dynamic, Eigen::Dynamic> > ConstMapMatrix;

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
extern tape_config_t tape_config;

Rcomplex ad2cplx(const ad &x);
ad cplx2ad(const Rcomplex &x);
ad* adptr(const Rcpp::ComplexVector &x);
bool is_advector (SEXP x);
bool is_adsparse (SEXP x);
bool is_adscalar (SEXP x);
bool valid(const ad &x);
bool valid(Rcpp::ComplexVector x);
bool ad_context();
Rcpp::ComplexVector& as_advector(Rcpp::ComplexVector &x);

#define CHECK_INPUT(x)                                                  \
if (!is_advector(x))                                                    \
  Rcpp::stop("'" #x "' must be 'advector' (lost class attribute?)" );   \
 if (!valid(Rcpp::ComplexVector(x)))                                    \
  Rcpp::stop("'" #x "' is not a valid 'advector' (constructed using illegal operation?)" );

Eigen::SparseMatrix<ad> SparseInput(Rcpp::S4 x);
Rcpp::S4 SparseOutput (const Eigen::SparseMatrix<ad> &S);
ad ScalarInput(SEXP x_);
