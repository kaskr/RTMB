// [[Rcpp::depends(TMB)]]

/*
  ====================================================================

  This file implements _complex AD_.

  WARNING: Not to confuse with the RComplex representation of AD types

  ====================================================================
*/

#include "RTMB.h"
#include <unsupported/Eigen/FFT>

// Some extra interface helper functions
typedef std::complex<ad> adcplx;

// [[Rcpp::export]]
Rcpp::ComplexVector Arith2_complex(const Rcpp::ComplexVector &x,
                                   const Rcpp::ComplexVector &y,
                                   std::string op) {
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  size_t nx = x.size(), ny = y.size();
  size_t n = (std::min(nx, ny) > 0 ? std::max(nx, ny) : 0);
  Rcpp::ComplexVector z(n);
  // AD complex
  adcplx* X = (adcplx*) adptr(x);
  adcplx* Y = (adcplx*) adptr(y);
  adcplx* Z = (adcplx*) adptr(z);
  size_t ncplx = n / 2;
#define CALL(OP)                                                \
  for (size_t i=0; i<ncplx; i++) Z[i] = X[i % nx] OP Y[i % ny];
  if (!op.compare("+")) CALL(+)
  else if (!op.compare("-")) CALL(-)
  else if (!op.compare("*")) CALL(*)
  else if (!op.compare("/")) CALL(/)
  else if (!op.compare("^")) {
    for (size_t i=0; i<n; i++) Z[i] = pow(X[i % nx] , Y[i % ny]);
  }
  else Rcpp::stop("'%s' not implemented", op.c_str());
#undef CALL
  return as_advector(z);
}

// [[Rcpp::export]]
Rcpp::ComplexVector Math1_complex(const Rcpp::ComplexVector &x, std::string op) {
  CHECK_INPUT(x);
  size_t n = x.size();
  Rcpp::ComplexVector y(n);
  // AD complex
  adcplx* X = (adcplx*) adptr(x);
  adcplx* Y = (adcplx*) adptr(y);
  size_t ncplx = n / 2;
#define CALL(OP) for (size_t i=0; i<ncplx; i++) Y[i] = OP ( X[i] )
#define CUMC(OP) for (size_t i=1; i<ncplx; i++) Y[i] = Y[i-1] OP X[i];
  if (!op.compare("abs")) CALL(std::abs);
  else if (!op.compare("sqrt")) CALL(sqrt);
  else if (!op.compare("exp")) CALL(exp);
  else if (!op.compare("log")) CALL(log);
  else if (!op.compare("cos")) CALL(cos);
  else if (!op.compare("sin")) CALL(sin);
  else if (!op.compare("tan")) CALL(tan);
  else if (!op.compare("acos")) CALL(acos);
  else if (!op.compare("asin")) CALL(asin);
  else if (!op.compare("atan")) CALL(atan);
  else if (!op.compare("cosh")) CALL(cosh);
  else if (!op.compare("sinh")) CALL(sinh);
  else if (!op.compare("tanh")) CALL(tanh);
  // FIXME:
  // else if (!op.compare("acosh")) CALL(acosh);
  // else if (!op.compare("asinh")) CALL(asinh);
  // else if (!op.compare("atanh")) CALL(atanh);
  // else if (!op.compare("lgamma")) CALL(lgamma);
  // else if (!op.compare("gamma")) CALL(rtmb_gamma);
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

namespace TMBad {

template<bool adjoint=false>
void fft_array(std::complex<double>* x,
               std::vector<size_t> dim) {
  Eigen::FFT<double> fft;
  fft.SetFlag(fft.Unscaled);
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Map<Matrix> Mat;
  vector<std::complex<double> > buf;
  size_t n = TMBad::prod_int(dim);
  for (size_t i=0; i<dim.size(); i++) {
    int nrow = dim[i];
    int ncol = n / dim[i];
    Mat X(x, nrow, ncol);
    buf.resize(nrow);
    for (int j=0; j<ncol; j++) {
      std::complex<double>* src = X.col(j).data();
      std::complex<double>* dest = buf.data();
      if (!adjoint)
        fft.fwd(dest, src, nrow);
      else
        fft.inv(dest, src, nrow);
      X.col(j).array() = buf;
    }
    if ((nrow != 1) && (ncol != 1)) {
      Matrix XT = X.transpose();
      XT.resize(nrow, ncol);
      X = XT;
    }
  }
}

template<bool adjoint=false>
struct FFTOp : global::DynamicOperator< -1 , -1 > {
  typedef std::complex<double> cplx;
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  size_t n;
  std::vector<size_t> dim;
  Index input_size()  const { return n; }
  Index output_size() const { return n; }
  FFTOp (size_t n, std::vector<size_t> dim) : n(n), dim(dim) { }
  void forward(ForwardArgs<double> &args) {
    args.y_segment(0, n) = args.x_segment(0, n);
    fft_array<adjoint>( (cplx*) args.y_ptr(0), dim);
  }
  void reverse(ReverseArgs<double> &args) {
    std::vector<double> buf = args.dy_segment(0, n);
    fft_array<!adjoint>( (cplx*) buf.data(), dim);
    args.dx_segment(0, n) += buf;
  }
  void reverse(ReverseArgs<Replay> &args) {
    std::vector<Replay> buf = args.dy_segment(0, n);
    args.dx_segment(0, n) += global::Complete<FFTOp<!adjoint> >(n, dim)(buf);
  }
  template <class Type> void forward(ForwardArgs<Type> &args) {
    TMBAD_ASSERT2(false, "FFT forward not implemented for this type");
  }
  template <class Type> void reverse(ReverseArgs<Type> &args) {
    TMBAD_ASSERT2(false, "FFT reverse not implemented for this type");
  }
  const char* op_name() {return adjoint ? "iFFT" : "FFT";}
};
}

// [[Rcpp::export]]
Rcpp::ComplexVector fft_complex(const Rcpp::ComplexVector &x, std::vector<size_t> dim, bool inverse=false) {
  CHECK_INPUT(x);
  size_t n = x.size();
  // Check dim
  if (TMBad::prod_int(dim) * 2 != n)
    Rcpp::stop("prod(dim) must equal length(x)/2");
  ad* X = adptr(x);
  // Add to tape
  std::vector<ad> X_(X, X + n);
  std::vector<ad> Y_;
  if (inverse)
    Y_ = TMBad::global::Complete<TMBad::FFTOp<true> >(n, dim) (X_);
  else
    Y_ = TMBad::global::Complete<TMBad::FFTOp<false> >(n, dim) (X_);
  // Pass to R
  Rcpp::ComplexVector y(n);
  for (size_t j=0; j < n; j++) {
    y[j] = ad2cplx(Y_[j]);
  }
  return as_advector(y);
}
