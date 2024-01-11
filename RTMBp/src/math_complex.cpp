// [[Rcpp::depends(TMB)]]

/*
  ====================================================================

  This file implements _complex AD_.

  WARNING: Not to confuse with the RComplex representation of AD types

  ====================================================================
*/

#include "RTMB.h"
#include <unsupported/Eigen/FFT>

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
