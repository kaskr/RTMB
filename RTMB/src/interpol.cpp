// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// [[Rcpp::export]]
Rcpp::XPtr<tmbutils::interpol2D<double> >
ip2D(Rcpp::NumericMatrix data,
     Rcpp::NumericVector x_range,
     Rcpp::NumericVector y_range,
     Rcpp::List con) {
  typedef tmbutils::interpol2D<double> ip2D_t;
  tmbutils::interpol2D_config<double> cfg;
  cfg.safe_check = false; // redundant (data is numeric)
  Rcpp::NumericVector R = con["R"];
  cfg.R = R[0];
  ip2D_t* ptr = new ip2D_t (asMatrix<double>(data),
                            asVector<double>(x_range),
                            asVector<double>(y_range),
                            cfg);
  return Rcpp::XPtr<ip2D_t> (ptr);
}

// [[Rcpp::export]]
Rcpp::NumericVector
ip2D_eval_num(Rcpp::XPtr<tmbutils::interpol2D<double> > ptr,
              Rcpp::NumericVector x,
              Rcpp::NumericVector y) {
  size_t nx = x.size();
  size_t ny = y.size();
  size_t n = std::max(nx, ny);
  Rcpp::NumericVector z(n);
  for (size_t i=0; i<n; i++) {
    z[i] = (*ptr)(x[i % nx], y[i % ny]);
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexVector
ip2D_eval_ad(Rcpp::XPtr<tmbutils::interpol2D<double> > ptr,
             Rcpp::ComplexVector x,
             Rcpp::ComplexVector y) {
  if (!ad_context())
    Rcpp::stop("'ip2D_eval_ad' requires an active tape");
  CHECK_INPUT(x);
  CHECK_INPUT(y);
  size_t nx = x.size();
  size_t ny = y.size();
  size_t n = std::max(nx, ny);
  Rcpp::ComplexVector z(n);
  for (size_t i=0; i<n; i++) {
    ad xi = cplx2ad(x[i % nx]);
    ad yi = cplx2ad(y[i % ny]);
    ad zi = (*ptr) (xi, yi);
    z[i] = ad2cplx(zi);
  }
  return as_advector(z);
}
