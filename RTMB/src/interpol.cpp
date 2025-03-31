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
ADrep
ip2D_eval_ad(Rcpp::XPtr<tmbutils::interpol2D<double> > ptr,
             ADrep x,
             ADrep y) {
  if (!ad_context())
    Rcpp::stop("'ip2D_eval_ad' requires an active tape");
  size_t nx = x.size();
  size_t ny = y.size();
  size_t n = std::max(nx, ny);
  ADrep z(n);
  ad* X = adptr(x);
  ad* Y = adptr(y);
  ad* Z = adptr(z);
  for (size_t i=0; i<n; i++) {
    Z[i] = (*ptr) (X[i % nx], Y[i % ny]);
  }
  return z;
}

/* ======================= SPLINE ==================================== */

// [[Rcpp::export]]
Rcpp::XPtr<tmbutils::splinefun<ad> >
splineptr(ADrep x,
          ADrep y,
          int method=3) {
  typedef tmbutils::splinefun<ad> spline_t;
  spline_t* ptr = new spline_t (x, y, method);
  return Rcpp::XPtr<spline_t> (ptr);
}

// [[Rcpp::export]]
ADrep
splineptr_eval(Rcpp::XPtr<tmbutils::splinefun<ad> > ptr,
               ADrep x) {
  vector<ad> y = (*ptr)(x);
  return y;
}
