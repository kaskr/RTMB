// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

/* -------------------------------------------------------------------------- */
/* --- Log space summation with hacks ! ------------------------------------- */
/* -------------------------------------------------------------------------- */

// function: log(sum(exp(x)))
// grad:     p = exp(x) / sum(exp(x))
// hessian:  diag(p) - p p^T

namespace TMBad {

std::vector<ad_plain> LSE (const std::vector<ad_plain> &x, int order); // Forward declare
template<class T>
std::vector<T> LSE (const std::vector<T> &x_, int order); // Forward declare

struct LSEOp : global::DynamicOperator< -1 , 1 > {
  size_t n;
  int order; // 0 or 1
  static const bool have_input_size_output_size = true; // FIXME: Should give compile time error if 'false'
  Index input_size()  const { return this->n; }
  Index output_size() const {
    return (order==0 ? 1 : n);
  }
  LSEOp (size_t n, int order) : n(n), order(order) { }
  void forward(ForwardArgs<Scalar> &args) {
    Scalar Max = -INFINITY;
    for (size_t i=0; i<n; i++) {
      if (Max < args.x(i)) Max = args.x(i) ;
    }
    // log(sum(exp(x)))
    if (order == 0) {
      args.y(0) = 0;
      for (size_t i=0; i<n; i++) {
        args.y(0) += exp( args.x(i) - Max );
      }
      args.y(0) = Max + log(args.y(0));
    }
    // p = exp(x) / sum(exp(x))
    if (order == 1) {
      Scalar s = 0;
      for (size_t i=0; i<n; i++) {
        args.y(i) = exp( args.x(i) - Max );
        s += args.y(i);
      }
      for (size_t i=0; i<n; i++) args.y(i) /= s;
    }   
  }
  void forward(ForwardArgs<Replay> &args) {
    // std::vector<ad_plain> x(input_size());
    // for (Index i=0; i<input_size(); i++)
    //   x[i] = args.x(i);
    // LSE(x, order);
    // args.y(0) = LSE(x, order);
    std::vector<Replay> x = args.x_segment(0, input_size());
    args.y_segment(0, output_size()) = LSE(x, order);
  }
  //template <class Type> void reverse(ReverseArgs<Type> &args) {
  void reverse(ReverseArgs<Scalar> &args) {
    if (order == 0) {
      for (size_t i=0; i<n; i++) {
        args.dx(i) += exp( args.x(i) - args.y(0) ) * args.dy(0);
      }
    }
    if (order == 1) {
      Eigen::Array<Scalar, -1, 1> p = args.y_segment(0, n);
      Eigen::Array<Scalar, -1, 1> w = args.dy_segment(0, n);
      // (diag(p) - p pT) w = (w-sum(p*w))*p
      Eigen::Array<Scalar, -1, 1> dx = ( w - (p*w).sum() ) * p;
      args.dx_segment(0, n) += dx;
    }
  }

  void reverse(ReverseArgs<Replay> &args) {
    if (order == 0) {
      std::vector<Replay> x = args.x_segment(0, n);
      std::vector<Replay> p = LSE(x, 1);
      for (size_t i=0; i<n; i++) {
        args.dx(i) += p[i] * args.dy(0);
      }
    }
    if (order == 1) {
      if (newton::get_alternative_hessian()) {
        // Report fake derivatives
        Eigen::Array<Replay, -1, 1> p = args.y_segment(0, n);
        Eigen::Array<Replay, -1, 1> w = args.dy_segment(0, n);
        // Dominant direction
        Eigen::Array<Replay, -1, 1> pu = p / sqrt((p*p).sum());
        // True Hessian = diag(p) - p pT
        // Approx Hessian = ( pu^T ( True Hessian ) pu * pu pu^T
        // 1. pu^T ( True Hessian ) pu
        //    = pu^T ( diag(p) - p p^T ) pu
        Replay tmp1 = (pu*p*pu).sum() - (pu*p).sum() * (pu*p).sum();
        Eigen::Array<Replay, -1, 1> dx = tmp1 * (pu*w).sum() * pu;
        args.dx_segment(0, n) += dx;
      } else {
        Eigen::Array<Replay, -1, 1> p = args.y_segment(0, n);
        Eigen::Array<Replay, -1, 1> w = args.dy_segment(0, n);
        Eigen::Array<Replay, -1, 1> dx = p * w - p * (p * w).sum();
        args.dx_segment(0, n) += dx;
      }
    }
  }

  
  template <class Type> void reverse(ReverseArgs<Type> &args) {
    Rcpp::stop("Not implemented");
  }
  
  const char* op_name() {return "LSSumOp";}
};
std::vector<ad_plain> LSE (const std::vector<ad_plain> &x, int order) {
  OperatorPure* pOp = get_glob()->getOperator<LSEOp>(x.size(), order);
  return get_glob()->add_to_stack<LSEOp>(pOp, x);
}
template<class T>
std::vector<T> LSE (const std::vector<T> &x_, int order) {
  std::vector<ad_plain> x(x_.begin(), x_.end());
  std::vector<ad_plain> y = LSE(x, order);
  std::vector<T> ans(y.begin(), y.end());
  return ans;
}
}

// [[Rcpp::export]]
Rcpp::ComplexVector LSE(const Rcpp::ComplexVector x) {
  CHECK_INPUT(x);
  size_t n = x.size();
  ad* X = adptr(x);
  // Add to tape
  std::vector<ad> X_(X, X + n);
  std::vector<ad> Y_;
  Y_ = TMBad::global::Complete<TMBad::LSEOp>(n, 0) (X_);
  // Pass to R
  size_t m = Y_.size();
  Rcpp::ComplexVector y(m);
  for (size_t j=0; j < m; j++) {
    y[j] = ad2cplx(Y_[j]);
  }
  return as_advector(y);
}
