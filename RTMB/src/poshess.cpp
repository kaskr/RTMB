// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

/* -------------------------------------------------------------------------- */
/* --- Positive Hessian ----------------------------------------------------- */
/* -------------------------------------------------------------------------- */

namespace TMBad {

// std::vector<ad_plain> PHE (const std::vector<ad_plain> &x, int order); // Forward declare
// template<class T>
// std::vector<T> PHE (const std::vector<T> &x_, int order); // Forward declare

struct PHEOp : global::DynamicOperator< -1 , 1 > {
  size_t n;
  int order; // 0 or 1
  std::shared_ptr<TMBad::ADFun<> > F;   // Function tape
  std::shared_ptr<TMBad::ADFun<> > G;   // Gradient tape
  std::shared_ptr<TMBad::ADFun<> > H2;  // *Alternative* Hessian tape
  static const bool have_input_size_output_size = true; // FIXME: Should give compile time error if 'false'
  Index input_size()  const { return this->n; }
  Index output_size() const {
    return (order==0 ? 1 : n);
  }
  PHEOp (TMBad::ADFun<> F, TMBad::ADFun<> G, TMBad::ADFun<> H2) :
    n(F.Domain()),
    order(0),
    F(std::make_shared<TMBad::ADFun<> >(F)),
    G(std::make_shared<TMBad::ADFun<> >(G)),
    H2(std::make_shared<TMBad::ADFun<> >(H2))
  {
    TMBAD_ASSERT(F.Range() == 1);
    TMBAD_ASSERT(F.Domain() == G.Domain());
    // if (F.Range() != 1) Rcpp::stop("oops");
    // if (F.Range() != 1) Rcpp::stop("oops");
  }
  void forward(ForwardArgs<Scalar> &args) {
    std::vector<Scalar> x = args.x_segment(0, n);
    if (order == 0) {
      args.y(0) = (*F)(x)[0];
    }
    if (order == 1) {
      args.y_segment(0, n) = (*F).Jacobian(x);
    }
  }
  // Forward replay => copy
  void forward(ForwardArgs<Replay> &args) {
    std::vector<Replay> x = args.x_segment(0, n);
    std::vector<Replay> y = TMBad::global::Complete<PHEOp>(*this)(x);
    args.y_segment(0, output_size()) = y;
  }
  //template <class Type> void reverse(ReverseArgs<Type> &args) {
  void reverse(ReverseArgs<Scalar> &args) {
    std::vector<Scalar> x  = args.x_segment(0, n);
    std::vector<Scalar> dy = args.dy_segment(0, output_size()); // n or 1    
    if (order == 0) {
      args.dx_segment(0, n) += (*F).Jacobian(x, dy); // length(dy)=1 length(F(x))=1
    }
    if (order == 1) {
      args.dx_segment(0, n) += (*G).Jacobian(x, dy); // length(dy)=n length(G(x))=n
    }
  }

  void reverse(ReverseArgs<Replay> &args) {
    if (order == 0) {
      PHEOp cpy(*this);
      cpy.order = 1;
      std::vector<Replay> x = args.x_segment(0, n);
      std::vector<Replay> g = TMBad::global::Complete<PHEOp>(cpy)(x);
      for (size_t i=0; i<n; i++) {
        args.dx(i) += g[i] * args.dy(0);
      }
    }
    if (order == 1) {
      if (newton::get_alternative_hessian()) {
        // Report fake derivatives
        std::vector<Replay> x = args.x_segment(0, n);
        std::vector<Replay> h2 = (*H2)(x);
        std::vector<Replay> w = args.dy_segment(0, n);
        for (size_t i=0; i<n*n; i++) {
          args.dx(i / n) += h2[i] * w[i % n];
        }
      } else {
        // Report real derivatives
        std::vector<Replay> x = args.x_segment(0, n);
        std::vector<Replay> w = args.dy_segment(0, n);
        args.dx_segment(0, n) += (*G).Jacobian(x, w);
      }
    }
  }

  
  template <class Type> void reverse(ReverseArgs<Type> &args) {
    Rcpp::stop("Not implemented");
  }
  
  const char* op_name() {return "PHEOp";}
};

// std::vector<ad_plain> PHE (const std::vector<ad_plain> &x, int order) {
//   OperatorPure* pOp = get_glob()->getOperator<PHEOp>(x.size(), order);
//   return get_glob()->add_to_stack<PHEOp>(pOp, x);
// }
// template<class T>
// std::vector<T> PHE (const std::vector<T> &x_, int order) {
//   std::vector<ad_plain> x(x_.begin(), x_.end());
//   std::vector<ad_plain> y = PHE(x, order);
//   std::vector<T> ans(y.begin(), y.end());
//   return ans;
// }

}

// [[Rcpp::export]]
Rcpp::ComplexVector PHE(Rcpp::XPtr<TMBad::ADFun<> > F,
                        Rcpp::XPtr<TMBad::ADFun<> > G,
                        Rcpp::XPtr<TMBad::ADFun<> > H2,
                        const Rcpp::ComplexVector x) {
  CHECK_INPUT(x);
  size_t n = x.size();
  ad* X = adptr(x);
  // Add to tape
  std::vector<ad> X_(X, X + n);
  std::vector<ad> Y_;
  Y_ = TMBad::global::Complete<TMBad::PHEOp>(*F, *G, *H2) (X_);
  // Pass to R
  size_t m = Y_.size();
  Rcpp::ComplexVector y(m);
  for (size_t j=0; j < m; j++) {
    y[j] = ad2cplx(Y_[j]);
  }
  return as_advector(y);
}
