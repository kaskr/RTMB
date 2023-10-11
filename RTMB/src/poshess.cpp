// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

/* -------------------------------------------------------------------------- */
/* --- POSDEF Term tag ------------------------------------------------------ */
/* -------------------------------------------------------------------------- */

namespace TMBad {

template<int order=0, bool setZero=false>
struct TermOp : global::Operator< 1 , 1 > {
  static const bool dynamic = true; // Prevent tape optimizer from remapping it
  template<class Type> void forward(TMBad::ForwardArgs<Type> &args) {
    if (!setZero)
      args.y(0) = (*this)(args.x(0));
    else
      args.y(0) = Type(0);
  }
  template<class Type> void reverse(TMBad::ReverseArgs<Type> &args) {
    if (order == 0)
      args.dx(0) += TermOp<1>()(args.dy(0));
    else
      args.dx(0) += args.dy(0);
  }
  const char* op_name() { return (order==0 ? "TermOp0" : "TermOp1") ; }
  ad operator()(ad x) {
    std::vector<ad_plain> xv(1,x);
    return global::Complete<TermOp>()(xv)[0];
  }
  double operator()(double x) {
    return x;
  }
  TMBad::Writer operator()(TMBad::Writer x) {
    return x;
  }
};

// Replay persistent InvOp
struct InvOp_ : global::InvOp {
  static const bool add_forward_replay_copy = true;
  const char* op_name() {return "InvOp_";}
};

}

// [[Rcpp::export]]
SEXP Term(const SEXP x_) {
  if (Rf_isNumeric(x_))
    return x_;
  if (!ad_context())
    return x_;
  Rcpp::ComplexVector x(x_);
  CHECK_INPUT(x);
  size_t n = x.size();
  ad* X = adptr(x);
  TMBad::TermOp<> F;
  Rcpp::ComplexVector y(n);
  for (size_t j=0; j < n; j++) {
    y[j] = ad2cplx(F(X[j]));
  }
  return as_advector(y);
}

// [[Rcpp::export]]
void TermsZero(Rcpp::XPtr<TMBad::ADFun<> > adf, bool setZero) {
  std::vector<TMBad::Index> nodes = find_op_by_name(adf->glob, "TermOp1");
  for (size_t i=0; i<nodes.size(); i++) {
    TMBad::OperatorPure* op;
    if (setZero)
      op = new TMBad::global::Complete<TMBad::TermOp<1, true> >();
    else
      op = new TMBad::global::Complete<TMBad::TermOp<1, false> >();
    std::swap(adf->glob.opstack[nodes[i]], op);
    op->deallocate();
  }
}

// [[Rcpp::export]]
void InvPersistent(Rcpp::XPtr<TMBad::ADFun<> > adf, bool setPers) {
  //std::vector<TMBad::Index> nodes = find_op_by_name(adf->glob, "InvOp");
  TMBad::OperatorPure* invop  = adf->glob.template getOperator<TMBad::global::InvOp>();
  TMBad::OperatorPure* invop_ = adf->glob.template getOperator<TMBad::InvOp_>();
  for (size_t i=0; i < adf->glob.opstack.size(); i++) {
    // Finds 'InvOp' or 'InvOp_' (The only operators with 'op_info::independent_variable' set)
    bool test = adf -> glob.opstack[i] -> info().test(TMBad::op_info::independent_variable);
    if (test) {
      if (setPers)
        adf -> glob.opstack[i] = invop_;
      else
        adf -> glob.opstack[i] = invop;
    }
  }
}
