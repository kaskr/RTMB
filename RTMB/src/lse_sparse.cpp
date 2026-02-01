// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

template<bool sparse=false>
struct LogSpaceSumOp : TMBad::global::DynamicOperator< -1 , 1 > {
  size_t n;
  static const bool have_input_size_output_size = true; // FIXME: Should give compile time error if 'false'
  static const bool add_forward_replay_copy = true;
  TMBad::Index input_size()  const { return this->n; }
  TMBad::Index output_size() const { return 1; }
  LogSpaceSumOp (size_t n) : n(n) { }
  void forward(TMBad::ForwardArgs<TMBad::Scalar> &args) {
    TMBad::Scalar Max = -INFINITY;
    for (size_t i=0; i<n; i++) {
      if (Max < args.x(i)) Max = args.x(i) ;
    }
    args.y(0) = 0;
    for (size_t i=0; i<n; i++) {
      args.y(0) += exp( args.x(i) - Max );
    }
    args.y(0) = Max + log(args.y(0));
  }
  template<class Type>
  void forward(TMBad::ForwardArgs<Type> &args) {
    Rcpp::stop("unsupported method");
  }
  // void forward(ForwardArgs<Replay> &args) {
  //   std::vector<ad_plain> x(input_size());
  //   for (TMBad::Index i=0; i<input_size(); i++)
  //     x[i] = args.x(i);
  //   args.y(0) = logspace_sum(x);
  // }
  template <class Type>
  Type SPDERIV(Type x) { return ( sparse ? Type(TMBad::sparseDeriv(x)) : x ); }
  template <class Type> void reverse(TMBad::ReverseArgs<Type> &args) {
    for (size_t i=0; i<n; i++) {
      args.dx(i) += SPDERIV ( exp( args.x(i) - args.y(0) ) * args.dy(0) );
    }
  }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) {
    Rcpp::stop("LSE source writer not available");
  }
  const char* op_name() {return sparse ? "LSE0" : "LSE"; }
};
// ad_plain logspace_sum (const std::vector<ad_plain> &x) {
//   OperatorPure* pOp = get_glob()->getOperator<LogSpaceSumOp>(x.size());
//   return get_glob()->add_to_stack<LogSpaceSumOp>(pOp, x)[0];
// }
// template<class T>
// T logspace_sum (const std::vector<T> &x_) {
//   std::vector<ad_plain> x(x_.begin(), x_.end());
//   return logspace_sum(x);
// }
                          
// [[Rcpp::export]]
ADrep LSE0(ADrep x) {
  std::vector<ad> x_(x.adptr(), x.adptr() + x.size());
  TMBad::global::Complete<LogSpaceSumOp<true> > op(x.size());
  std::vector<ad> y_ = op(x_);
  return ADrep(y_.data(), y_.data() + y_.size());
}
