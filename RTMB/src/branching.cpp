// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

namespace TMBad {

/** \brief Fixed derivative table used by `AtomOp` */
template<class ADFun, bool packed_ = false>
struct derivative_table : std::vector<ADFun> {
  static const bool packed = packed_;
  /** \brief Add derivatives up to this order. */
  void requireOrder(size_t n) {
    while ( (*this).size() <= n) {
      (*this).push_back((*this).back().JacFun());
      (*this).back().optimize();
    }
  }
  /** \brief Retaping this derivative table has no effect. */
  template<class ARGS>
  void retape(ARGS &args) { }
  /** \brief Set zero order function of this derivative table. */
  derivative_table(const ADFun &F) :
    std::vector<ADFun>(1, F) { }
};

struct BranchOp : global::DynamicOperator< -1, -1> {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  typedef AtomOp<derivative_table<ADFun<> > > Atomic;
  Atomic f;
  Atomic g; // same Domain and Range as f
  BranchOp(Atomic f, Atomic g) : f(f), g(g) {}
  Index input_size()  const {  return f.input_size() + 1; }
  Index output_size() const {  return f.output_size(); }
  void forward(ForwardArgs<Scalar> &args) {
    bool select_f = (args.x(0) != 0);
    args.ptr.first++; // Skip first input for f and g
    if (select_f) {
      f.forward(args);
    } else {
      g.forward(args);
    }
    args.ptr.first--; // Restore
  }
  void reverse(ReverseArgs<Scalar> &args) {
    bool select_f = (args.x(0) != 0);
    args.ptr.first++; // Skip first input for f and g
    if (select_f) {
      f.reverse(args);
    } else {
      g.reverse(args);
    }
    args.ptr.first--; // Restore
  }
  void reverse(ReverseArgs<global::Replay> &args) {
    // Number of inputs / outputs for this order
    size_t n = input_size();
    size_t m = output_size();
    // Concatenation (x, dy)
    std::vector<global::Replay> x = args.x_segment(0, n);
    std::vector<global::Replay> w = args.dy_segment(0, m);
    // Eval
    BranchOp cpy(*this);
    cpy.f.order++;
    cpy.g.order++;
    (cpy.f.dtab)->requireOrder(cpy.f.order);
    (cpy.g.dtab)->requireOrder(cpy.g.order);
    // Jacobian **by-row**
    n--;
    std::vector<global::Replay> J = global::Complete<BranchOp>(cpy)(x);
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        args.dx(j+1) += w[i] * J[j +  i * n];
      }
    }
  }
  // Not yet implemented
  template<class T>
  void forward(ForwardArgs<T> &args) { ASSERT(false); }
  void reverse(ReverseArgs<Writer> &args) { ASSERT(false); }

  // Any of the two 'dtab' can be used to uniqely identify this operator.
  // However, unlike AtomOp, we must assume that identifier is order
  // dependent (because input/output dimension doesn't necessarily
  // change for higher order JacFun() unlike WgtJacFun() ).
  static const bool have_custom_identifier = true;
  void* custom_identifier() { return &((*(f.dtab))[f.order]); }

  const char* op_name() { return "BranchOp"; }

  void print(global::print_config cfg) {
    // f
    Rcout << cfg.prefix;
    Rcout << "Left branch: order=" << f.order << " ";
    Rcout << "(*dtab).size()=" << (f.dtab)->size() << " ";
    Rcout << "dtab=" << &(*(f.dtab)) << "\n";
    (*(f.dtab))[f.order].print(cfg);
    // g
    Rcout << cfg.prefix;
    Rcout << "Right branch: order=" << g.order << " ";
    Rcout << "(*dtab).size()=" << (g.dtab)->size() << " ";
    Rcout << "dtab=" << &(*(g.dtab)) << "\n";
    (*(g.dtab))[g.order].print(cfg);
  }

};

}


// [[Rcpp::export]]
ADrep Branch(Rcpp::XPtr<TMBad::ADFun<> > f, Rcpp::XPtr<TMBad::ADFun<> > g, ADrep x) {
  CHECK_INPUT(x);
  size_t n = x.size();
  TMBad::global::Complete<TMBad::BranchOp> F(*f, *g);
  ad* px = adptr(x);
  std::vector<ad> X(px, px + n);
  std::vector<ad> Y = F(X);
  return ADrep (Y.data(), Y.data() + Y.size());
}
