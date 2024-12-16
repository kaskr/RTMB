// Modified checkpoint operator
#include "RTMB.h"

namespace TMBad {

/* Copy-pasted from TMBad, then renamed and tweaked */
template<class DerivativeTable>
struct RTMB_AtomOp : global::DynamicOperator< -1, -1> {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;

  // derivatives table
  TMBAD_SHARED_PTR<DerivativeTable> dtab;

  // Order of this instance
  int order;

  // CTOR from derivative table ctor
  template<class T1>
  RTMB_AtomOp(const T1 &F) :
    dtab(std::make_shared<DerivativeTable>(F)),
    order(0) { }
  template<class T1, class T2>
  RTMB_AtomOp(const T1 &F, const T2 &x) :
    dtab(std::make_shared<DerivativeTable>(F, x)),
    order(0) { }
  template<class T1, class T2, class T3>
  RTMB_AtomOp(const T1 &F, const T2 &x, const T3 &t) :
    dtab(std::make_shared<DerivativeTable>(F, x, t)),
    order(0) { }

  Index input_size()  const {  return (*dtab)[order].Domain(); }
  Index output_size() const {  return (*dtab)[order].Range(); }

  void forward(ForwardArgs<Scalar> &args) {
    // Conditionally retape
    (*dtab).retape(args);
    // Make sure order is available
    (*dtab).requireOrder(order);
    // Number of inputs / outputs for this order
    size_t n = input_size();
    size_t m = output_size();
    // Get vector of inputs
    auto x = args.x_segment(0, n);
    // Eval
    args.y_segment(0, m) = (*dtab)[order](x);
  }

  void reverse(ReverseArgs<Scalar> &args) {
    // Number of inputs / outputs for this order
    size_t n = input_size();
    size_t m = output_size();
    // Get vector of inputs
    auto x = args.x_segment(0, n);
    auto w = args.dy_segment(0, m);
    // Eval
    args.dx_segment(0, n) += (*dtab)[order].Jacobian(x, w);
  }

  void reverse(ReverseArgs<global::Replay> &args) {
    // Number of inputs / outputs for this order
    size_t n = input_size();
    size_t m = output_size();
    // Concatenation (x, dy)
    std::vector<global::Replay> x = args.x_segment(0, n);
    if (DerivativeTable::packed) x = repack(x);
    std::vector<global::Replay> w = args.dy_segment(0, m);
    std::vector<global::Replay> xw;
    xw.insert(xw.end(), x.begin(), x.end());
    xw.insert(xw.end(), w.begin(), w.end());
    // Eval
    (*dtab).requireOrder(order + 1);
    RTMB_AtomOp cpy(*this);
    cpy.order++;
    args.dx_segment(0, n) += global::Complete<RTMB_AtomOp>(cpy)(xw);
  }

  // Not yet implemented
  // void forward(ForwardArgs<Writer> &args) { ASSERT(false); }
  template<class T>
  void forward(ForwardArgs<T> &args) { ASSERT(false); }
  void reverse(ReverseArgs<Writer> &args) { ASSERT(false); }

  const char* op_name() { return "AtomOp"; }

  void print(global::print_config cfg) {
    std::cout << cfg.prefix;
    std::cout << "order=" << order << " ";
    std::cout << "(*dtab).size()=" << (*dtab).size() << " ";
    std::cout << "dtab=" << &(*dtab) << "\n";
    (*dtab)[order].print(cfg);
  }

  // ------------ Tweaks

  // Boolean. Use defaults with optimization tweak for packed case
  static const bool have_forward_mark_reverse_mark = true;
  void forward(ForwardArgs<bool>& args) {
    args.mark_dense(*this);
  }
  void reverse(ReverseArgs<bool>& args) {
    if (DerivativeTable::packed) {
      for (size_t i=0; i<output_size(); i++) std::cout << args.y(i);
      std::cout << "\n";
    }
    args.mark_dense(*this);
  }
  // TMBad FIXME: Ideally shouldn't have to add this...
  static const bool have_dependencies = true;
  void dependencies(Args<> &args, Dependencies &dep) const {
    // 'dep' is empty on input and contains the result on output
    Index ninput_ = this->input_size();
    for (Index j=0; j<ninput_; j++) dep.push_back(args.input(j));
  }
};

}

// 'atomic_transform' helpers
// Pack single argument
std::vector<ad> pack1(std::vector<ad> x) {
  TMBad::ad_segment xs(x.data(), x.size());
  TMBad::ad_segment xp = TMBad::pack(xs);
  return TMBad::concat(std::vector<TMBad::ad_segment>(1, xp));
}
// unpack single argument
std::vector<ad> unpack1(std::vector<ad> x) {
  TMBad::ad_segment xs(x.data(), x.size());
  TMBad::ad_segment x0 = TMBad::unpack(xs);
  return TMBad::concat(std::vector<TMBad::ad_segment>({x0}));
}
// unpack two arguments
std::vector<ad> unpack2(std::vector<ad> x) {
  TMBad::ad_segment x0 = TMBad::unpack(x, 0);
  TMBad::ad_segment x1 = TMBad::unpack(x, 1);
  return TMBad::concat(std::vector<TMBad::ad_segment>({x0, x1}));
}
// Pack inputs/outputs of a tape
struct PackedTape {
  TMBad::ADFun<>* adf;
  PackedTape(TMBad::ADFun<>* adf) : adf(adf) { }
  std::vector<ad> operator() (const std::vector<ad> &x) {
    return pack1((*adf)(unpack2(x)));
  }
};

void rtmb_atomic_transform(TMBad::ADFun<>* adf) {
  TMBad::ADFun<> F;
  std::vector<double> xd = (*adf).DomainVec();
  // resolve refs
  std::vector<ad> outer_vars = (*adf).resolve_refs();
  // No outer vars?
  if (outer_vars.size() == 0) {
    *adf = (*adf).atomic();
    return;
  }
  TMBad::forceContiguous(outer_vars);
  // Start new context
  F.glob.ad_start();
  std::vector<ad> xad(xd.begin(), xd.end());
  TMBad::Independent(xad);
  // Copy outer vars to new context (adds RefOps)
  for (size_t i=0; i<outer_vars.size(); i++)
    outer_vars[i] = outer_vars[i].copy();
  // pack inner params
  std::vector<ad> xad_pack = pack1(xad);
  std::vector<ad> outer_vars_pack = pack1(outer_vars);
  xad_pack.insert(xad_pack.end(), outer_vars_pack.begin(), outer_vars_pack.end());
  TMBad::ADFun<> Tape(PackedTape(adf), xad_pack);
  typedef TMBad::standard_derivative_table< TMBad::ADFun<>, /*packed*/ true > DTab;
  TMBad::global::Complete<TMBad::RTMB_AtomOp<DTab> > Fatom(Tape);
  std::vector<ad> yp = Fatom(xad_pack);
  std::vector<ad> y = unpack1(yp);
  TMBad::Dependent(y);
  F.glob.ad_stop();
  *adf = F;
}
