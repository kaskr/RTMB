// Modified checkpoint operator
#include "RTMB.h"

namespace TMBad {

typedef standard_derivative_table< ADFun<>, /*packed*/ true > RTMB_DTab;

template<class DerivativeTable>
struct RTMB_AtomOp : AtomOp<DerivativeTable> {
  typedef AtomOp<DerivativeTable> Base;
  INHERIT_CTOR(RTMB_AtomOp, Base)

  template<class Type>
  void forward(ForwardArgs<Type>& args) { Base::forward(args); }
  template<class Type>
  void reverse(ReverseArgs<Type>& args) { Base::reverse(args); }
    
  // Boolean. Use defaults with optimization tweak for packed case
  static const bool have_forward_mark_reverse_mark = true;
  void forward(ForwardArgs<bool>& args) {
    args.mark_dense(*this);
  }
  void reverse(ReverseArgs<bool>& args) {
    if (DerivativeTable::packed) {
      for (size_t i=0; i<Base::output_size(); i++) std::cout << args.y(i);
      std::cout << "\n";
    }
    args.mark_dense(*this);
  }
  // TMBad FIXME: Ideally, shouldn't have to add this member...
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
