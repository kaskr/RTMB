// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

namespace TMBad {

/* Construct a 'value remapping vector' (vr) such that the long vector of `values` can be accessed during a forward pass via the short version `short_values`.

   `values[i] = short_values[vr[i]]`

   The purpose is to optimize memory locality for source code written to the GPU (applicable to forward pass only!).
   Note that `vr[i]` is the ith values position in the 'register' during its entire lifetime.
   Once `value[i]` is not needed anymore, its position can be marked 'available' and assigned to a new variable.

   Forward pass 1:
   - Count the number of times variable 'i' is used as input.
   - Add extra 1-count to selected operators (ConstOp and probably InvOp) to give them infinite lifetime.

   Forward pass 2:
   - Assign new 'stack index' for variable 'i': vr[i] = get_new_index()
   - Remove inputs that are no longer needed:
     count[op.dependencies()] -= 1;
     elim = which(count[op.dependencies()]==0);
     free_index(vr[elim]); // Free these registers
*/
std::vector<Index> remap_values(global glob) {
  std::vector<Index> vr(glob.values.size()); // ans
  std::vector<Index> count(glob.values.size());
  struct increment_counter {
    std::vector<Index>& count;
    increment_counter(std::vector<Index>& count) : count(count) {}
    void operator()(Index var) {
      count[var]++;
    }
  } inc_counter(count);
  struct stack_t {
    size_t n;
    std::vector<bool> b;
    std::vector<Index> avail; // Previously used - now free
    stack_t ( size_t max_capacity ) : n(0) { b.resize(max_capacity, false); }
    Index get_new_index() {
      if (false) {
        std::cout << "get_new_index() ";
        std::cout << "b=" << b << " ";
        std::cout << "avail=" << avail << " ";
        std::cout << "n=" << n;
        std::cout << "\n";
      }
      n++;
      Index ans;
      if (avail.size() > 0) {
        ans = avail.back();
        avail.pop_back();
      } else {
        // avail is empty
        ans = n-1;
      }
      b[ans] = true;
      return ans;
    }
    void free_index(Index i) {
      n--;
      b[i] = false;
      avail.push_back(i);
    }
  } stack(glob.values.size());
  struct decrement_counter {
    std::vector<Index>& count;
    stack_t &stack;
    std::vector<Index>& vr;
    decrement_counter(std::vector<Index>& count, stack_t &stack, std::vector<Index>& vr) : count(count), stack(stack), vr(vr) {}
    void operator()(Index var) {
      count[var]--;
      if (count[var] == 0)
        stack.free_index(vr[var]);
    }
  } dec_counter(count, stack, vr);
  // Forward pass 1
  Args<> args(glob.inputs);
  Dependencies dep;
  for (size_t i=0; i<glob.opstack.size(); i++) {
    // Add up inputs to this node
    dep.resize(0);
    glob.opstack[i]->dependencies(args, dep);
    // count inputs
    dep.apply(inc_counter);
    // Increment pointers
    glob.opstack[i]->increment(args.ptr);
  }
  // Forward pass 2
  args = Args<>(glob.inputs);
  for (size_t i=0; i<glob.opstack.size(); i++) {
    // Get inputs to this node
    dep.resize(0);
    glob.opstack[i]->dependencies(args, dep);
    // remove inputs
    dep.apply(dec_counter);
    // Get new index
    vr[i] = stack.get_new_index();
    // Increment pointers
    glob.opstack[i]->increment(args.ptr);
  }
  return vr;
}

}


// [[Rcpp::export]]
Rcpp::IntegerVector remap_values(Rcpp::XPtr<TMBad::ADFun<> > adf) {
  std::vector<TMBad::Index> rv = TMBad::remap_values(adf->glob);
  return Rcpp::IntegerVector(rv.data(), rv.data() + rv.size());
}
