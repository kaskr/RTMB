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

     FIXME: Many ignored details, such as
     - contiguous variable requirements (pointer ops, updating ops, output_size()>1)
     - ConstOp protection
*/
std::vector<Index> remap_values(global glob) {
  std::vector<Index> vr(glob.values.size()); // ans
  std::vector<Index> count(glob.values.size(), 0);
  struct stack_t {
    // Fictive container of n consecutive elements {0,...,n-1}.
    size_t n;
    // Previously used elements that are now free.
    std::vector<Index> avail;
    stack_t () : n(0) { }
    // Get new unoccupied index
    Index get_new_index(bool use_avail = true) {
      Index ans;
      if (use_avail && avail.size() > 0) {
        ans = avail.back();
        avail.pop_back();
      } else {
        ans = n++;
      }
      return ans;
    }
    void free_index(Index i) {
      avail.push_back(i);
    }
  } stack;
  struct increment_counter {
    std::vector<Index>& count;
    increment_counter(std::vector<Index>& count) : count(count) {}
    void operator()(Index var) {
      count[var]++;
    }
  } inc_counter(count);
  struct decrement_counter {
    std::vector<Index>& count;
    stack_t &stack;
    std::vector<Index>& vr;
    decrement_counter(std::vector<Index>& count,
                      stack_t &stack,
                      std::vector<Index>& vr) : count(count), stack(stack), vr(vr) {}
    void operator()(Index var) {
      count[var]--;
      if (count[var] == 0)
        stack.free_index(vr[var]);
    }
  } dec_counter(count, stack, vr);
  // Forward pass 1
  Args<> args(glob.inputs);
  Dependencies dep;
  // Example: how to mark infinite lifetime for InvOp
  // dep.Base::operator=(glob.inv_index);
  // dep.apply(inc_counter);
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
    if (glob.opstack[i]->output_size() != 1)
      Rcpp::stop("FIXME: This loop assumes 'output_size()==1' for all opstack[i]");
    // Get new index
    // If DepOp get a brand new index. Otherwise reuse available indices.
    bool isDepOp = glob.opstack[i] -> info().test(TMBad::op_info::dependent_variable);
    if (isDepOp) {
      vr[i] = stack.get_new_index(false); // Get brand new index, and don't free inputs
    } else {
      // Get inputs to this node
      dep.resize(0);
      glob.opstack[i]->dependencies(args, dep);
      // remove inputs
      dep.apply(dec_counter);
      vr[i] = stack.get_new_index(true); // Re-use available
    }
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

// Patched souce code writers
// - We need a bit more control with it
namespace TMBad {
void rtmb_write_forward(global &glob, code_config cfg = code_config()) {
  using std::setw; using std::left; using std::endl;
  std::ostream& cout = *cfg.cout;
  cfg.write_header_comment();
  std::string type = (cfg.gpu ? "access" : "double*");
  std::string void_str = cfg.void_str();
  cfg.gpu=false;
  cout << void_str << " forward(" << type << " v) {" << endl;
  cfg.init_code();
  ForwardArgs<Writer> args(glob.inputs, glob.values);
  for (size_t i=0; i<glob.opstack.size(); i++) {
    std::ostringstream buffer;
    Writer::cout = &buffer;
    glob.opstack[i]->forward(args); // FIXME: pass args copy
    write_common(buffer, cfg, i);
    glob.opstack[i]->increment(args.ptr);
  }
  cout << "}" << endl;
}
void rtmb_write_reverse(global &glob, code_config cfg = code_config()) {
  using std::setw; using std::left; using std::endl;
  std::ostream& cout = *cfg.cout;
  cfg.write_header_comment();
  std::string type = (cfg.gpu ? "access" : "double*");
  std::string void_str = cfg.void_str();
  cfg.gpu=false;
  cout << void_str << " reverse(" << type << " v, " << type << " d) {" << endl;
  cfg.init_code();
  ReverseArgs<Writer> args(glob.inputs, glob.values);
  for (size_t i=glob.opstack.size(); i>0; ) {
    i--;
    glob.opstack[i]->decrement(args.ptr);
    std::ostringstream buffer;
    Writer::cout = &buffer;
    glob.opstack[i]->reverse(args); // FIXME: pass args copy
    write_common(buffer, cfg, i);
  }
  cout << "}" << endl;
}
}

// [[Rcpp::export]]
void src_transform(Rcpp::XPtr<TMBad::ADFun<> > adf, Rcpp::List config) {
  TMBad::code_config cfg;
  cfg.gpu = config["gpu"];
  cfg.asm_comments = false;
  cfg.cout = &Rcout;
  TMBad::global glob = adf->glob; // Invoke deep copy
  TMBad::compress(glob);
  rtmb_write_forward(glob, cfg);
  rtmb_write_reverse(glob, cfg);
}

// [[Rcpp::export]]
void reorder_depth_first(Rcpp::XPtr<TMBad::ADFun<> > adf) {
  TMBad::reorder_depth_first(adf->glob);
}
