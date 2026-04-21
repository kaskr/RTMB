#include "RTMB.h"

/* Patched version of TMB's sparse matrix exponential

   - Extends with rescaling option to avoid overflow/underflow.
   - Adds extra configuration member 'rescale_freq' to control how often to rescale.
   - Adds operator to calculate L1 norm used to rescale. Derivatives are not needed.
*/

namespace TMBad {
/* Add operator to evaluate L1 norm without derivatives */
struct VL1Op : global::DynamicOperator< 1 , 1 > {
  size_t n;
  VL1Op (size_t n) : n(n) { }
  template<class Type> void forward(ForwardArgs<Type> &args) {
    const Type* x = args.x_ptr(0);
    Type& y = args.y(0);
    y = 0;
    for (size_t i=0; i<n; i++) y += fabs(x[i]);
  }
  template <class Type> void reverse(ReverseArgs<Type> &args) { }
  // ---- Dependencies ----
  void dependencies(Args<> &args, Dependencies &dep) const {
    dep.add_segment(args.input(0), n);
  }
  static const bool have_dependencies = true;
  /** \brief This operator **has** implicit dependencies */
  static const bool implicit_dependencies = true;
  /** \brief It is **not* safe to remap the inputs of this operator */
  static const bool allow_remap = false; // FIXME: Deprecated?
  void forward(ForwardArgs<Writer> &args) { ASSERT(false); }
  void reverse(ReverseArgs<Writer> &args) { ASSERT(false); }
  const char* op_name() {return "VL1Op";}
};
ad_aug L1_norm(ad_segment x) {
  global::Complete<VL1Op> F(x.size());
  return F(x)[0];
}
} // end namespace TMBad

namespace sparse_matrix_exponential {
/* Extend original config with new member */
template<class T>
struct config2 : config<T> {
  typedef config<T> Base;
  /** \brief Rescale frequency */
  int rescale_freq;
  config2() : Base(), rescale_freq(Base::Nmax) {}
};
/* Modified matrix exponential with rescaling */
template <class T>
struct expm_series_rescaled {
  typedef vectorize::vector<T> vec;
  /** Number of terms */
  T N;
  /** Scale by exp(C) */
  T C;
  /** Vector of sparse matrix non zeros */
  vec A_values;
  /** Operator that multiplies sparse matrix with vector */
  TMBad::global::Complete<SpAxOp<T> > multiply;
  /** ADFun object that holds the entire matexp tape */
  TMBad::ADFun_packed<> F;
  /** Configuration */
  config2<T> cfg;
  /** Helper to update generator */
  void update(Eigen::SparseMatrix<T> &A) {
    // FIXME: Assert same pattern
    A_values = vec(A.valuePtr(), A.nonZeros());
  }
  expm_series_rescaled() {}
  expm_series_rescaled(Eigen::SparseMatrix<T> &A, T N, T C, config2<T> cfg=config2<T>()) :
    N(N), C(C), A_values(A.valuePtr(), A.nonZeros()), multiply(A), cfg(cfg)
  { }
  /** \brief Evaluate `x^T*exp(A)` */
  vec operator()(vec x) {
    N = min(N, T(cfg.Nmax));
    std::vector<TMBad::ad_segment> args = {A_values, x, vec(N, 1), vec(C, 1)};
    if (! F.initialized() ) {
      struct Test {
        config2<T> cfg;
        TMBad::Scalar Nold;
        bool operator() (const std::vector<TMBad::Scalar*> &x) {
          using TMBad::operator<<;
          TMBad::Scalar N = x[2][0];
          if ( (int) N == cfg.Nmax) {
            if (cfg.warn)
              Rf_warning("expm: N terms reduced to Nmax (%i)", (int) cfg.Nmax);
          }
          bool change = (Nold != N);
          if (cfg.trace && change) {
            Rcout << "Retaping:" << " Nold=" << Nold << " Nnew=" << N << "\n";
            Nold = N;
          }
          return change;
        }
      };
      Test N_changed = {cfg, N.Value()};
      F = TMBad::ADFun_retaping(*this, args, N_changed);
    }
    return F(args);
  }
private:
  friend class TMBad::PackWrap<expm_series_rescaled>;
  // Packed: (A, x) -> exp(A) * x
  TMBad::ad_segment operator() (const std::vector<TMBad::ad_segment> &args) {
    // Unpack input
    vec A = args[0];
    vec x = args[1];
    vec N_= args[2];
    vec C = args[3];
    int N = (int) N_[0].Value();
    // Evaluate series
    vec term(x), y(x);
    T log_scale = 0;
    for (int n=1; n<N; n++) {
      term = multiply(A, term) / n;
      y += term;
      if (!(n % cfg.rescale_freq)) {
        // rescale all previous terms
        T s = L1_norm(y);
        y = y / s;
        term = term / s;
        log_scale += log(s);
      }
    }
    y = y * exp(C + log_scale);
    return y;
  }
};
} // End namespace sparse_matrix_exponential

// [[Rcpp::export]]
ADrep expATv (Rcpp::RObject AT,
              ADrep v,
              ADrep N,
              ADrep C,
              Rcpp::List cfg,
              Rcpp::RObject cache) {
  if (!is_adsparse(AT)) Rcpp::stop("Expecting adsparse 'AT'");
  if (!is_adscalar(N)) Rcpp::stop("Expecting adscalar 'N'");
  // Inputs
  Eigen::SparseMatrix<ad> AT_ = SparseInput(AT);
  matrix<ad> v_ = MatrixInput(v);
  ad N_ = ScalarInput(N);
  ad C_ = ScalarInput(C);
  // Configuration parameters
  sparse_matrix_exponential::config2<ad> cfg_;
#define SET_CONFIG(XXX) if (!Rf_isNull(cfg[#XXX]))      \
  cfg_.XXX = Rcpp::IntegerVector((SEXP) cfg[#XXX])[0]
  SET_CONFIG(Nmax);
  SET_CONFIG(trace);
  SET_CONFIG(warn);
  SET_CONFIG(rescale_freq);
#undef SET_CONFIG
  // Set or get cached 'expm_series' object
  typedef sparse_matrix_exponential::expm_series_rescaled<ad> expm_t;
  expm_t* F;
  if (!cache.hasAttribute("expm_series")) {
    F = new expm_t(AT_, N_, C_, cfg_);
    SEXP ptr = Rcpp::XPtr<expm_t>(F);
    if (!Rf_isNull(cache))
      cache.attr("expm_series") = ptr;
  } else {
    SEXP ptr = cache.attr("expm_series");
    F = Rcpp::XPtr<expm_t> (ptr);
    F->update(AT_);
    F->N = N_;
    F->C = C_;
    F->cfg = cfg_;
  }
  // Evaluate
  matrix<ad> ans(v_.rows(), v_.cols());
  for (int j=0; j<ans.cols(); j++) {
    vector<ad> vec = v_.col(j).array();
    vector<ad> out = (*F)(vec);
    ans.col(j).array() = out;
  }
  return MatrixOutput(ans);
}
