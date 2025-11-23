#include "RTMB.h"

// [[Rcpp::export]]
ADrep bisect_atom(Rcpp::XPtr<TMBad::ADFun<> > adf, ADrep x_, Rcpp::List cfg) {
  struct bisect {
    TMBad::ADFun<> f;
    // Configurable absolute error tolerance
    double tol;
    // Test if a term is known to machine tolerance, but relaxed a bit just in case.
    double machine_tolerance;
    // Test if Wynn series has converged
    double wynn_convergence_tolerance;
    // Configurable max number of subdivisions
    int subdivisions;
    // Configurable flag
    bool stop_on_error;
    // CTOR
    bisect(Rcpp::XPtr<TMBad::ADFun<> > adf, Rcpp::List cfg) {
      f = *adf;
      tol = Rcpp::NumericVector(cfg["abs.tol"])[0];
      machine_tolerance = 1e-14;
      wynn_convergence_tolerance = 1e-10;
      subdivisions = Rcpp::IntegerVector(cfg["subdivisions"])[0];
      stop_on_error = Rcpp::LogicalVector(cfg["stop.on.error"])[0];
    }
    // Current number of subdivisions
    int n;
    // Input:  x = { a, b, user parameters }
    // Output: y = { integrate(f, a, b) , absolute error, subdivisions }
    std::vector<ad> operator()(const std::vector<ad> &x) {
      this->n = 1;
      double tol = this->tol;
      // Initialize 1st eval
      std::vector<ad> y = f(x);
      ad result = y[0];
      ad abserr = asDouble(y[1]);
      if (abserr > tol) {
        subdiv(x, y, tol, true, true);
      }
      y.push_back(this->n);
      if (stop_on_error) {
        if (this->n >= subdivisions) {
          TMBad::get_glob()->ad_stop();
          Rcpp::stop("Max subdivisions exceeded");
        }
      }
      return y;
    }
    // Input:  x = { a, b, user parameters }
    // Output: integrate(f, a, b) satisfying abs tolerance 'tol'
    // NOTE: Argument 'y' is overwritten with output!
    void subdiv(const std::vector<ad> &x,
                std::vector<ad> &y,
                double tol,
                bool left_bound,  // Does sub-interval contain left endpoint of original interval?
                bool right_bound  // Does sub-interval contain right endpoint of original interval?
                ) {
      if (n >= subdivisions) return;
      // One more subdivision
      n++;
      // Mid point
      ad mid = .5 * (asDouble(x[0]) + asDouble(x[1]));
      // Interval empty to machine tolerance
      if (asDouble(x[0]) == asDouble(mid) || asDouble(x[1]) == asDouble(mid)) {
        y[0] = 0;
        y[1] = 0;
        return;
      }
      // Left interval
      std::vector<ad> xl(x);
      xl[1] = mid;
      // Right interval
      std::vector<ad> xr(x);
      xr[0] = mid;
      // Eval
      std::vector<ad> yl = f(xl);
      std::vector<ad> yr = f(xr);
      // Abs errors
      double el = asDouble(yl[1]);
      double er = asDouble(yr[1]);
      // Handle left bound singularity
      if (left_bound && !right_bound) { // Approaching left limit
        if (er < machine_tolerance && el > tol) {   // Is error of right half negligible?
          y = wynn_extrapolate(xl, yr, tol, left_bound, right_bound);
          return;
        }
      }
      // Handle right bound singularity
      if (!left_bound && right_bound) { // Approaching right limit
        if (el < machine_tolerance && er > tol) {   // Is error of left half negligible?
          y = wynn_extrapolate(xr, yl, tol, left_bound, right_bound);
          return;
        }
      }
      // Continure subdividing?
      if (el + er > tol) {
        // Subdivide largest error first
        if (el > er) {
          double ratio = el / (el+er);
          subdiv(xl, yl, tol * ratio, left_bound, false);
          double err = asDouble(yl[1]);
          if (err + asDouble(yr[1]) > tol)
            subdiv(xr, yr, tol - err, false, right_bound);
        } else {
          double ratio = er / (el+er);
          subdiv(xr, yr, tol * ratio, false, right_bound);
          double err = asDouble(yr[1]);
          if (err + asDouble(yl[1]) > tol)
            subdiv(xl, yl, tol - err, left_bound, false);
        }
      }
      // Add up
      y = { yl[0] + yr[0],
            yl[1] + yr[1] };
    }
    // Extrapolation toward interval bounds (singular case)
    std::vector<ad> wynn_extrapolate(std::vector<ad> x,
                                     std::vector<ad> y,
                                     double tol,
                                     bool left_bound,  // Does sub-interval contain left endpoint of original interval?
                                     bool right_bound  // Does sub-interval contain right endpoint of original interval?
                                     ) {
      std::vector<ad> eps;
      eps = wynn_update(eps, y[0]);
      for (size_t i=0; i<20; i++) {
        // Mid point
        ad mid = .5 * (asDouble(x[0]) + asDouble(x[1]));
        // Left interval
        std::vector<ad> xl(x);
        xl[1] = mid;
        // Right interval
        std::vector<ad> xr(x);
        xr[0] = mid;
        // Eval
        if (left_bound)  y = f(xr);
        if (right_bound) y = f(xl);
        // Must know term to machine tolerance
        if (asDouble(y[1]) > machine_tolerance) { // Fail => Reject extraplation
          // Could now be that:
          // 1. This term really is imprecise due to a difficult integrand. Then we must drop the extrapolation.
          // 2. OR, the term is good enough, but the reported abs error y[1] is inaccurate. Then we should carry on extrapolating.
          if (false) { // This would seem like the right thing to do if in case 1:
            // Skip this term and fall back on normal sub division for remaining terms
            subdiv(x, y, tol, left_bound, right_bound);
            // First element of eps table holds the plain sum of terms known without error
            y[0] += eps[0];
            // Done
            return y;
          }
          // Hard to tell if (1) or (2), so we report NaN
          if (stop_on_error) {
            TMBad::get_glob()->ad_stop();
            Rcpp::stop("Extrapolation failed");
          }
          y[0] = R_NaN;
          y[1] = R_NaN;
          return y;
        }
        eps = wynn_update(eps, y[0]);
        if (wynn_converged(eps))
          break;
        // Set next interval
        if (left_bound)  x = xl;
        if (right_bound) x = xr;
      }
      std::vector<ad> ans(2);
      ans[0] = eps.back();
      ans[1] = 0;
      return ans;
    }
    // Wynn's series acceleration: Update with new term
    std::vector<ad> wynn_update(std::vector<ad> eps, ad term) {
      if (eps.size() == 0)
        return std::vector<ad> (1, term);
      std::vector<ad> y(eps.size() + 1);
      y[0] = eps[0] + term;
      for (size_t i=1; i < y.size(); i++) {
        if (i==1)
          y[i] = 1. / term;
        else
          y[i] = eps[i-2] + 1. / (y[i-1] - eps[i-1]);
      }
      return y;
    }
    // Check Wynn series for convergence
    bool wynn_converged(const std::vector<ad> &eps) {
      size_t n = eps.size();
      if (n % 2 == 0) return false;
      if (n < 3) return false;
      return (std::abs(asDouble(eps[n-1]) - asDouble(eps[n-3])) < wynn_convergence_tolerance);
    }
  } integrate(adf, cfg);
  std::vector<ad> x(x_.adptr(), x_.adptr() + x_.size());
  TMBad::ADFun<> A = TMBad::ADFun_retaping(integrate, x);
  std::vector<ad> y = A(x);
  ADrep ans(y.data(), y.data() + y.size());
  return ans;
}
