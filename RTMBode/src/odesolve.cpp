#include <RcppEigen.h>

typedef void * (*FUN_PTR)(SEXP f, const Eigen::VectorXd &x, Eigen::VectorXd &y);

static SEXP PTR;
static Eigen::VectorXd X;
static Eigen::VectorXd Y;

extern "C"
SEXP set_pointers(SEXP ptr, SEXP x, SEXP y) {
  if (!R_ExternalPtrAddr(ptr)) Rf_error("oops!");
  PTR = ptr;
  X = Eigen::Map<Eigen::VectorXd> (REAL(x), LENGTH(x));
  Y = Eigen::Map<Eigen::VectorXd> (REAL(y), LENGTH(y));
  return R_NilValue;
}

/* Derivatives (calldef in RTMB.h) */
extern "C"
void desolve_derivs (int *neq, double *t, double *y, double *ydot,
                     double *yout, int *ip) {
  static FUN_PTR tmb_forward = (FUN_PTR) R_GetCCallable("RTMB", "tmb_forward");
  X[0] = *t;
  for (int i=0; i < *neq; i++) X[i+1] = y[i];
  tmb_forward(PTR, X, Y);
  for (int i=0; i < *neq; i++) ydot[i] = Y[i];
}

static const R_CallMethodDef CallEntries[] = {
    {"desolve_derivs", (DL_FUNC) &desolve_derivs, 6},
    {"set_pointers",   (DL_FUNC) &set_pointers, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_RTMBode(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
