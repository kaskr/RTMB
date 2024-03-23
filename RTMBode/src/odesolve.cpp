#include <RcppEigen.h>

typedef void * (*FUN_PTR)(void*);

static void* F;
static double* X;
static double* Y;
static FUN_PTR run_forward;

extern "C"
SEXP set_pointers(SEXP tape, SEXP x, SEXP y) {
  F = R_ExternalPtrAddr(tape);
  X = (double*) R_ExternalPtrAddr(x);
  Y = (double*) R_ExternalPtrAddr(y);
  run_forward = (FUN_PTR) R_GetCCallable("RTMB", "ptr_forward");
  return R_NilValue;
}

/* Derivatives (calldef in RTMB.h) */
extern "C"
void desolve_derivs (int *neq, double *t, double *y, double *ydot,
                     double *yout, int *ip) {
  X[0] = *t;
  for (int i=0; i < *neq; i++) X[i+1] = y[i];
  run_forward(F);
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
