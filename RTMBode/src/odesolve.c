#include <Rinternals.h>
#include <R_ext/Rdynload.h>

typedef void * (*FUN_PTR)(void*);
static FUN_PTR run_forward;

static void* F;
static double* X;
static double* Y;
static int nstate;
static int nparms;

SEXP set_pointers(SEXP tape, SEXP x, SEXP y, SEXP nstate_, SEXP nparms_) {
  F = R_ExternalPtrAddr(tape);
  X = (double*) R_ExternalPtrAddr(x);
  Y = (double*) R_ExternalPtrAddr(y);
  nstate = INTEGER(nstate_)[0];
  nparms = INTEGER(nparms_)[0];
  return R_NilValue;
}

SEXP set_parms(SEXP x) {
  if (LENGTH(x) != nparms) Rf_error("Wrong parameter length");
  double* px = REAL(x);
  for (int i=0; i<nparms; i++) X[1 + nstate + i] = px[i];
  return R_NilValue;
}

/* Derivatives (calldef in RTMB.h) */
void desolve_derivs (int *neq, double *t, double *y, double *ydot,
                     double *yout, int *ip) {
  X[0] = *t;
  for (int i=0; i < *neq; i++) X[i+1] = y[i];
  run_forward(F);
  for (int i=0; i < *neq; i++) ydot[i] = Y[i];
}

static const R_CallMethodDef CallEntries[] = {
    {"desolve_derivs", (DL_FUNC) &desolve_derivs, 6},
    {"set_pointers",   (DL_FUNC) &set_pointers, 5},
    {"set_parms",      (DL_FUNC) &set_parms, 1},
    {NULL, NULL, 0}
};

void R_init_RTMBode(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    run_forward = (FUN_PTR) R_GetCCallable("RTMB", "ptr_forward");
}
