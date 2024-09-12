// Pass shared pointers RTMB => Other Package
#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
void rtmb_set_shared_pointers() {
  typedef SEXP(*funptr)(SEXP);
  funptr fun = (funptr) R_FindSymbol("getSetGlobalPtr", "RTMB", NULL);
  SEXP ptr = fun(R_NilValue);
  TMBad::global_ptr = (TMBad::global**) R_ExternalPtrAddr(ptr);
}
