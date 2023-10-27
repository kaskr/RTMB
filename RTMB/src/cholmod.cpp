#include <Matrix/stubs.c>
#include "cholmod.h"

extern "C" {
  CHM_DN cholmod_solve(int sys, CHM_FR L, CHM_DN B, CHM_CM Common) {
    return M_cholmod_solve(sys, L, B, Common);
  }
  CHM_FR cholmod_analyze(CHM_SP A, CHM_CM Common) {
    return M_cholmod_analyze(A, Common);
  }
  int cholmod_free_dense(CHM_DN *A, CHM_CM Common) {
    return M_cholmod_free_dense(A, Common);
  }
  int cholmod_finish(CHM_CM Common) {
    return M_cholmod_finish(Common);
  }
  int cholmod_start(CHM_CM Common) {
    return M_cholmod_start(Common);
  }
  int cholmod_factorize_p(CHM_SP A, double beta[2], int *fset,
                          size_t fsize, CHM_FR L, CHM_CM Common) {
    Common -> try_catch = TRUE; // Disable 'warning: not positive definite'
    Common -> quick_return_if_not_posdef = TRUE;
    return M_cholmod_factorize_p(A, beta, fset, fsize, L, Common);
  }
  int cholmod_free_factor(CHM_FR *L, CHM_CM Common) {
    return M_cholmod_free_factor(L, Common);
  }
}
