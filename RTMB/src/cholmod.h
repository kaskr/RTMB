// From the Matrix package
#include <Matrix/cholmod.h>

// Eigen's CHOLMOD support uses these
extern "C" {
  int cholmod_start(CHM_CM);
  int cholmod_l_start(CHM_CM);
  int cholmod_finish(CHM_CM);
  int cholmod_l_finish(CHM_CM);
  int cholmod_free_factor(CHM_FR *L, CHM_CM Common);
  int cholmod_l_free_factor(CHM_FR *L, CHM_CM Common);
  int cholmod_free_dense(CHM_DN *A, CHM_CM Common);
  int cholmod_l_free_dense(CHM_DN *A, CHM_CM Common);
  int cholmod_free_sparse(CHM_SP *A, CHM_CM Common);
  int cholmod_l_free_sparse(CHM_SP *A, CHM_CM Common);
  CHM_FR cholmod_analyze(CHM_SP A, CHM_CM Common);
  CHM_FR cholmod_l_analyze(CHM_SP A, CHM_CM Common);
  CHM_DN cholmod_solve(int sys, CHM_FR L, CHM_DN B, CHM_CM Common);
  CHM_DN cholmod_l_solve(int sys, CHM_FR L, CHM_DN B, CHM_CM Common);
  CHM_SP cholmod_spsolve(int sys, CHM_FR L, CHM_SP B, CHM_CM Common);
  CHM_SP cholmod_l_spsolve(int sys, CHM_FR L, CHM_SP B, CHM_CM Common);
  int cholmod_factorize_p(CHM_SP A, double beta[2], int *fset,
                          size_t fsize, CHM_FR L, CHM_CM Common);
  int cholmod_l_factorize_p(CHM_SP A, double beta[2], SuiteSparse_long *fset,
                            size_t fsize, CHM_FR L, CHM_CM Common);
  // FIXME
  void cholmod_print_common(const char*, CHM_CM Common);
}
