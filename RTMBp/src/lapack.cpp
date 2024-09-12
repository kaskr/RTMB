//#include <Rcpp.h>
#include "RTMB.h"
#define USE_FC_LEN_T
#include <R.h>
#include <Rversion.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif

/* sytrisol: Given lower L and lower W, find symmetric S such that

   TRIL(S*L)=TRIL(W)
*/
static void sytrisol_recursion(double* L, double* W, int n, int isub, int nsub) {
  double* S = W;
  if (nsub == 1) {
    int k = n * isub + isub;
    S[k] = W[k] / L[k];
    return;
  }
  // Split current block in left <= right
  int nsub1  = nsub / 2; // floor of division
  int nsub2 = nsub - nsub1;
  sytrisol_recursion(L, W, n, isub + nsub1, nsub2);
  double ONE=1.0, MONE=-1.0;
#define BLK11(X) (X + n*isub + isub)
#define BLK21(X) (X + n*isub + isub + nsub1)
#define BLK22(X) (X + n*(isub+nsub1) + isub + nsub1)
  F77_CALL(dsymm)("L", "L",
                  &nsub2, &nsub1, &MONE, BLK22(S), &n, BLK21(L), &n, &ONE, BLK21(W), &n
                  FCONE FCONE);
  F77_CALL(dtrsm)("R", "L", "N", "N",
                  &nsub2, &nsub1, &ONE, BLK11(L), &n, BLK21(S), &n
                  FCONE FCONE FCONE FCONE);
  F77_CALL(dgemm)("T", "N",
                  &nsub1, &nsub1, &nsub2, &MONE, BLK21(S), &n, BLK21(L), &n, &ONE, BLK11(W), &n
                  FCONE FCONE);
  sytrisol_recursion(L, W, n, isub, nsub1);
#undef BLK11
#undef BLK21
#undef BLK22
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sytrisol(Rcpp::NumericMatrix L, Rcpp::NumericMatrix W) {
  int n = L.rows();
  double* Lptr = &(L[0]);
  Rcpp::NumericMatrix W_(Rcpp::clone(W));
  double* Wptr = &(W_[0]);
  sytrisol_recursion(Lptr, Wptr, n, 0, n);
  // Zero upper triangle
  for(int j=0; j<n; j++) {
    for(int i=0; i<j; i++) {
      Wptr[i+j*n] = 0;
    }
  }
  return W_;
}
