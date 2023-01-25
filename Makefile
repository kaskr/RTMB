install:
	R CMD INSTALL RTMB

rcpp:
	echo 'Rcpp::compileAttributes("RTMB", verbose=TRUE)' | R --slave
	sed -i '/RcppEigen/d' RTMB/src/RcppExports.cpp
	sed -i '/R_CallMethodDef/ s/$$/\n    TMB_CALLDEFS,/' RTMB/src/RcppExports.cpp
	sed -i '/include.*Rcpp/ s/$$/\n#include "TMB.h"/' RTMB/src/RcppExports.cpp
	sed -i '/R_useDynamicSymbols/ s/$$/\n    TMB_CCALLABLES("RTMB");/' RTMB/src/RcppExports.cpp

test-all: linreg spatial mvrw spde

linreg:
	cd tmb_examples; R --slave < linreg.R

spatial:
	cd tmb_examples; R --slave < spatial.R

mvrw:
	cd tmb_examples; R --slave < mvrw.R

spde:
	cd tmb_examples; R --slave < spde.R
