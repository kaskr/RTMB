install:
	R CMD INSTALL RTMB

doc-update:
	echo "library(roxygen2);suppressWarnings(roxygenize(\"RTMB\",roclets = c(\"collate\", \"rd\"), load_code=load_installed))" | R --slave
	sed -i '/RoxygenNote/d' RTMB/DESCRIPTION

unexport TEXINPUTS
pdf:
	rm -f RTMB.pdf
	R CMD Rd2pdf --no-preview RTMB

distributions:
	R --slave < distr.R

rcpp:
	echo 'Rcpp::compileAttributes("RTMB", verbose=TRUE)' | R --slave
	sed -i '/RcppEigen/d' RTMB/src/RcppExports.cpp
	sed -i '/R_CallMethodDef/ s/$$/\n    TMB_CALLDEFS,/' RTMB/src/RcppExports.cpp
	sed -i '/include.*Rcpp/ s/$$/\n#include "TMB.h"/' RTMB/src/RcppExports.cpp
	sed -i '/R_useDynamicSymbols/ s/$$/\n    TMB_CCALLABLES("RTMB");/' RTMB/src/RcppExports.cpp

test-all: linreg spatial mvrw spde sdv_multi

linreg:
	cd tmb_examples; R --slave < linreg.R

spatial:
	cd tmb_examples; R --slave < spatial.R

mvrw:
	cd tmb_examples; R --slave < mvrw.R

spde:
	cd tmb_examples; R --slave < spde.R

sdv_multi:
	cd tmb_examples; R --slave < sdv_multi.R

cran-version:
	sed -i 's/-DTMB_SAFEBOUNDS//g' RTMB/src/Makevars
	R CMD build RTMB
	git checkout RTMB/src/Makevars

cran-check:
	R CMD check --as-cran RTMB*.tar.gz
