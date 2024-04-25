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
	sed -i '/inst.include.RTMB/d' RTMB/src/RcppExports.cpp
	sed -i '/R_CallMethodDef/ s/$$/\n    TMB_CALLDEFS,/' RTMB/src/RcppExports.cpp
	sed -i '/include.*Rcpp/ s/$$/\n#include "RTMB.h"/' RTMB/src/RcppExports.cpp
	sed -i '/R_useDynamicSymbols/ s/$$/\n    RTMB_CCALLABLES;/' RTMB/src/RcppExports.cpp

rtmb-stubs:
	cat RTMB/src/TMB.cpp > RTMB/inst/include/RTMB_stubs.cpp
	cat RTMB/src/RTMB.cpp >> RTMB/inst/include/RTMB_stubs.cpp
	echo '#include "rtmb_set_shared_pointers.cpp"' >> RTMB/inst/include/RTMB_stubs.cpp
	echo '#define TMB_SKINNY' > RTMB/inst/include/TMB.h
	cat RTMB/src/TMB.h >> RTMB/inst/include/TMB.h
	cp RTMB/src/config.h RTMB/inst/include
	cp RTMB/src/RTMB.h RTMB/inst/include

test:
	R -s -e "tinytest::test_package('RTMB')"

test-all: linreg spatial mvrw spde sdv_multi test-sparse

test-sparse:
	R --slave < tests/test-sparse.R

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
	R CMD build RTMB

cran-check:
	R CMD check --as-cran RTMB*.tar.gz

parallel-version:
	rm -rf RTMBp
	cp -r RTMB RTMBp
	sed -i 's/RTMB::/RTMBp::/g' RTMBp/R/*.R
	sed -i 's/"RTMB"/"RTMBp"/g' RTMBp/R/*.R
	sed -i 's/library(RTMB)/library(RTMBp);TMB::openmp(1,DLL="RTMBp")/g' RTMBp/vignettes/*.rmd
	sed -i 's/RTMB/RTMBp/g' RTMBp/inst/tinytest/*.R
	sed -i 's/RTMB/RTMBp/g' RTMBp/DESCRIPTION
	sed -i 's/RTMB/RTMBp/g' RTMBp/NAMESPACE
	sed -i 's/"RTMB"/"RTMBp"/g' RTMBp/src/RTMB.h
	sed -i 's/R_init_RTMB/R_init_RTMBp/g' RTMBp/src/RcppExports.cpp
	echo 'SystemRequirements: GNU make' >> RTMBp/DESCRIPTION
	echo 'PKG_LIBS+=$$(SHLIB_OPENMP_CXXFLAGS)' >> RTMBp/src/Makevars
	echo 'PKG_CXXFLAGS+=$$(SHLIB_OPENMP_CXXFLAGS)' >> RTMBp/src/Makevars
	echo 'PKG_LIBS+=$$(SHLIB_OPENMP_CXXFLAGS)' >> RTMBp/src/Makevars.win
	echo 'PKG_CXXFLAGS+=$$(SHLIB_OPENMP_CXXFLAGS)' >> RTMBp/src/Makevars.win
	echo '.onLoad <- function(libname, pkgname) { TMB::openmp(parallel::detectCores(), autopar=TRUE, DLL="RTMBp") }' >> RTMBp/R/zzz.R
	git branch -D RTMBp
	git checkout -b RTMBp
	git add `git ls-files RTMB | sed 's/^RTMB/RTMBp/g'`
	git commit -m "autogenerate"
	git checkout master

## vignettes
%.html: %.rmd
	cd RTMB/vignettes; echo "rmarkdown::render(basename(\"$<\"))" | R --slave
vignettes-build: RTMB/vignettes/RTMB-introduction.html RTMB/vignettes/RTMB-advanced.html

rcpp-rtmbXtra:
	echo 'Rcpp::compileAttributes("rtmbXtra", verbose=TRUE)' | R --slave
	sed -i '/RcppExport void R_init/ s/^/void rtmb_set_shared_pointers();\n/' rtmbXtra/src/RcppExports.cpp
	sed -i '/R_useDynamicSymbols/ s/$$/\n    rtmb_set_shared_pointers();/' rtmbXtra/src/RcppExports.cpp
