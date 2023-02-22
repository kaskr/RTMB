rfiles = $(basename $(wildcard *.R))
efiles = $(basename  $(basename $(wildcard *.expected.RData)) )
all_examples = $(filter $(rfiles), $(efiles))
excludes = lr_test checkConsistency
examples = $(filter-out $(excludes), $(all_examples))
outputfiles = $(examples:=.output.RData)
profiletargets = $(examples:=.profile)
logpidtargets = $(examples:=.logpid)

# Allow e.g. "make R=R-devel"
R=R

%.output.RData : %.R
	unset MAKEFLAGS; example=$(basename $<) $(R) --slave < tools/unittest.R

all: $(outputfiles)
	make report

$(all_examples) : % : %.output.RData
	

report:
	report_level=0 $(R) --slave < tools/unittest.R

report-full:
	report_level=1 $(R) --slave < tools/unittest.R

clean :
	rm -rf *.o *.so *.dll *~ core *.output.RData *.pdf *.profile *.memprofile
	rm -rf */*.o */*.so */*.dll */*~ core */*.output.RData */*.pdf */*.profile */*.memprofile

expected_output:
	cp ../../adcomp/tmb_examples/*.expected.RData .
	cp ../../adcomp/tmb_examples/validation/*.expected.RData .
