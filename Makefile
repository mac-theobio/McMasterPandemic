## This is McMaster pandemic
## https://github.com/bbolker/McMasterPandemic

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard tests/*.R)
Sources += $(wildcard man/*.Rd) NAMESPACE

######################################################################

style: misc/macpan_style.R misc/macpan_lint.R
	Rscript misc/macpan_style.R
	Rscript misc/macpan_lint.R

## The pipeR transition!

Sources += $(wildcard sandbox/*.R)

run_sim_mre.Rout: sandbox/run_sim_mre.R sandbox/ON.short.breaks.RDS
	$(run-R)

Sources += sandbox/sim.RData

jaggy.Rout: sandbox/sim.RData sandbox/jaggy.R
	$(run-R)

sandbox/kernel_test.Rout: sandbox/kernel_test.R

time_varying_mre.Rout: sandbox/time_varying_mre.R
	$(pipeR)

######################################################################

tests/moments.Rout:

break_test.Rout: notes/break_test.R
	$(run-R)

######################################################################

## Ali calibration

tests/test_calibrate.Rout: tests/test_calibrate.R

######################################################################

Sources += $(wildcard R/*.R)

Sources += dottarget.mk

######################################################################

subdirs += ontario notes

Ignore += $(subdirs)
ontario/%:
	cd ontario &&$(MAKE) Makefile
	$(makethere)

alldirs += $(subdirs)
Ignore += $(subdirs)

######################################################################

## Deprecated 2020 Apr 16 (Thu)
output/ontario_nbfit.Rout: output/ontario_nbfit.R output/ontario_clean.RData
	$(run-R)

output/ontario_cal_plots.Rout: output/ontario_cal_plots.R output/ontario_calibration.RData output/ontario_clean.RData

output/epiestim_plot.Rout: output/ontario_clean.RData output/ontario_calibration.RData output/epiestim.RData output/epiestim_plot.R

notes/ontario_calibration_report.html: output/ontario_clean.RData output/ontario_calibration.RData notes/ontario_calibration_report.Rmd output/epiestim.RData output/epiestim_plot.Rout output/ontario_cal_plots.Rout

comb_calib.Rout: notes/comb_calib.R
	$(run-R)

%.html: %.Rmd
	Rscript -e 'library("rmarkdown"); render("$<", output_format="html_document")'

%.html: %.md
	Rscript -e 'library("rmarkdown"); render("$<", output_format="html_document")'

######################################################################

# Try to break out aggfun stuff
aggfuns.Rout: notes/aggfuns.R
	$(run-R)

######################################################################

## Test before pushing 2020 Aug 03 (Mon) UNTESTED
Ignore += dtest.log
dtest:
	env SKIP_SLOW_TESTS=true Rscript -e 'devtools::test()' > $@.log

## Simple install

lpackage:
	R CMD INSTALL .

package:
	sudo R CMD INSTALL .

newpackage: pull package

######################################################################

### package-building stuff, copied from glmmTMB

R=R
# -> you can do    R=R-devel  make ....
PACKAGE=McMasterPandemic
# get VERSION from glmmTMB/DESCRIPTION
## ("::" = expand only  once, but doesn't work in make <= 3.81)
VERSION := $(shell sed -n '/^Version: /s///p' ./DESCRIPTION)
SPECVERSION := $(shell cat inst/tmb/recommended_spec_version)

testversion:
	echo "${VERSION}"
	echo "${SPECVERSION}"

Ignore += McMasterPandemic*.tar.gz
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

doc-update: R/*.R
	echo "suppressWarnings(roxygen2::roxygenize(\".\",roclets = c(\"collate\", \"rd\")))" | $(R) --slave
	@touch $@

namespace-update: R/*.R
	echo "suppressWarnings(roxygen2::roxygenize('.',roclets = 'namespace'))" | $(R) --slave
	@touch $@

pkgall: clean doc-update namespace-update install pkgcheck

pkgtest:
	echo "devtools::test('.')" | $(R) --slave

pkgcheck:
	echo "devtools::check('.')" | $(R) --slave

test-tmb:
	echo "testthat::test_file('tests/testthat/test-tmb.R')" | $(R) --slave

test-tmb-struc:
	echo "testthat::test_file('tests/testthat/test-tmb-struc.R')" | $(R) --slave

test-tmb-make-state:
	echo "testthat::test_file('tests/testthat/test-tmb-make-state.R')" | $(R) --slave

test-tmb-calibrate:
	echo "testthat::test_file('tests/testthat/test-tmb-calibrate.R')" | $(R) --slave

test-tmb-timevar:
	echo "testthat::test_file('tests/testthat/test-tmb-timevar.R')" | $(R) --slave

test-tmb-forecast:
	echo "testthat::test_file('tests/testthat/test-tmb-forecast.R')" | $(R) --slave

test-tmb-formula:
	echo "testthat::test_file('tests/testthat/test-tmb-forecast.R')" | $(R) --slave

test-tmb-simulate:
	echo "testthat::test_file('tests/testthat/test-tmb-simulate.R')" | $(R) --slave

test-tmb-all: test-tmb test-tmb-struc test-tmb-make-state test-tmb-calibrate test-tmb-timevar test-tmb-forecast test-tmb-formula test-tmb-simulate

src/McMasterPandemic.cpp: inst/tmb/**/* inst/tmb/recommended_spec_version cleanobjects
	cp inst/tmb/$(SPECVERSION)/macpan.cpp src/McMasterPandemic.cpp

cpp-sync-diff:
	diff inst/tmb/$(SPECVERSION)/macpan.cpp src/McMasterPandemic.cpp

vignettes/flex_specs.html: vignettes/flex_specs.rmd
	echo "rmarkdown::render('vignettes/flex_specs.rmd')" | $(R) --slave

vignettes/flex_intro.html: vignettes/flex_intro.rmd
	echo "rmarkdown::render('vignettes/flex_intro.rmd')" | $(R) --slave

cleanobjects:
	rm inst/tmb/*/macpan.o inst/tmb/*/macpan.so inst/tmb/*/macpan.dll || true

clean:
	find . \( -name "\.#*" -o -name "*~" -o -name ".Rhistory" \) -exec rm {} \;

CPP_SRC=

dependencies:
	Rscript misc/dependencies.R

## FIXME: depend on ??
## added $(BUILDARGS) so that this is possible:
## make install BUILDARGS="--no-build-vignettes"
build-package: $(TARBALL)
$(TARBALL): ./NAMESPACE
	$(info spec version: $(SPECVERSION))
	$(R) CMD build $(BUILDARGS) .
	mv $@ ..

## JD: Why is TARBALL pushed up a directory in the user's tree?
## Is it OK to bring it back here?
install: $(TARBALL)
	export NOT_CRAN=true; $(R) CMD INSTALL --preclean ../$<
	@touch $@

newinstall:
	make install BUILDARGS="--no-build-vignettes"

######################################################################

## Looks cool; clashes with current Bolker rules.
Ignore += maker
maker:
	git clone https://github.com/ComputationalProteomicsUnit/maker.git
## -include maker/Makefile

## Why the hell this doesn't work?
## mr_build:
mr_%:
	make $* -f maker/Makefile MAKERMAKEFILE=maker/Makefile PKGDIR=.

######################################################################

### Makestuff

Sources += Makefile

## Sources += content.mk
## include content.mk

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls $@

localstuff:
	ln -s ../makestuff .
	ls $@

-include makestuff/os.mk

-include makestuff/pipeR.mk
-include makestuff/rmd.mk

-include makestuff/git.mk
-include makestuff/visual.mk
