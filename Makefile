## This is McMaster pandemic

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard */*.R)

Sources += $(wildcard man/*.Rd) NAMESPACE

Sources += dottarget.mk

######################################################################

Sources += $(wildcard sandbox/*.R) sandbox/sim.RData

jaggy.Rout: sandbox/sim.RData sandbox/jaggy.R
	$(run-R)

sandbox/kernel_test.Rout: sandbox/kernel_test.R

tests/moments.Rout:

break_test.Rout: notes/break_test.R
	$(run-R)

######################################################################

subdirs += ontario notes

ontario/%:
	cd ontario &&$(MAKE) Makefile
	$(makethere)

alldirs += $(subdirs)

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

######################################################################

# Try to break out aggfun stuff
aggfuns.Rout: notes/aggfuns.R
	$(run-R)

lpackage:
	R CMD INSTALL .

package:
	sudo R CMD INSTALL .

######################################################################

### package-building stuff, copied from glmmTMB

R=R
# -> you can do    R=R-devel  make ....
PACKAGE=McMasterPandemic
# get VERSION from glmmTMB/DESCRIPTION  
## ("::" = expand only  once, but doesn't work in make <= 3.81)
VERSION := $(shell sed -n '/^Version: /s///p' ./DESCRIPTION)

testversion:
	echo "${VERSION}"

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

clean:
	find . \( -name "\.#*" -o -name "*~" -o -name ".Rhistory" \) -exec rm {} \;

CPP_SRC=

## FIXME: depend on ??
build-package: $(TARBALL)
$(TARBALL): ./NAMESPACE
	$(R) CMD build --resave-data=no .
	mv $@ ..

install: $(TARBALL)
	export NOT_CRAN=true; $(R) CMD INSTALL --preclean ../$<
	@touch $@

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

-include makestuff/wrapR.mk
-include makestuff/rmd.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
