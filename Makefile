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

######################################################################

sandbox/kernel_test.Rout: sandbox/kernel_test.R

tests/moments.Rout:

break_test.Rout: notes/break_test.R
	$(run-R)

######################################################################

%.html: %.Rmd
	Rscript -e 'library("rmarkdown"); render("$<", output_format="html_document")'

######################################################################

# Try to break out aggfun stuff
aggfuns.Rout: notes/aggfuns.R
	$(run-R)

ontario_clean.Rout: notes/ontario_clean.R
	$(run-R)

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
