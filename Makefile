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

######################################################################

package:
	sudo R CMD INSTALL .

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

## -include makestuff/wrapR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

### package-building stuff, copied from glmmTMB

R=R
# -> you can do    R=R-devel  make ....
PACKAGE=glmmTMB
# get VERSION from glmmTMB/DESCRIPTION  
## ("::" = expand only  once, but doesn't work in make <= 3.81)
VERSION := $(shell sed -n '/^Version: /s///p' DESCRIPTION)

testversion:
	echo "${VERSION}"

TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

doc-update: R/*.R
	echo "suppressWarnings(roxygen2::roxygenize(\".\",roclets = c(\"collate\", \"rd\")))" | $(R) --slave
	@touch doc-update

pkgtest:
	echo "devtools::test('.')" | $(R) --slave

pkgcheck:
	echo "devtools::check('.')" | $(R) --slave
