## This is McMaster pandemic:  notes subdir

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard */*.R)

pkgdesc = $(../DESCRIPTION)

######################################################################

stoch_forecasts.Rout: ../ontario/ontario_calibration_2brks_ndt.RData stoch_forecasts.R

plot_caltest.Rout: run_caltest.RData

######################################################################

## Mike's nobreak stuff

Sources += $(wildcard *.R)

run_caltest_nobreak.Rout: run_caltest_nobreak.R
	$(run-R)

## Mike's one break stuff

run_caltest_1break.Rout: run_caltest_1break.R
	$(run-R)

## Plot many sims at once
plot_caltest_nobreak.Rout: run_caltest_nobreak.Rout plot_caltest_nobreak.R
	$(run-R)

%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$<", output_format="html_document")'

%.pdf: %.Rnw
	Rscript -e 'knitr::knit2pdf("$<")'

### Makestuff

Sources += Makefile

## Sources += content.mk
## include content.mk

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	ln -s ../makestuff .
	ls $@

-include makestuff/os.mk

-include makestuff/wrapR.mk
-include makestuff/rmd.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
