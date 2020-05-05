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

run_caltest.Rout: run_params_breaks_DE.Rout run_caltest.R
	$(run-R)

plot_caltest.Rout: run_caltest.RData plot_caltest.R
	$(run-R)

test0.Rout: test0.R

######################################################################

## Mike's nobreak stuff

Sources += $(wildcard *.R)

run_caltest_nobreak.Rout: run_params_nobreak.Rout run_caltest.R
	$(run-R)

run_params_breaks.Rout: run_params_breaks.R
	$(run-R)

run_caltest_breaks.Rout: run_params_breaks.Rout run_caltest.R
	$(run-R)

run_DEoptim_breaks.Rout: run_params_breaks.Rout run_DEoptim.R
	$(run-R)

## Plot many sims at once
plot_caltest_nobreak.Rout: run_caltest_nobreak.Rout plot_caltest.R
	$(run-R)

explore_nobreaks.Rout: run_caltest_nobreak.Rout explore_nobreaks.R
	$(run-R)

plot_caltest_breaks.Rout: run_caltest_breaks.Rout plot_caltest.R
	$(run-R)

plot_DEoptim_breaks.Rout: run_DEoptim_breaks.Rout plot_caltest.R
	$(run-R)


alt_cal.Rout:  ../ontario/ontario_calibration.RData ../ontario/ontario_calibration_noICU_2brks_prior.RData alt_cal.R

explore_ontcal.Rout: alt_cal.RData ../ontario/ontario_de_cal.RData ../ontario/ontario_calibration_noICU_2brks_prior.RData explore_ontcal.R


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
