## This is McMaster pandemic:  notes subdir

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard *.R) Makefile

######################################################################

## Playing with spline shapes

rshape.Rout: rshape.R
	$(makeR)

### Separating two parts

## These scripts seem to be deprecated; should we rename things and be clearer?
spline_fit.Rout: spline_fit.R
	$(makeR)

spline_projection.Rout: spline_projection.R spline_fit.rda
	$(makeR)

## Ali's second attempt to refactor
spline_projection1.Rout: spline_projection1.R
	$(makeR)

spline_projection2.Rout: spline_projection2.R
	$(makeR)

######################################################################

## Mike's development stuff
explore_testify.Rout: explore_testify.R ../inst/params/PHAC_testify.csv
	$(makeR)

## Structural pictures from Ben
## refactor.html: refactor.Rmd

## Testify sim

Sources += ../inst/params/PHAC_testify.csv
testify_sim.Rout: testify_sim.R ../inst/params/PHAC_testify.csv
	$(makeR)

testify_simplot.Rout: testify_simplot.R testify_sim.rda
	$(makeR)

## 2020 Aug 03 (Mon)
## Testing a small devtools::test() glitch that BB did not reproduce
chesterfield.Rout: chesterfield.R
	$(makeR)

######################################################################

stoch_forecasts.Rout: ../ontario/ontario_calibration_2brks_ndt.RData stoch_forecasts.R

run_caltest.Rout: run_params_mobility_DE.Rout run_caltest.R
	$(run-R)

plot_caltest.Rout: run_caltest.RData plot_caltest.R
	$(run-R)

test0.Rout: test0.R

######################################################################

## Jonathan is causing trouble

nonsense.Rout: nonsense.R
	$(makeR)

bad_nonsense.Rout: bad_nonsense.R
	$(makeR)

######################################################################

## Mike's nobreak stuff

Sources += $(wildcard *.R)

run_params_nobreak.Rout: run_params_nobreak.R
	$(run-R)

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

######################################################################
## Ali: run_sim, simple no stochasticity

ali_simple_run1.Rout: ali_simple_run1.R
	$(makeR)

######################################################################

%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$<", output_format="html_document")'

%.pdf: %.Rnw
	Rscript -e 'knitr::knit2pdf("$<")'

######################################################################

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

-include makestuff/makeR.mk
-include makestuff/rmd.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

