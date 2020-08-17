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

ignore += cachestuff/

cachestuff:
	git clone https://github.com/wzmli/MacPan_cache.git $@


######################################################################

## Structure notes

## refactor.html: refactor.Rmd
%.html: %.Rmd
	$(rmdh_r)

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

######################################################################

## Modularizing simplots

## This is no longer supposed to run on its own! Except interactively.
Sources += ../inst/params/PHAC_testify.csv
## testify_sim.Rout: testify_sim.R ../inst/params/PHAC_testify.csv; $(makeR)
## testify_simplot.Rout: testify_simplot.R testify_sim.rda; $(makeR)

Sources += slow.csv

wrap_makeR=yes

## testify_funs.R

impmakeR += sims
# testwt_N.sims.Rout: testify_sim.R testwt_N.R
%.sims.Rout: testify_sim.R %.rda testify_funs.rda sims.csv
	$(makeR)

# testwt_N.simplots.Rout: testify_simplot.R testwt_N.R
%.simplots.Rout: testify_simplot.R %.sims.rda
	$(makeR)

# testwt_N.testify_calibration.Rout: testify_calibration.R
# basic.testify_calibration.Rout: basic.R testify_calibration.R
%.testify_calibration.Rout: testify_calibration.R %.rda %.sims.rda
	$(makeR)

# basic.testcatplots.Rout: testify_testcatplot.R
%.testcatplots.Rout: testify_testcatplot.R %.sims.rda
	$(makeR)

# ode.simplots.Rout: ode.R

## Absolute testing rates (oops, that's what basic is already doing?)
## Play here anyway
# absolute.simplots.Rout: absolute.R testify_sim.R
# absolute.sims.Rout: absolute.R testify_sim.R

######################################################################

## Nested testify sims

testify_funs.Rout: testify_funs.R
	$(makeR)

batchtools.Rout: batchtools.R
	$(makeR)

## testwt_N.nested_testify.Rout: nested_testify.R testify_funs.R
%.nested_testify.Rout: nested_testify.R batchtools.rda testify_funs.rda %.rda 
	$(makeR)

## Temporarily suppressing slow dependency (bad form?)
## testwt_N.simcalib_tidy.Rout: simcalib_tidy.R
## %.simcalib_tidy.Rout: simcalib_tidy.R %.nested_testify.rda 

testwt_N.simcalib_tidy.Rout: simcalib_tidy.R
%.simcalib_tidy.Rout: simcalib_tidy.R %.rda 
	$(makeR)

testwt_N.simcalib_plot.Rout: simcalib_plot.R
%.simcalib_plot.Rout: simcalib_plot.R %.simcalib_tidy.rda
	$(makeR)

## Side pipe (faster, lighter sims)

testsim.Rout: testsim.R simple.rda slow.csv
	$(makeR)

testsim.plots.Rout: simplots.R testsim.rda
	$(makeR)

slowpars.olddiff: ../inst/params/PHAC_testify.csv slow.csv
	$(diff)

######################################################################

stoch_forecasts.Rout: ../ontario/ontario_calibration_2brks_ndt.RData stoch_forecasts.R

run_caltest.Rout: run_params_mobility_DE.Rout run_caltest.R
	$(run-R)

plot_caltest.Rout: run_caltest.RData plot_caltest.R
	$(run-R)

test0.Rout: test0.R

test_hosp.Rout: test_hosp.R

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


## Ali: simple SIR model
test_sir.Rout: test_sir.R
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
-include makestuff/yushan.mk
-include makestuff/nodes.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

