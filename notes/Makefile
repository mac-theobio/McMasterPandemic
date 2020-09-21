## This is McMaster pandemic:  notes subdir

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

vim_session:
	bash -cl "vmt"

wrap_makeR=yes

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

## Create a phenomenological curve and
## compare two different spline fits directly to the curve
## could be modularized
spline_fit.Rout: spline_fit.R
	$(makeR)

## Takes parameters from the ns fit and feeds it to MacPan to make data
Sources += spline.csv
spline_sim.Rout: spline_sim.R spline_fit.rda spline.csv
	$(makeR)

spline_calib.Rout: spline_calib.R spline_sim.rda spline.csv
	$(makeR)

spline_compare.Rout: spline_compare.R cachestuff/spline_calib.rda spline_sim.rda
	$(makeR)

spline_recalib.Rout: spline_recalib.R batchtools.rda cachestuff/spline_calib.rda 
	$(makeR)

collect_recalib.Rout: collect_recalib.R 
	$(makeR)

spline_recalib_compare.Rout: spline_recalib_compare.R cachestuff/spline_recalib.rda cachestuff/spline_calib.rda spline_sim.rda

######################################################################

## Older attempts from Ali

spline_projection.Rout: spline_projection.R spline_fit.rda
	$(makeR)

## Ali's second attempt to refactor
spline_projection1.Rout: spline_projection1.R
	$(makeR)

spline_projection2.Rout: spline_projection2.R
	$(makeR)

## parameter synthesis project

midas_estimates.Rout: midas_estimates.R
	$(makeR)

######################################################################

## Structural pictures from Ben
## refactor.html: refactor.Rmd

######################################################################

## Modularizing simplots

## This is no longer supposed to run on its own! Except interactively.
Sources += ../inst/params/PHAC_testify.csv
Sources += slow.csv

## How many parameter files do we have?
## testwt_N.R basic.R simple.R ode.R

## testify_funs.R

impmakeR += sims
# testwt_N.sims.Rout: testify_sim.R testwt_N.R
%.sims.Rout: testify_sim.R %.rda testify_funs.rda
	$(makeR)

# testwt_N.simplots.Rout: testify_simplot.R testwt_N.R
%.simplots.Rout: testify_simplot.R %.sims.rda
	$(makeR)


## BMB's new code

# testwt_N.testify_simreportplot.Rout: 
%.testify_simreportplot.Rout: testify_simreportplot.R %.sims.rda
	$(makeR)


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

ali_test_sir.html: ali_test_sir.Rmd

## Side pipe (faster, lighter sims)
testsim.Rout: testsim.R simple.rda slow.csv
	$(makeR)

testsim.plots.Rout: simplots.R testsim.rda
	$(makeR)

slowpars.olddiff: ../inst/params/PHAC_testify.csv slow.csv
	$(diff)

######################################################################

## Old stuff, do we still want? keep for now, mostly exploration and debugging!!! 

stoch_forecasts.Rout: ../ontario/ontario_calibration_2brks_ndt.RData stoch_forecasts.R

run_caltest.Rout: run_params_mobility_DE.Rout run_caltest.R
	$(run-R)

plot_caltest.Rout: run_caltest.RData plot_caltest.R
	$(run-R)

test_hosp.Rout: test_hosp.R

######################################################################

## Old simulation exploration pipeline!

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

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	ln -s ../makestuff .
	ls $@

-include makestuff/os.mk

-include makestuff/cacheflow.mk
-include makestuff/cacherepo.mk

-include makestuff/makeR.mk
-include makestuff/rmd.mk
-include makestuff/yushan.mk
-include makestuff/nodes.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

