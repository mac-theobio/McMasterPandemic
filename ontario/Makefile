## This is McMaster pandemic: ontario subdir
## Seems outdated 2020 Aug 06 (Thu)
## Most PHAC fitting is done in PHAC_covid

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

vim_session:
	bash -cl "vmt"

######################################################################

Sources += .gitignore
Sources += $(wildcard *.R)

Ignore += $(wildcard *.rda)

Ignore += cachestuff
cachestuff:
	git clone https://github.com/wzmli/MacPan_cache.git $@

######################################################################
## Example

example.pdf: example.rmd
	$(rmdpdf)



## Ontario_calibration_setup.rda.cache
Ontario_calibration_setup.Rout: Ontario_calibration_setup.R
	$(makeR)

## Ontario_calibration_setup_nomake.R

## Ontario_PH.rda.cache:
Ontario_PH.Rout: Ontario_PH.R Ontario_calibration_setup.rda
	$(makeR)

Ontario_Mobility.Rout: Ontario_Mobility.R Ontario_calibration_setup.rda
	$(makeR)

Ontario_splines.Rout: Ontario_splines.R Ontario_calibration_setup.rda
	$(makeR)

Ontario_current.Rout: Ontario_current.R cachestuff/Ontario_calibration_setup.rda
	$(makeR)

reference_mod_forecast.Rout: reference_mod_forecast.R
	$(makeR)

#####################################################################
ontario_nbfit.Rout: ontario_clean.Rout ontario_nbfit.R

ontario_clean.Rout: ontario_clean.R
	$(run-R)

ont_keep.Rout: ontario_clean.Rout ont_keep.R
	$(run-R)

ont_cal.Rout:  ont_keep.Rout ont_cal.R
	$(run-R)

ont_cal_H.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_H.R

ont_cal_2brks.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_2brks.R

ont_cal_1brks.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_1brks.R

ont_cal_noICU.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_noICU.R

ont_cal_noICU_1brks.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_noICU_1brks.R

ont_cal_noICU_2brks.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_noICU_2brks.R

ont_cal_HD_2brks.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_HD_2brks.R

ont_cal_2brks_ndt.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_2brks_ndt.R

ont_cal_noICU_2brks_prior.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_noICU_2brks_prior.R

ontario_de_cal.Rout:  ontario_clean.Rout ont_cal.Rout ont_cal_noICU_2brks_prior.Rout ontario_de_cal.R

ont_cal_mob1.Rout: ontario_clean.Rout ont_cal_mob1.R

ont_cal_mob1_HD.Rout: ontario_clean.Rout ont_cal_mob1.Rout onb_cal_mob1_HD.R

ont_cal_spline.Rout: ontario_clean.Rout ont_cal_spline1.R

compare_calib.Rout:  ontario_clean.Rout ont_cal.RData ont_cal_noICU_1brks.RData ont_cal_noICU_2brks.RData ont_cal_noICU.RData ont_cal_2brks.RData ont_cal_1brks.RData ont_cal_H.RData ont_cal_HD_2brks.RData ont_cal_noICU_2brks_prior.RData mobility1.RData mobility1_HD.RData

ont_cal_factorial.Rout: ont_keep.Rout ont_cal_factorial.R
	$(run-R)

ont_cal_pred_Rt.Rout: ont_cal_pred_Rt.R ont_cal_factorial.Rout

ont_cal_factorial_plots.html: ont_cal_factorial_plots.Rmd ont_cal_factorial.Rout ont_cal_pred_Rt.Rout

epiestim.Rout: ontario_clean.Rout
epiestim_plot.Rout: ontario_clean.Rout ont_cal.Rout epiestim.Rout epiestim_plot.R

ontario_calibration_report.html: ontario_clean.RData ont_cal.RData epiestim.RData epiestim_plot.Rout.png ontario_nbfit.Rout.png ont_cal1.png ontario_calibration_report.Rmd

ont_cal_splinecomp_plots.Rout: ont_keep.Rout ont_cal_splinecomp.RData

Ignore += *.datatime
%.datatime: ; touch $@
ontario_clean.Rout: 0601.datatime

## NOT the current report
Sources += McMasterReport.Rnw
McMasterReport.pdf: McMasterReport.Rnw cachestuff/ont_cal_noICU_2brks_prior.rda cachestuff/ontario_clean.rda cachestuff/ont_cal_mob1.rda

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

-include makestuff/cacheflow.mk
-include makestuff/cacherepo.mk
-include makestuff/makeR.mk
-include makestuff/rmd.mk


-include makestuff/pandoc.mk
-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
